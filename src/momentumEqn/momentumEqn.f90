! ************************************************************************************** !
!    TBFsolver - DNS turbulent bubbly flow solver
!    Copyright (C) 2018  University of Twente.
!
!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with this program.  If not, see <https://www.gnu.org/licenses/>.
! ************************************************************************************** !

module momentumEqnMod
	
	use timeMod
	
	implicit none
	
	!storage old fluxes for AB2 scheme
	type(vfield) :: gPhi0, phi0


	type, public :: momentumEqn
	
	
		!keep a pointer to grid
		type(grid), pointer :: ptrMesh_ => NULL() 
	
		real(DP), allocatable, dimension(:,:,:) :: phiX_, phiY_, phiZ_
		real(DP), allocatable, dimension(:,:,:) :: phiPrevX_, phiPrevY_, phiPrevZ_
		
		!volume stag cells ux
		real(DP), allocatable, dimension(:,:,:) :: Vsx_
		
		!momentum equation bounds
		integer :: isx_, iex_, jsx_, jex_, ksx_, kex_
		integer :: isy_, iey_, jsy_, jey_, ksy_, key_
		integer :: isz_, iez_, jsz_, jez_, ksz_, kez_
		
		!convection scheme parFile
		integer :: scheme_
		
		!keep a pointer to time
		type(time), pointer :: ptrTime_ => NULL() 
		
		!uniform source term 
		real(DP) :: fs_, g_, gCH_
		
		!flow rate
		real(DP) :: Q0_
		
		!theoretical shear Re
		real(DP) :: Ret_
		
		!flow control
		integer :: flowCtrl_

	end type
	
	private :: updateConveDiff
	private :: compute_limiters
	private :: compute_CD_limiter
	private :: compute_UD_limiter
	private :: compute_QUICK_limiter
	private :: compute_VanLeer_limiter
	private :: compute_Superbee_limiter
	private :: compute_localVF_limiter
	private :: addPressureGrad
	private :: addSource
	private :: addConvDiff
	private :: setEqnBounds
	private :: initFluxes
	private :: resetFluxes
	private :: copyOldFluxes
	private :: info
	private :: initOldTimeFlux
#ifdef FAST_MODE	
	private :: computeOldPressGrad
#endif
	
	public :: momentumEqnCTOR
	public :: makeVelocityDivFree
	public :: solveMomentumEqn
	public :: setFlowRate
	public :: printOldTimeFlux

  	
contains


!========================================================================================!
	subroutine momentumEqnCTOR(this,gMesh,mesh,gu,u,rt)
		type(momentumEqn) :: this
		type(grid), intent(in), target :: gMesh,mesh
        type(vfield), intent(in) :: gu,u
        type(time), intent(in), target :: rt
        type(parFile) :: pfile_conv, pfile_flow, pfile_g
        logical :: read_fr
        
        call parFileCTOR(pfile_conv,'schemes','specs')
        !read convection scheme
        call readParameter(pfile_conv,this%scheme_,'convection_scheme')
        
        call parFileCTOR(pfile_flow,'flowControl','specs')
        call readParameter(pfile_flow,this%flowCtrl_,'flowCtrl')

		this%ptrMesh_ => mesh
		this%ptrTime_ => rt
        
        call setEqnBounds(this,u)
        
        !allocate fields
        !ux
        call allocateArray(this%phiX_,this%isx_,this%iex_,this%jsx_,this%jex_,this%ksx_,this%kex_)
        call allocateArray(this%phiPrevX_,this%isx_,this%iex_,this%jsx_,this%jex_,this%ksx_,this%kex_)
        !uy
        call allocateArray(this%phiY_,this%isy_,this%iey_,this%jsy_,this%jey_,this%ksy_,this%key_)
        call allocateArray(this%phiPrevY_,this%isy_,this%iey_,this%jsy_,this%jey_,this%ksy_,this%key_)
        !uz
        call allocateArray(this%phiZ_,this%isz_,this%iez_,this%jsz_,this%jez_,this%ksz_,this%kez_)
        call allocateArray(this%phiPrevZ_,this%isz_,this%iez_,this%jsz_,this%jez_,this%ksz_,this%kez_)
        
        !init flux to zero
        call initFluxes(this)
        
        !init source term
        this%fs_ = 0.d0
        
        if (this%flowCtrl_==1) then
        	call readParameter(pfile_flow,this%Ret_,'Ret')
        end if
        
        if (this%flowCtrl_==2) then
        	call readParameter(pfile_flow,read_fr,'read_flow_rate')
        	if (read_fr) then
        		call readParameter(pfile_flow,this%Q0_,'flow_rate')
        		this%Q0_=this%Q0_*mesh%Lxg_*mesh%Lzg_
        	else
        		call computeFlowRate(u,this%Q0_)
        	end if
        end if
        
        !read gravity
		call parFileCTOR(pfile_g,'parameters','specs')
        call readParameter(pfile_g,this%g_,'g')
        call readParameter(pfile_g,this%gCH_,'gCH')
        
        !init u0 if AB2
        call initOldTimeFlux(this,gMesh,mesh,gu,rt)
					    
	end subroutine
!========================================================================================!

!========================================================================================!
    subroutine solveMomentumEqn(this,u,p,mu,rho,st,c)
    	type(momentumEqn), intent(inout) :: this
    	type(vfield), intent(inout) :: u
    	type(field), intent(in) :: p, mu, rho, c
    	type(vfield), intent(in) :: st
    	real(DP) :: start, finish
    	

    	start = MPI_Wtime() 
    	
    	!copy fluxes to Prev arrays
		call copyOldFluxes(this)
			
    	!reset fluxes for computation of new fluxes
    	call resetFluxes(this)
    	
    	!update fluxes
    	call updateConveDiff(this,u,c,mu,rho)
		call addConvDiff(this,u)
    	
    	!add sources
    	call addSource(this,u,rho,st)
    	
    	!update boundaries
    	call updateBoundariesV(u)
    	
    	finish = MPI_Wtime()
    	
    	!print-out info
    	call info(this,finish-start)
    	
        
    end subroutine
!========================================================================================!

!========================================================================================!
    subroutine updateConveDiff(this,u,c,mu,rho)
    	type(momentumEqn), intent(inout) :: this
    	type(vfield), intent(in) :: u
    	type(field), intent(in) :: c, mu, rho
    	integer :: im, imm, jm, jmm, km, kmm
    	integer :: ip, ipp, jp, jpp, kp, kpp
    	real(DP) :: qp, Ap, Bp, Fp
    	real(DP) :: qm, Am, Bm, Fm
    	real(DP) :: mur, mul, dxx, dyy, dzz, dxxt, dyyt, dzzt
    	real(DP) :: psi_hpup, psi_hpum, psi_hmup, psi_hmum
    	integer :: i, j, k, scheme
    	real(DP) :: invrho

    	
		!select scheme
		scheme=this%scheme_


		! x comp
		
 		!$OMP PARALLEL DO DEFAULT(none) &
		!$OMP SHARED(this,u,c,mu,rho,scheme) &
		!$OMP PRIVATE(i,j,k) &
		!$OMP PRIVATE(ip,im,ipp,imm) &
		!$OMP PRIVATE(jp,jm,jpp,jmm) &
		!$OMP PRIVATE(kp,km,kpp,kmm) &
		!$OMP PRIVATE(qp,Ap,Bp,Fp) &
		!$OMP PRIVATE(qm,Am,Bm,Fm) &
		!$OMP PRIVATE(dxx,dyy,dzz,dxxt,dyyt,dzzt) &
		!$OMP PRIVATE(psi_hpup,psi_hpum,psi_hmup,psi_hmum) &
		!$OMP PRIVATE(mur,mul,invrho)
		do k = this%ksx_,this%kex_
		
			kp = k + 1
			km = k - 1
			kpp = kp + 1
			kmm = km - 1
			
				do j = this%jsx_,this%jex_
			
					jp = j + 1
					jm = j - 1
					jpp = jp + 1
					jmm = jm - 1
				
						do i = this%isx_,this%iex_
				
							ip = i + 1
							im = i - 1
							ipp = ip + 1
							imm = im - 1
		
							!***********************  convection  ***********************!
							!d (uu) / dx
							call compute_limiters(u%ux_%f_(i,j,k),u%ux_%f_(im,j,k),u%ux_%f_(ip,j,k),&
											      u%ux_%f_(ipp,j,k),u%ux_%f_(imm,j,k),&
											      c%f_(ip,j,k),c%f_(i,j,k),c%f_(i,jp,k),c%f_(i,jm,k),&
											      c%f_(i,j,kp),c%f_(i,j,km),&
											      psi_hpup,psi_hpum,psi_hmup,psi_hmum,scheme)

							qp = 0.5d0*(u%ux_%f_(i,j,k)+u%ux_%f_(ip,j,k))
							Ap = u%ux_%f_(i,j,k)+0.5d0*psi_hpup*(u%ux_%f_(ip,j,k)-u%ux_%f_(i,j,k))
							Bp = u%ux_%f_(ip,j,k)+0.5d0*psi_hpum*(u%ux_%f_(i,j,k)-u%ux_%f_(ip,j,k))
							Fp = max(qp,0.d0)*Ap + min(qp,0.d0)*Bp
					
							qm = 0.5d0*(u%ux_%f_(i,j,k)+u%ux_%f_(im,j,k))
							Am = u%ux_%f_(im,j,k)+0.5d0*psi_hmup*(u%ux_%f_(i,j,k)-u%ux_%f_(im,j,k))
							Bm = u%ux_%f_(i,j,k)+0.5d0*psi_hmum*(u%ux_%f_(im,j,k)-u%ux_%f_(i,j,k))
							Fm = max(qm,0.d0)*Am + min(qm,0.d0)*Bm
					
							this%phiX_(i,j,k) = this%phiX_(i,j,k) - (Fp-Fm)/this%ptrMesh_%dxc_(ip)
					
							!d (vu) / dy
							call compute_limiters(u%ux_%f_(i,j,k),u%ux_%f_(i,jm,k),u%ux_%f_(i,jp,k),&
											      u%ux_%f_(i,jpp,k),u%ux_%f_(i,jmm,k),&
											      c%f_(ip,j,k),c%f_(i,j,k),c%f_(i,jp,k),c%f_(i,jm,k),&
											      c%f_(i,j,kp),c%f_(i,j,km),&
											      psi_hpup,psi_hpum,psi_hmup,psi_hmum,scheme)

							qp = 0.5d0*(u%uy_%f_(ip,j,k)+u%uy_%f_(i,j,k))
							Ap = u%ux_%f_(i,j,k)+0.5d0*psi_hpup*(u%ux_%f_(i,jp,k)-u%ux_%f_(i,j,k))
							Bp = u%ux_%f_(i,jp,k)+0.5d0*psi_hpum*(u%ux_%f_(i,j,k)-u%ux_%f_(i,jp,k))
							Fp = max(qp,0.d0)*Ap + min(qp,0.d0)*Bp
					
							qm = 0.5d0*(u%uy_%f_(ip,jm,k)+u%uy_%f_(i,jm,k))
							Am = u%ux_%f_(i,jm,k)+0.5d0*psi_hmup*(u%ux_%f_(i,j,k)-u%ux_%f_(i,jm,k))
							Bm = u%ux_%f_(i,j,k)+0.5d0*psi_hmum*(u%ux_%f_(i,jm,k)-u%ux_%f_(i,j,k))
							Fm = max(qm,0.d0)*Am + min(qm,0.d0)*Bm
					
							this%phiX_(i,j,k) = this%phiX_(i,j,k) - (Fp-Fm)/this%ptrMesh_%dyf_(j)
					
							!d (wu) / dz
							call compute_limiters(u%ux_%f_(i,j,k),u%ux_%f_(i,j,km),u%ux_%f_(i,j,kp),&
											      u%ux_%f_(i,j,kpp),u%ux_%f_(i,j,kmm),&
											      c%f_(ip,j,k),c%f_(i,j,k),c%f_(i,jp,k),c%f_(i,jm,k),&
											      c%f_(i,j,kp),c%f_(i,j,km),&
											      psi_hpup,psi_hpum,psi_hmup,psi_hmum,scheme)							
							
							qp = 0.5d0*(u%uz_%f_(i,j,k)+u%uz_%f_(ip,j,k))
							Ap = u%ux_%f_(i,j,k)+0.5d0*psi_hpup*(u%ux_%f_(i,j,kp)-u%ux_%f_(i,j,k))
							Bp = u%ux_%f_(i,j,kp)+0.5d0*psi_hpum*(u%ux_%f_(i,j,k)-u%ux_%f_(i,j,kp))
							Fp = max(qp,0.d0)*Ap + min(qp,0.d0)*Bp
					
							qm = 0.5d0*(u%uz_%f_(i,j,km)+u%uz_%f_(ip,j,km))
							Am = u%ux_%f_(i,j,km)+0.5d0*psi_hmup*(u%ux_%f_(i,j,k)-u%ux_%f_(i,j,km))
							Bm = u%ux_%f_(i,j,k)+0.5d0*psi_hmum*(u%ux_%f_(i,j,km)-u%ux_%f_(i,j,k))
							Fm = max(qm,0.d0)*Am + min(qm,0.d0)*Bm
					
							this%phiX_(i,j,k) = this%phiX_(i,j,k) - (Fp-Fm)/this%ptrMesh_%dzf_(k)
							
							!***********************  diffusion  ***********************!
							mur = mu%f_(ip,j,k)
							mul = mu%f_(i,j,k)
							dxx = (  mur*(u%ux_%f_(ip,j,k)-u%ux_%f_(i,j,k))/(this%ptrMesh_%dxf_(ip))		&
								   - mul*(u%ux_%f_(i,j,k)-u%ux_%f_(im,j,k))/(this%ptrMesh_%dxf_(i))		&
								  )/(this%ptrMesh_%dxc_(ip))
							dxxt = dxx
							
							mur = 0.25d0*(mu%f_(ip,j,k)+mu%f_(i,j,k)+mu%f_(ip,jp,k)+mu%f_(i,jp,k))
							mul = 0.25d0*(mu%f_(ip,j,k)+mu%f_(i,j,k)+mu%f_(ip,jm,k)+mu%f_(i,jm,k)) 	  
							dyy = (  mur*(u%ux_%f_(i,jp,k)-u%ux_%f_(i,j,k))/(this%ptrMesh_%dyc_(jp))		&
								   - mul*(u%ux_%f_(i,j,k)-u%ux_%f_(i,jm,k))/(this%ptrMesh_%dyc_(j))		&
								  )/(this%ptrMesh_%dyf_(j))
							dyyt = (  mur*(u%uy_%f_(ip,j,k)-u%uy_%f_(i,j,k))/(this%ptrMesh_%dxc_(ip))		&
								    - mul*(u%uy_%f_(ip,jm,k)-u%uy_%f_(i,jm,k))/(this%ptrMesh_%dxc_(ip))		&
								   )/(this%ptrMesh_%dyf_(j))
						
							mur = 0.25d0*(mu%f_(ip,j,k)+mu%f_(i,j,k)+mu%f_(ip,j,kp)+mu%f_(i,j,kp)) 
							mul = 0.25d0*(mu%f_(ip,j,k)+mu%f_(i,j,k)+mu%f_(ip,j,km)+mu%f_(i,j,km))	  
							dzz = (  mur*(u%ux_%f_(i,j,kp)-u%ux_%f_(i,j,k))/(this%ptrMesh_%dzc_(kp))		&
								   - mul*(u%ux_%f_(i,j,k)-u%ux_%f_(i,j,km))/(this%ptrMesh_%dzc_(k))		&
								  )/(this%ptrMesh_%dzf_(k))
							dzzt = (  mur*(u%uz_%f_(ip,j,k)-u%uz_%f_(i,j,k))/(this%ptrMesh_%dxc_(ip))		&
								    - mul*(u%uz_%f_(ip,j,km)-u%uz_%f_(i,j,km))/(this%ptrMesh_%dxc_(ip))		&
								   )/(this%ptrMesh_%dzf_(k))
							
							
							invrho = 0.5d0*(1.d0/rho%f_(ip,j,k)+1.d0/rho%f_(i,j,k))
							this%phiX_(i,j,k) = this%phiX_(i,j,k) + (dxx + dyy + dzz +	&
											   						 dxxt + dyyt + dzzt)*invrho	
					
						end do
				end do
		end do
		!$OMP END PARALLEL DO 
		
			
		! y comp
		
 		!$OMP PARALLEL DO DEFAULT(none) &
		!$OMP SHARED(this,u,c,mu,rho,scheme) &
		!$OMP PRIVATE(i,j,k) &
		!$OMP PRIVATE(ip,im,ipp,imm) &
		!$OMP PRIVATE(jp,jm,jpp,jmm) &
		!$OMP PRIVATE(kp,km,kpp,kmm) &
		!$OMP PRIVATE(qp,Ap,Bp,Fp) &
		!$OMP PRIVATE(qm,Am,Bm,Fm) &
		!$OMP PRIVATE(dxx,dyy,dzz,dxxt,dyyt,dzzt) &
		!$OMP PRIVATE(psi_hpup,psi_hpum,psi_hmup,psi_hmum) &
		!$OMP PRIVATE(mur,mul,invrho)
		do k = this%ksy_,this%key_
		
			kp = k + 1
			km = k - 1
			kpp = kp + 1
			kmm = km - 1
			
				do j = this%jsy_,this%jey_
			
					jp = j + 1
					jm = j - 1
					jpp = jp + 1
					jmm = jm - 1
				
						do i = this%isy_,this%iey_
				
							ip = i + 1
							im = i - 1
							ipp = ip + 1
							imm = im - 1
					
							!***********************  convection  ***********************!
							!d (uv) / dx
							call compute_limiters(u%uy_%f_(i,j,k),u%uy_%f_(im,j,k),u%uy_%f_(ip,j,k),&
											      u%uy_%f_(ipp,j,k),u%uy_%f_(imm,j,k),&
											      c%f_(ip,j,k),c%f_(im,j,k),c%f_(i,jp,k),c%f_(i,j,k),&
											      c%f_(i,j,kp),c%f_(i,j,km),&
											      psi_hpup,psi_hpum,psi_hmup,psi_hmum,scheme)

							qp = 0.5d0*(u%ux_%f_(i,jp,k)*this%ptrMesh_%dyf_(jp)+ &
							            u%ux_%f_(i,j,k)*this%ptrMesh_%dyf_(j))/this%ptrMesh_%dyc_(jp)
							Ap = u%uy_%f_(i,j,k)+0.5d0*psi_hpup*(u%uy_%f_(ip,j,k)-u%uy_%f_(i,j,k))
							Bp = u%uy_%f_(ip,j,k)+0.5d0*psi_hpum*(u%uy_%f_(i,j,k)-u%uy_%f_(ip,j,k))
							Fp = max(qp,0.d0)*Ap + min(qp,0.d0)*Bp
					
							qm = 0.5d0*(u%ux_%f_(im,jp,k)*this%ptrMesh_%dyf_(jp)+ &
							            u%ux_%f_(im,j,k)*this%ptrMesh_%dyf_(j))/this%ptrMesh_%dyc_(jp)
							Am = u%uy_%f_(im,j,k)+0.5d0*psi_hmup*(u%uy_%f_(i,j,k)-u%uy_%f_(im,j,k))
							Bm = u%uy_%f_(i,j,k)+0.5d0*psi_hmum*(u%uy_%f_(im,j,k)-u%uy_%f_(i,j,k))
							Fm = max(qm,0.d0)*Am + min(qm,0.d0)*Bm
					
							this%phiY_(i,j,k) = this%phiY_(i,j,k) - (Fp-Fm)/this%ptrMesh_%dxf_(i)					
				
							!d (vv) / dy
							call compute_limiters(u%uy_%f_(i,j,k),u%uy_%f_(i,jm,k),u%uy_%f_(i,jp,k),&
											      u%uy_%f_(i,jpp,k),u%uy_%f_(i,jmm,k),&
											      c%f_(ip,j,k),c%f_(im,j,k),c%f_(i,jp,k),c%f_(i,j,k),&
											      c%f_(i,j,kp),c%f_(i,j,km),&
											      psi_hpup,psi_hpum,psi_hmup,psi_hmum,scheme)
							
							qp = 0.5d0*(u%uy_%f_(i,j,k)+u%uy_%f_(i,jp,k))
							Ap = u%uy_%f_(i,j,k)+0.5d0*psi_hpup*(u%uy_%f_(i,jp,k)-u%uy_%f_(i,j,k))
							Bp = u%uy_%f_(i,jp,k)+0.5d0*psi_hpum*(u%uy_%f_(i,j,k)-u%uy_%f_(i,jp,k))
							Fp = max(qp,0.d0)*Ap + min(qp,0.d0)*Bp
					
							qm = 0.5d0*(u%uy_%f_(i,j,k)+u%uy_%f_(i,jm,k))
							Am = u%uy_%f_(i,jm,k)+0.5d0*psi_hmup*(u%uy_%f_(i,j,k)-u%uy_%f_(i,jm,k))
							Bm = u%uy_%f_(i,j,k)+0.5d0*psi_hmum*(u%uy_%f_(i,jm,k)-u%uy_%f_(i,j,k))
							Fm = max(qm,0.d0)*Am + min(qm,0.d0)*Bm
					
							this%phiY_(i,j,k) = this%phiY_(i,j,k) - (Fp-Fm)/this%ptrMesh_%dyc_(jp)
					
							!d (wv) / dz
							call compute_limiters(u%uy_%f_(i,j,k),u%uy_%f_(i,j,km),u%uy_%f_(i,j,kp),&
											      u%uy_%f_(i,j,kpp),u%uy_%f_(i,j,kmm),&
											      c%f_(ip,j,k),c%f_(im,j,k),c%f_(i,jp,k),c%f_(i,j,k),&
											      c%f_(i,j,kp),c%f_(i,j,km),&
											      psi_hpup,psi_hpum,psi_hmup,psi_hmum,scheme)
											      
							qp = 0.5d0*(u%uz_%f_(i,jp,k)*this%ptrMesh_%dyf_(jp)+ &
									    u%uz_%f_(i,j,k)*this%ptrMesh_%dyf_(j))/this%ptrMesh_%dyc_(jp)
							Ap = u%uy_%f_(i,j,k)+0.5d0*psi_hpup*(u%uy_%f_(i,j,kp)-u%uy_%f_(i,j,k))
							Bp = u%uy_%f_(i,j,kp)+0.5d0*psi_hpum*(u%uy_%f_(i,j,k)-u%uy_%f_(i,j,kp))
							Fp = max(qp,0.d0)*Ap + min(qp,0.d0)*Bp
					
							qm = 0.5d0*(u%uz_%f_(i,jp,km)*this%ptrMesh_%dyf_(jp)+ &
									    u%uz_%f_(i,j,km)*this%ptrMesh_%dyf_(j))/this%ptrMesh_%dyc_(jp)
							Am = u%uy_%f_(i,j,km)+0.5d0*psi_hmup*(u%uy_%f_(i,j,k)-u%uy_%f_(i,j,km))
							Bm = u%uy_%f_(i,j,k)+0.5d0*psi_hmum*(u%uy_%f_(i,j,km)-u%uy_%f_(i,j,k))
							Fm = max(qm,0.d0)*Am + min(qm,0.d0)*Bm
					
							this%phiY_(i,j,k) = this%phiY_(i,j,k) - (Fp-Fm)/this%ptrMesh_%dzf_(k)
							
							!***********************  diffusion  ***********************!
							mur = 0.25d0*(mu%f_(i,j,k)+mu%f_(ip,j,k)+mu%f_(i,jp,k)+mu%f_(ip,jp,k))
							mul = 0.25d0*(mu%f_(i,j,k)+mu%f_(im,j,k)+mu%f_(i,jp,k)+mu%f_(im,jp,k))
							dxx = (  mur*(u%uy_%f_(ip,j,k)-u%uy_%f_(i,j,k))/(this%ptrMesh_%dxc_(ip))		&
								   - mul*(u%uy_%f_(i,j,k)-u%uy_%f_(im,j,k))/(this%ptrMesh_%dxc_(i))		&
								  )/(this%ptrMesh_%dxf_(i))
							dxxt = (  mur*(u%ux_%f_(i,jp,k)-u%ux_%f_(i,j,k))/(this%ptrMesh_%dyc_(jp))		&
								    - mul*(u%ux_%f_(im,jp,k)-u%ux_%f_(im,j,k))/(this%ptrMesh_%dyc_(jp))		&
								   )/(this%ptrMesh_%dxf_(i))
							
							mur = mu%f_(i,jp,k)
							mul = mu%f_(i,j,k)	  
							dyy = (  mur*(u%uy_%f_(i,jp,k)-u%uy_%f_(i,j,k))/(this%ptrMesh_%dyf_(jp))		&
								   - mul*(u%uy_%f_(i,j,k)-u%uy_%f_(i,jm,k))/(this%ptrMesh_%dyf_(j))		&
								  )/(this%ptrMesh_%dyc_(jp))
							dyyt = dyy
						
							mur = 0.25d0*(mu%f_(i,j,k)+mu%f_(i,j,kp)+mu%f_(i,jp,k)+mu%f_(i,jp,kp)) 
							mul = 0.25d0*(mu%f_(i,j,k)+mu%f_(i,j,km)+mu%f_(i,jp,k)+mu%f_(i,jp,km))	  
							dzz = (  mur*(u%uy_%f_(i,j,kp)-u%uy_%f_(i,j,k))/(this%ptrMesh_%dzc_(kp))		&
								   - mul*(u%uy_%f_(i,j,k)-u%uy_%f_(i,j,km))/(this%ptrMesh_%dzc_(k))		&
								  )/(this%ptrMesh_%dzf_(k))
							dzzt = (  mur*(u%uz_%f_(i,jp,k)-u%uz_%f_(i,j,k))/(this%ptrMesh_%dyc_(jp))		&
								    - mul*(u%uz_%f_(i,jp,km)-u%uz_%f_(i,j,km))/(this%ptrMesh_%dyc_(jp))		&
								   )/(this%ptrMesh_%dzf_(k))
								  
							
							invrho = 0.5d0*(1.d0/rho%f_(i,jp,k)+1.d0/rho%f_(i,j,k))
							this%phiY_(i,j,k) = this%phiY_(i,j,k) + (dxx + dyy + dzz +	&
											   						 dxxt + dyyt + dzzt)*invrho	
					
					
						end do
				end do
		end do
		!$OMP END PARALLEL DO 
		
					
		! z comp
		
 		!$OMP PARALLEL DO DEFAULT(none) &
		!$OMP SHARED(this,u,c,mu,rho,scheme) &
		!$OMP PRIVATE(i,j,k) &
		!$OMP PRIVATE(ip,im,ipp,imm) &
		!$OMP PRIVATE(jp,jm,jpp,jmm) &
		!$OMP PRIVATE(kp,km,kpp,kmm) &
		!$OMP PRIVATE(qp,Ap,Bp,Fp) &
		!$OMP PRIVATE(qm,Am,Bm,Fm) &
		!$OMP PRIVATE(dxx,dyy,dzz,dxxt,dyyt,dzzt) &
		!$OMP PRIVATE(psi_hpup,psi_hpum,psi_hmup,psi_hmum) &
		!$OMP PRIVATE(mur,mul,invrho)
		do k = this%ksz_,this%kez_
		
			kp = k + 1
			km = k - 1
			kpp = kp + 1
			kmm = km - 1
			
				do j = this%jsz_,this%jez_
			
					jp = j + 1
					jm = j - 1
					jpp = jp + 1
					jmm = jm - 1
				
						do i = this%isz_,this%iez_
				
							ip = i + 1
							im = i - 1
							ipp = ip + 1
							imm = im - 1

							!***********************  convection  ***********************!
							!d (uw) / dx
							call compute_limiters(u%uz_%f_(i,j,k),u%uz_%f_(im,j,k),u%uz_%f_(ip,j,k),&
											      u%uz_%f_(ipp,j,k),u%uz_%f_(imm,j,k),&
											      c%f_(ip,j,k),c%f_(im,j,k),c%f_(i,jp,k),c%f_(i,jm,k),&
											      c%f_(i,j,kp),c%f_(i,j,k),&
											      psi_hpup,psi_hpum,psi_hmup,psi_hmum,scheme)			
						
							qp = 0.5d0*(u%ux_%f_(i,j,k)+u%ux_%f_(i,j,kp))
							Ap = u%uz_%f_(i,j,k)+0.5d0*psi_hpup*(u%uz_%f_(ip,j,k)-u%uz_%f_(i,j,k))
							Bp = u%uz_%f_(ip,j,k)+0.5d0*psi_hpum*(u%uz_%f_(i,j,k)-u%uz_%f_(ip,j,k))
							Fp = max(qp,0.d0)*Ap + min(qp,0.d0)*Bp
					
							qm = 0.5d0*(u%ux_%f_(im,j,k)+u%ux_%f_(im,j,kp))
							Am = u%uz_%f_(im,j,k)+0.5d0*psi_hmup*(u%uz_%f_(i,j,k)-u%uz_%f_(im,j,k))
							Bm = u%uz_%f_(i,j,k)+0.5d0*psi_hmum*(u%uz_%f_(im,j,k)-u%uz_%f_(i,j,k))
							Fm = max(qm,0.d0)*Am + min(qm,0.d0)*Bm
					
							this%phiZ_(i,j,k) = this%phiZ_(i,j,k) - (Fp-Fm)/this%ptrMesh_%dxf_(i)	
					
							!d (vw) / dy
							call compute_limiters(u%uz_%f_(i,j,k),u%uz_%f_(i,jm,k),u%uz_%f_(i,jp,k),&
											      u%uz_%f_(i,jpp,k),u%uz_%f_(i,jmm,k),&
											      c%f_(ip,j,k),c%f_(im,j,k),c%f_(i,jp,k),c%f_(i,jm,k),&
											      c%f_(i,j,kp),c%f_(i,j,k),&
											      psi_hpup,psi_hpum,psi_hmup,psi_hmum,scheme)
											      
							qp = 0.5d0*(u%uy_%f_(i,j,k)+u%uy_%f_(i,j,kp))
							Ap = u%uz_%f_(i,j,k)+0.5d0*psi_hpup*(u%uz_%f_(i,jp,k)-u%uz_%f_(i,j,k))
							Bp = u%uz_%f_(i,jp,k)+0.5d0*psi_hpum*(u%uz_%f_(i,j,k)-u%uz_%f_(i,jp,k))
							Fp = max(qp,0.d0)*Ap + min(qp,0.d0)*Bp
					
							qm = 0.5d0*(u%uy_%f_(i,jm,k)+u%uy_%f_(i,jm,kp))
							Am = u%uz_%f_(i,jm,k)+0.5d0*psi_hmup*(u%uz_%f_(i,j,k)-u%uz_%f_(i,jm,k))
							Bm = u%uz_%f_(i,j,k)+0.5d0*psi_hmum*(u%uz_%f_(i,jm,k)-u%uz_%f_(i,j,k))
							Fm = max(qm,0.d0)*Am + min(qm,0.d0)*Bm
					
							this%phiZ_(i,j,k) = this%phiZ_(i,j,k) - (Fp-Fm)/this%ptrMesh_%dyf_(j)	
					
							!d (ww) / dz
							call compute_limiters(u%uz_%f_(i,j,k),u%uz_%f_(i,j,km),u%uz_%f_(i,j,kp),&
											      u%uz_%f_(i,j,kpp),u%uz_%f_(i,j,kmm),&
											      c%f_(ip,j,k),c%f_(im,j,k),c%f_(i,jp,k),c%f_(i,jm,k),&
											      c%f_(i,j,kp),c%f_(i,j,k),&
											      psi_hpup,psi_hpum,psi_hmup,psi_hmum,scheme)
							qp = 0.5d0*(u%uz_%f_(i,j,k)+u%uz_%f_(i,j,kp))
							Ap = u%uz_%f_(i,j,k)+0.5d0*psi_hpup*(u%uz_%f_(i,j,kp)-u%uz_%f_(i,j,k))
							Bp = u%uz_%f_(i,j,kp)+0.5d0*psi_hpum*(u%uz_%f_(i,j,k)-u%uz_%f_(i,j,kp))
							Fp = max(qp,0.d0)*Ap + min(qp,0.d0)*Bp
					
							qm = 0.5d0*(u%uz_%f_(i,j,k)+u%uz_%f_(i,j,km))
							Am = u%uz_%f_(i,j,km)+0.5d0*psi_hmup*(u%uz_%f_(i,j,k)-u%uz_%f_(i,j,km))
							Bm = u%uz_%f_(i,j,k)+0.5d0*psi_hmum*(u%uz_%f_(i,j,km)-u%uz_%f_(i,j,k))
							Fm = max(qm,0.d0)*Am + min(qm,0.d0)*Bm
					
							this%phiZ_(i,j,k) = this%phiZ_(i,j,k) - (Fp-Fm)/this%ptrMesh_%dzc_(kp)
							
							!***********************  diffusion  ***********************!
							mur = 0.25d0*(mu%f_(i,j,k)+mu%f_(i,j,kp)+mu%f_(ip,j,k)+mu%f_(ip,j,kp))
							mul = 0.25d0*(mu%f_(i,j,k)+mu%f_(i,j,kp)+mu%f_(im,j,k)+mu%f_(im,j,kp))
							dxx = (  mur*(u%uz_%f_(ip,j,k)-u%uz_%f_(i,j,k))/(this%ptrMesh_%dxc_(ip))		&
								   - mul*(u%uz_%f_(i,j,k)-u%uz_%f_(im,j,k))/(this%ptrMesh_%dxc_(i))		&
								  )/(this%ptrMesh_%dxf_(i))
							dxxt = (  mur*(u%ux_%f_(i,j,kp)-u%ux_%f_(i,j,k))/(this%ptrMesh_%dzc_(kp))		&
								    - mul*(u%ux_%f_(im,j,kp)-u%ux_%f_(im,j,k))/(this%ptrMesh_%dzc_(kp))		&
								   )/(this%ptrMesh_%dxf_(i))
								  
							mur = 0.25d0*(mu%f_(i,j,k)+mu%f_(i,j,kp)+mu%f_(i,jp,k)+mu%f_(i,jp,kp)) 
							mul = 0.25d0*(mu%f_(i,j,k)+mu%f_(i,j,kp)+mu%f_(i,jm,k)+mu%f_(i,jm,kp))	  
							dyy = (  mur*(u%uz_%f_(i,jp,k)-u%uz_%f_(i,j,k))/(this%ptrMesh_%dyc_(jp))		&
								   - mul*(u%uz_%f_(i,j,k)-u%uz_%f_(i,jm,k))/(this%ptrMesh_%dyc_(j))		&
								  )/(this%ptrMesh_%dyf_(j))
							dyyt = (  mur*(u%uy_%f_(i,j,kp)-u%uy_%f_(i,j,k))/(this%ptrMesh_%dzc_(kp))		&
								    - mul*(u%uy_%f_(i,jm,kp)-u%uy_%f_(i,jm,k))/(this%ptrMesh_%dzc_(kp))		&
								   )/(this%ptrMesh_%dyf_(j))
							
							mur = mu%f_(i,j,kp)
							mul = mu%f_(i,j,k)	  
							dzz = (  mur*(u%uz_%f_(i,j,kp)-u%uz_%f_(i,j,k))/(this%ptrMesh_%dzf_(kp))		&
								   - mul*(u%uz_%f_(i,j,k)-u%uz_%f_(i,j,km))/(this%ptrMesh_%dzf_(k))		&
								  )/(this%ptrMesh_%dzc_(kp))
							dzzt = dzz
										   		
							invrho = 0.5d0*(1.d0/rho%f_(i,j,kp)+1.d0/rho%f_(i,j,k))
							this%phiZ_(i,j,k) = this%phiZ_(i,j,k) + (dxx + dyy + dzz +	&
											   						 dxxt + dyyt + dzzt)*invrho					
					
					
						end do
				end do
		end do
		!$OMP END PARALLEL DO 
	
				
    end subroutine
!========================================================================================!

!========================================================================================!
    subroutine compute_limiters(vc,vm,vp,vpp,vmm,cxp,cxm,cyp,cym,czp,czm,&
    							psi_hpup,psi_hpum,psi_hmup,psi_hmum,scheme)
    	real(DP), intent(in) ::  vc,vm,vp,vpp,vmm,cxp,cxm,cyp,cym,czp,czm
    	real(DP), intent(out) ::  psi_hpup,psi_hpum,psi_hmup,psi_hmum
    	integer, intent(in) :: scheme
    	real(DP) :: r_hpup,r_hpum,r_hmup,r_hmum

    	select case(scheme)
    		case(0)
    			call compute_CD_limiter(psi_hpup)
    			call compute_CD_limiter(psi_hpum)
    			call compute_CD_limiter(psi_hmup)
    			call compute_CD_limiter(psi_hmum)
    			return
    		case(1)
    			call compute_UD_limiter(psi_hpup)
    			call compute_UD_limiter(psi_hpum)
    			call compute_UD_limiter(psi_hmup)
    			call compute_UD_limiter(psi_hmum)
    			return  		
    		case default
    	end select
  
		r_hpup = (vc-vm)/approx_zero(vp-vc)
		r_hpum = (vpp-vp)/approx_zero(vp-vc)
		r_hmup = (vm-vmm)/approx_zero(vc-vm)
		r_hmum = (vp-vc)/approx_zero(vc-vm)
    	
    	select case(scheme)
    		case(2)
    			call compute_QUICK_limiter(r_hpup,psi_hpup)
    			call compute_QUICK_limiter(r_hpum,psi_hpum)
    			call compute_QUICK_limiter(r_hmup,psi_hmup)
    			call compute_QUICK_limiter(r_hmum,psi_hmum)	
    		case(3)
    			call compute_VanLeer_limiter(r_hpup,psi_hpup)
    			call compute_VanLeer_limiter(r_hpum,psi_hpum)
    			call compute_VanLeer_limiter(r_hmup,psi_hmup)
    			call compute_VanLeer_limiter(r_hmum,psi_hmum)
    		case(4)
    			call compute_Superbee_limiter(r_hpup,psi_hpup)
    			call compute_Superbee_limiter(r_hpum,psi_hpum)
    			call compute_Superbee_limiter(r_hmup,psi_hmup)
    			call compute_Superbee_limiter(r_hmum,psi_hmum)
    		case(5)
    			call compute_localVF_limiter(r_hpup,psi_hpup,cxp,cxm,cyp,cym,czp,czm)
    			call compute_localVF_limiter(r_hpum,psi_hpum,cxp,cxm,cyp,cym,czp,czm)
    			call compute_localVF_limiter(r_hmup,psi_hmup,cxp,cxm,cyp,cym,czp,czm)
    			call compute_localVF_limiter(r_hmum,psi_hmum,cxp,cxm,cyp,cym,czp,czm)
    		case default
    	end select
    	
    end subroutine
!========================================================================================!

!========================================================================================!
    subroutine compute_CD_limiter(psi)
    	real(DP), intent(out) :: psi
			
		psi = 1.d0
			
    end subroutine
!========================================================================================!

!========================================================================================!
    subroutine compute_UD_limiter(psi)
    	real(DP), intent(out) :: psi
			
		psi = 0.d0
			
    end subroutine
!========================================================================================!

!========================================================================================!
    subroutine compute_QUICK_limiter(r,psi)
    	real(DP), intent(in) :: r
    	real(DP), intent(out) :: psi
			
		psi = (r+3.d0)/4.d0
			
    end subroutine
!========================================================================================!

!========================================================================================!
    subroutine compute_VanLeer_limiter(r,psi)
    	real(DP), intent(inout) :: r
    	real(DP), intent(out) :: psi
			
		r = max(0.d0,r)
    	psi = (r+abs(r))/(1.d0+abs(r))
    		
    end subroutine
!========================================================================================!

!========================================================================================!
    subroutine compute_Superbee_limiter(r,psi)
    	real(DP), intent(inout) :: r
    	real(DP), intent(out) :: psi
			
		r = max(0.d0,r)
    	psi = maxval((/0.d0,min(2.d0*r,1.d0),min(r,2.d0)/))
    		
    end subroutine
!========================================================================================!

!========================================================================================!
    subroutine compute_localVF_limiter(r,psi,cxp,cxm,cyp,cym,czp,czm)
    	real(DP), intent(inout) :: r
    	real(DP), intent(out) :: psi
    	real(DP), intent(in) :: cxp,cxm,cyp,cym,czp,czm
    	real(DP) :: dcx,dcy,dcz,small
    	
    	dcx=cxp-cxm
    	dcy=cyp-cym
    	dcz=czp-czm
    	small=epsilon(0.d0)
    	
    	if ((dcx>small).OR.(dcy>small).OR.(dcz>small)) then
    		call compute_Superbee_limiter(r,psi)
    	else
    		call compute_CD_limiter(psi)
    	end if
    	
    end subroutine
!========================================================================================!

!========================================================================================!
    subroutine addConvDiff(this,u)
    	type(momentumEqn), intent(in) :: this
    	type(vfield), intent(inout) :: u 
    	real(DP) :: dt, gamma, xi, alpha
    	integer :: i,j,k
    	
    	dt = this%ptrTime_%dt_
    	alpha = alphaRKS(this%ptrTime_)
    	gamma = gammaRKS(this%ptrTime_)
    	xi = xiRKS(this%ptrTime_)
    		
 		!$OMP PARALLEL DO DEFAULT(none) &
		!$OMP SHARED(this,u) &
		!$OMP SHARED(dt,alpha,gamma,xi) &
		!$OMP PRIVATE(i,j,k)   	
    	do k=this%ksx_,this%kex_
    		do j=this%jsx_,this%jex_
    			do i=this%isx_,this%iex_
    				u%ux_%f_(i,j,k) = u%ux_%f_(i,j,k) + &
    							      dt*(gamma*this%phiX_(i,j,k) + xi*this%phiPrevX_(i,j,k)) 								
    			end do
    		end do
    	end do
    	!$OMP END PARALLEL DO
    	
 		!$OMP PARALLEL DO DEFAULT(none) &
		!$OMP SHARED(this,u) &
		!$OMP SHARED(dt,alpha,gamma,xi) &
		!$OMP PRIVATE(i,j,k) 
    	do k=this%ksy_,this%key_
    		do j=this%jsy_,this%jey_
    			do i=this%isy_,this%iey_
 					u%uy_%f_(i,j,k) = u%uy_%f_(i,j,k) + &
    								  dt*(gamma*this%phiY_(i,j,k) + xi*this%phiPrevY_(i,j,k))    								
    			end do
    		end do
    	end do   
    	!$OMP END PARALLEL DO 
    	
 		!$OMP PARALLEL DO DEFAULT(none) &
		!$OMP SHARED(this,u) &
		!$OMP SHARED(dt,alpha,gamma,xi) &
		!$OMP PRIVATE(i,j,k) 
    	do k=this%ksz_,this%kez_
    		do j=this%jsz_,this%jez_
    			do i=this%isz_,this%iez_
					u%uz_%f_(i,j,k)  = u%uz_%f_(i,j,k) + &
    								   dt*(gamma*this%phiZ_(i,j,k) + xi*this%phiPrevZ_(i,j,k))  								
    			end do
    		end do
    	end do
    	!$OMP END PARALLEL DO 	
    	
    end subroutine
!========================================================================================!

!========================================================================================!
    subroutine addPressureGrad(this,u,p,rho)
    	type(momentumEqn), intent(in) :: this
    	type(vfield), intent(inout) :: u
    	type(field), intent(in) :: p, rho
    	real(DP) :: dt, alpha
    	real(DP) :: invrho
    	integer :: i, j, k
    	
    	dt = this%ptrTime_%dt_
    	alpha = alphaRKS(this%ptrTime_)
		
		!x comp
		
 		!$OMP PARALLEL DO DEFAULT(none) &
		!$OMP SHARED(this,u,p,rho) &
		!$OMP SHARED(dt,alpha) &
		!$OMP PRIVATE(i,j,k) &
		!$OMP PRIVATE(invrho)	
		do k = this%ksx_,this%kex_
			do j = this%jsx_,this%jex_
				do i = this%isx_,this%iex_
						
					invrho = 0.5d0*(1.d0/rho%f_(i+1,j,k)+1.d0/rho%f_(i,j,k))
							
					u%ux_%f_(i,j,k) = u%ux_%f_(i,j,k) - &
						dt*alpha*invrho*(p%f_(i+1,j,k)-p%f_(i,j,k))/(this%ptrMesh_%dxc_(i+1))
							
				end do
			end do
		end do
		!$OMP END PARALLEL DO
		
		!y-comp
		
 		!$OMP PARALLEL DO DEFAULT(none) &
		!$OMP SHARED(this,u,p,rho) &
		!$OMP SHARED(dt,alpha) &
		!$OMP PRIVATE(i,j,k) &
		!$OMP PRIVATE(invrho)	
		do k = this%ksy_,this%key_
			do j = this%jsy_,this%jey_
				do i = this%isy_,this%iey_
						
					invrho = 0.5d0*(1.d0/rho%f_(i,j+1,k)+1.d0/rho%f_(i,j,k))
						
					u%uy_%f_(i,j,k) = u%uy_%f_(i,j,k) - &
						dt*alpha*invrho*(p%f_(i,j+1,k)-p%f_(i,j,k))/(this%ptrMesh_%dyc_(j+1))
							
				end do
			end do
		end do
		!$OMP END PARALLEL DO
			
		!z-comp
		
 		!$OMP PARALLEL DO DEFAULT(none) &
		!$OMP SHARED(this,u,p,rho) &
		!$OMP SHARED(dt,alpha) &
		!$OMP PRIVATE(i,j,k) &
		!$OMP PRIVATE(invrho)
		do k = this%ksz_,this%kez_
			do j = this%jsz_,this%jez_
				do i = this%isz_,this%iez_
						
					invrho = 0.5d0*(1.d0/rho%f_(i,j,k+1)+1.d0/rho%f_(i,j,k))
							
					u%uz_%f_(i,j,k) = u%uz_%f_(i,j,k) - &
						dt*alpha*invrho*(p%f_(i,j,k+1)-p%f_(i,j,k))/(this%ptrMesh_%dzc_(k+1))
							
				end do
			end do
		end do
		!$OMP END PARALLEL DO

    end subroutine
    	
!========================================================================================!

!========================================================================================!
    subroutine addSource(this,u,rho,st)
    	type(momentumEqn), intent(in) :: this
    	type(vfield), intent(inout) :: u
    	type(field), intent(in) :: rho
    	type(vfield), intent(in) :: st
    	type(grid), pointer :: mesh
    	real(DP) :: dt, alpha
    	real(DP) :: invrho,r
    	integer :: i, j, k
    	
    	dt = this%ptrTime_%dt_
    	alpha = alphaRKS(this%ptrTime_)
    	
    	mesh => rho%ptrMesh_
		
		!x comp
 		!$OMP PARALLEL DO DEFAULT(none) &
		!$OMP SHARED(this,mesh,u,rho,st) &
		!$OMP SHARED(dt,alpha) &
		!$OMP PRIVATE(i,j,k) &
		!$OMP PRIVATE(invrho,r)	
		do k = this%ksx_,this%kex_
			do j = this%jsx_,this%jex_
				do i = this%isx_,this%iex_
						
					invrho = 0.5d0*(1.d0/rho%f_(i+1,j,k)+1.d0/rho%f_(i,j,k))
					r = invrho*dt*alpha
							
					u%ux_%f_(i,j,k) = u%ux_%f_(i,j,k) + r*st%ux_%f_(i,j,k) + r*this%fs_ &
									  + dt*alpha*this%gCH_
							
				end do
			end do
		end do
		!$OMP END PARALLEL DO
		
		!y comp
 		!$OMP PARALLEL DO DEFAULT(none) &
		!$OMP SHARED(this,mesh,u,rho,st) &
		!$OMP SHARED(dt,alpha) &
		!$OMP PRIVATE(i,j,k) &
		!$OMP PRIVATE(invrho,r)
		do k = this%ksy_,this%key_
			do j = this%jsy_,this%jey_
				do i = this%isy_,this%iey_
						
					invrho = 0.5d0*(1.d0/rho%f_(i,j+1,k)+1.d0/rho%f_(i,j,k))
					r = invrho*dt*alpha
				
					u%uy_%f_(i,j,k) = u%uy_%f_(i,j,k) + r*st%uy_%f_(i,j,k) + dt*alpha*this%g_
							
				end do
			end do
		end do
		!$OMP END PARALLEL DO
		
		!z comp
 		!$OMP PARALLEL DO DEFAULT(none) &
		!$OMP SHARED(this,mesh,u,rho,st) &
		!$OMP SHARED(dt,alpha) &
		!$OMP PRIVATE(i,j,k) &
		!$OMP PRIVATE(invrho,r)
		do k = this%ksz_,this%kez_
			do j = this%jsz_,this%jez_
				do i = this%isz_,this%iez_
						
					invrho = 0.5d0*(1.d0/rho%f_(i,j,k)+1.d0/rho%f_(i,j,k+1))
					r = invrho*dt*alpha
						
					u%uz_%f_(i,j,k) = u%uz_%f_(i,j,k) + r*st%uz_%f_(i,j,k)
							
				end do
			end do
		end do  
		!$OMP END PARALLEL DO  	

    end subroutine
    	
!========================================================================================!

!========================================================================================!
#ifdef MG_MODE

    subroutine makeVelocityDivFree(this,u,psi,rho)
    	type(momentumEqn), intent(in) :: this
    	type(vfield), intent(inout) :: u 
    	type(field), intent(in) :: psi, rho
    	real(DP) :: dt, alpha 
    	
    	dt = this%ptrTime_%dt_
    	alpha = alphaRKS(this%ptrTime_)
    	
		call addPressureGrad(this,u,psi,rho)	
    	!update boundaries
    	call updateBoundariesV(u)
    	
    end subroutine
    
#endif
    
#ifdef FAST_MODE

    subroutine makeVelocityDivFree(this,u,psi,rho,rho0,nl)
    	type(momentumEqn), intent(in) :: this
    	type(vfield), intent(inout) :: u 
    	type(field), intent(in) :: psi, rho
    	real(DP), intent(in) :: rho0
    	integer, intent(in) :: nl
    	type(grid), pointer :: mesh
    	real(DP) :: dt,alpha,r,invrho0,invrho,dpr,dpr0
    	integer :: i,j,k
    	
    	dt = this%ptrTime_%dt_
    	alpha = alphaRKS(this%ptrTime_)
    	r=dt*alpha
    	invrho0=1.d0/rho0
    	
    	mesh => psi%ptrMesh_
    	
		!x comp
		
 		!$OMP PARALLEL DO DEFAULT(none) &
		!$OMP SHARED(this,mesh,u,rho,psi) &
		!$OMP SHARED(r,invrho0,nl) &
		!$OMP PRIVATE(i,j,k) &
		!$OMP PRIVATE(invrho,dpr,dpr0)	
		do k = this%ksx_,this%kex_
			do j = this%jsx_,this%jex_
				do i = this%isx_,this%iex_
						
					invrho = 0.5d0*(1.d0/rho%f_(i+1,j,k)+1.d0/rho%f_(i,j,k))
							
					call computeOldPressGrad(dpr,dpr0,i,j,k,1,psi,mesh,nl)
					
					u%ux_%f_(i,j,k) = u%ux_%f_(i,j,k) - &
									  r*(invrho0*dpr+(invrho-invrho0)*dpr0)					
							
				end do
			end do
		end do
		!$OMP END PARALLEL DO
		
		!y-comp
		
 		!$OMP PARALLEL DO DEFAULT(none) &
		!$OMP SHARED(this,mesh,u,rho,psi) &
		!$OMP SHARED(r,invrho0,nl) &
		!$OMP PRIVATE(i,j,k) &
		!$OMP PRIVATE(invrho,dpr,dpr0)	
		do k = this%ksy_,this%key_
			do j = this%jsy_,this%jey_
				do i = this%isy_,this%iey_
						
					invrho = 0.5d0*(1.d0/rho%f_(i,j+1,k)+1.d0/rho%f_(i,j,k))
					
					call computeOldPressGrad(dpr,dpr0,i,j,k,2,psi,mesh,nl)	
					
					u%uy_%f_(i,j,k) = u%uy_%f_(i,j,k) - &
									  r*(invrho0*dpr+(invrho-invrho0)*dpr0)	
							
				end do
			end do
		end do
		!$OMP END PARALLEL DO
			
		!z-comp
		
 		!$OMP PARALLEL DO DEFAULT(none) &
		!$OMP SHARED(this,mesh,u,rho,psi) &
		!$OMP SHARED(r,invrho0,nl) &
		!$OMP PRIVATE(i,j,k) &
		!$OMP PRIVATE(invrho,dpr,dpr0)
		do k = this%ksz_,this%kez_
			do j = this%jsz_,this%jez_
				do i = this%isz_,this%iez_
						
					invrho = 0.5d0*(1.d0/rho%f_(i,j,k+1)+1.d0/rho%f_(i,j,k))
					
					call computeOldPressGrad(dpr,dpr0,i,j,k,3,psi,mesh,nl)
							
					u%uz_%f_(i,j,k) = u%uz_%f_(i,j,k) - &
									  r*(invrho0*dpr+(invrho-invrho0)*dpr0)
							
				end do
			end do
		end do
		!$OMP END PARALLEL DO
		
			
    	!update boundaries
    	call updateBoundariesV(u)
    	
    end subroutine
    
    subroutine computeOldPressGrad(dpr,dpr0,i,j,k,dir,psi,mesh,n)
    	real(DP), intent(out) :: dpr,dpr0
    	integer, intent(in) :: i,j,k,n,dir
    	type(field), intent(in) :: psi
    	type(grid), intent(in) :: mesh
    	
    	select case(n)
    		case(2)
    			select case(dir)
    				case(1)
    				    dpr=(psi%f_(i+1,j,k)-psi%f_(i,j,k))/mesh%dxc_(1)
    					dpr0=(2.d0*(psi%ptrf_%f_(i+1,j,k)-psi%ptrf_%f_(i,j,k))-&
    					    (psi%ptrf_%ptrf_%f_(i+1,j,k)-psi%ptrf_%ptrf_%f_(i,j,k)))/mesh%dxc_(1)
    				case(2)
    				    dpr=(psi%f_(i,j+1,k)-psi%f_(i,j,k))/mesh%dyc_(j+1)
    					dpr0=(2.d0*(psi%ptrf_%f_(i,j+1,k)-psi%ptrf_%f_(i,j,k))-&
    					    (psi%ptrf_%ptrf_%f_(i,j+1,k)-psi%ptrf_%ptrf_%f_(i,j,k)))/mesh%dyc_(j+1)
    				case(3)
    				    dpr=(psi%f_(i,j,k+1)-psi%f_(i,j,k))/mesh%dzc_(1)
    					dpr0=(2.d0*(psi%ptrf_%f_(i,j,k+1)-psi%ptrf_%f_(i,j,k))-&
    					    (psi%ptrf_%ptrf_%f_(i,j,k+1)-psi%ptrf_%ptrf_%f_(i,j,k)))/mesh%dzc_(1)
    				case default
    			end select
    		case default 					
    	end select

    end subroutine

#endif

!========================================================================================!

!========================================================================================!
    subroutine setEqnBounds(this,u)
    	type(momentumEqn), intent(inout) :: this
    	type(vfield), intent(in) :: u

		!ux
    	this%isx_ = u%ux_%is_
    	this%iex_ = u%ux_%ie_  	
    	this%jsx_ = u%ux_%js_
    	this%jex_ = u%ux_%je_
    	this%ksx_ = u%ux_%ks_
    	this%kex_ = u%ux_%ke_
    	
       	!check for fixed value bc
       	if (u%ux_%bLeft_%bType_ == s_fixedValue) then
       		this%isx_ = this%isx_ + 1
       	end if
       
       	if (u%ux_%bRight_%bType_ == s_fixedValue) then
       		this%iex_ = this%iex_ - 1
       	end if 
       	
		!uy
    	this%isy_ = u%uy_%is_
    	this%iey_ = u%uy_%ie_   	
    	this%jsy_ = u%uy_%js_
    	this%jey_ = u%uy_%je_
    	this%ksy_ = u%uy_%ks_
    	this%key_ = u%uy_%ke_
    	
       	!check for fixed value bc
       	if (u%uy_%bBottom_%bType_ == s_fixedValue) then
       		this%jsy_ = this%jsy_ + 1
       	end if
       
       	if (u%uy_%bTop_%bType_ == s_fixedValue) then
       		this%jey_ = this%jey_ - 1
       	end if 
       	
		!uz
    	this%isz_ = u%uz_%is_
    	this%iez_ = u%uz_%ie_    	
    	this%jsz_ = u%uz_%js_
    	this%jez_ = u%uz_%je_
    	this%ksz_ = u%uz_%ks_
    	this%kez_ = u%uz_%ke_
    	
       	!check for fixed value bc
       	if (u%uz_%bBack_%bType_ == s_fixedValue) then
       		this%ksz_ = this%ksz_ + 1
       	end if
       
       	if (u%uz_%bFront_%bType_ == s_fixedValue) then
       		this%kez_ = this%kez_ - 1
       	end if 
    	
        
    end subroutine
!========================================================================================!

!========================================================================================!
    subroutine initFluxes(this)
    	type(momentumEqn), intent(inout) :: this
    	
		this%phiX_ = 0.d0
		this%phiPrevX_ = 0.d0
		
		this%phiY_ = 0.d0
		this%phiPrevY_ = 0.d0
		
		this%phiZ_ = 0.d0
		this%phiPrevZ_ = 0.d0
		
        
    end subroutine
!========================================================================================!

!========================================================================================!
    subroutine resetFluxes(this)
    	type(momentumEqn), intent(inout) :: this
    	
    	call set2zero_omp(this%phiX_)
    	call set2zero_omp(this%phiY_)
    	call set2zero_omp(this%phiZ_)
			 
    end subroutine
!========================================================================================!

!========================================================================================!
    subroutine copyOldFluxes(this)
    	type(momentumEqn), intent(inout) :: this
    	  	
    	call assign_omp(this%phiPrevX_,this%phiX_)
    	call assign_omp(this%phiPrevY_,this%phiY_)
    	call assign_omp(this%phiPrevZ_,this%phiZ_)
    	    	
		 
    end subroutine
!========================================================================================!

!========================================================================================!
    subroutine info(this,cpuTime)
    	type(momentumEqn), intent(in) :: this
    	real(DP), intent(in) :: cpuTime
    	type(mpiControl), pointer :: comm
    	real(DP) :: cpuTime_max
    	integer :: ierror
    	
    	comm => this%ptrMesh_%ptrMPIC_
    	
    	call Mpi_Reduce(cpuTime, cpuTime_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, &
    			        comm%cartComm_, ierror)

		if (IS_MASTER) then
			write(*,'(A,'//s_outputFormat(2:9)//')') '	U Eqn: CPU time = ', cpuTime_max
		end if
		
    	
    end subroutine
!========================================================================================!

!========================================================================================!
    subroutine initOldTimeFlux(this,gMesh,mesh,gu,rt)
    	type(momentumEqn), intent(inout) :: this
    	type(grid), intent(in) :: gMesh,mesh
    	type(vfield), intent(in) :: gu
    	type(time), intent(in) :: rt

		if ((IS_MASTER) .AND. (rt%scheme_==s_AB2)) then
			call vfieldCTOR(gPhi0,'phi0',gMesh,'sx','sy','sz',1,initOpt=3,&
						    nFolder=rt%inputFold_)
			call copyBoundaryV(gPhi0,gu)
		end if
		
		!decompose fluxes
		if (rt%scheme_==s_AB2) then
			call vfieldCTOR(phi0,'phi0',mesh,'sx','sy','sz',1,initOpt=-1)
			call decomposeFieldV(gPhi0,phi0)
			this%phiX_ = phi0%ux_%f_(this%isx_:this%iex_,this%jsx_:this%jex_,this%ksx_:this%kex_)
			this%phiY_ = phi0%uy_%f_(this%isy_:this%iey_,this%jsy_:this%jey_,this%ksy_:this%key_)
			this%phiZ_ = phi0%uz_%f_(this%isz_:this%iez_,this%jsz_:this%jez_,this%ksz_:this%kez_)
		end if
		
    	
    end subroutine
!========================================================================================!

!========================================================================================!
    subroutine printOldTimeFlux(this,scheme,output_folder)
    	type(momentumEqn), intent(in) :: this
    	integer, intent(in) :: scheme,output_folder
		
		if (scheme==s_AB2) then
			phi0%ux_%f_(this%isx_:this%iex_,this%jsx_:this%jex_,this%ksx_:this%kex_) = this%phiX_
			phi0%uy_%f_(this%isy_:this%iey_,this%jsy_:this%jey_,this%ksy_:this%key_) = this%phiY_
			phi0%uz_%f_(this%isz_:this%iez_,this%jsz_:this%jez_,this%ksz_:this%kez_) = this%phiZ_
			call reconstructAndWriteFieldV(phi0,gPhi0,output_folder)
		end if
			
    end subroutine
!========================================================================================!



end module momentumEqnMod

