module momentumEqnMod
	
	use timeMod
	
	implicit none
	
	!storage old fluxes for AB2 scheme
	type(vectorField) :: gPhi0, phi0


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
		
		!convection scheme dictionary
		real(DP) :: k_
		
		!keep a pointer to time
		type(time), pointer :: ptrTime_ => NULL() 
		
		!uniform source term 
		real(DP) :: fs_, g_, gCH_
		
		!flow rate
		real(DP) :: Q0_
		
		!flow control
		integer :: flowCtrl_

	end type
	
	private :: updateConveDiff
	private :: addPressureGrad
	private :: addSource
	private :: addConvDiff
	private :: setEqnBounds
	private :: initFluxes
	private :: resetFluxes
	private :: copyOldFluxes
	private :: info
	private :: initOldTimeFlux
!DIR$ IF DEFINED (FAST_MODE)	
	private :: computeOldPressGrad
!DIR$ ENDIF
	
	public :: momentumEqnCTOR
	public :: makeVelocityDivFree
	public :: solveMomentumEqn
	public :: setFlowRate
	public :: printOldTimeFlux

  	
contains


!========================================================================================!
	subroutine momentumEqnCTOR(this,gMesh,mesh,u,rt)
		type(momentumEqn) :: this
		type(grid), intent(in), target :: gMesh,mesh
        type(vectorField), intent(in) :: u
        type(time), intent(in), target :: rt
        type(dictionary) :: dict_conv, dict_flow, dict_g
        
        call dictionaryCTOR(dict_conv,'schemes','specs')
        !read k convection scheme parameter
        call readParameter(dict_conv,this%k_,'k')
        
        call dictionaryCTOR(dict_flow,'flowControl','specs')
        call readParameter(dict_flow,this%flowCtrl_,'flowCtrl')

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
        
        if (this%flowCtrl_==2) then
        	call computeFlowRate(u,this%Q0_)      	
        end if
        
        !read gravity
		call dictionaryCTOR(dict_g,'parameters','specs')
        call readParameter(dict_g,this%g_,'g')
        call readParameter(dict_g,this%gCH_,'gCH')
        
        !init u0 if AB2
        call initOldTimeFlux(this,gMesh,mesh,rt)
					    
	end subroutine
!========================================================================================!

!========================================================================================!
    subroutine solveMomentumEqn(this,u,p,mu,rho,st,c)
    	type(momentumEqn), intent(inout) :: this
    	type(vectorField), intent(inout) :: u
    	type(scalarField), intent(in) :: p, mu, rho, c
    	type(vectorField), intent(in) :: st
    	real(DP) :: start, finish
    	

    	start = MPI_Wtime() 
    	
    	!copy fluxes to Prev arrays
		call copyOldFluxes(this)
			
    	!reset fluxes for computation of new fluxes
    	call resetFluxes(this)
    	
    	!update fluxes
    	call updateConveDiff(this,u,mu,rho)
		call addConvDiff(this,u)
    		
    	!add pressure grad
    	!call addPressureGrad(this,u,p,rho)
    	
    	!add sources
    	call addSource(this,u,rho,st,c)
    	
    	!update boundaries
    	call updateBoundariesV(u)
    	
    	finish = MPI_Wtime()
    	
    	!print-out info
    	call info(this,finish-start)
    	
        
    end subroutine
!========================================================================================!

!========================================================================================!
    subroutine updateConveDiff(this,u,mu,rho)
    	type(momentumEqn), intent(inout) :: this
    	type(vectorField), intent(in) :: u
    	type(scalarField), intent(in) :: mu, rho
    	integer :: im, imm, jm, jmm, km, kmm
    	integer :: ip, ipp, jp, jpp, kp, kpp
    	real(DP) :: qp, Ap, Bp, Fp
    	real(DP) :: qm, Am, Bm, Fm
    	real(DP) :: mur, mul, dxx, dyy, dzz, dxxt, dyyt, dzzt
    	integer :: i, j, k
    	real(DP) :: r,invrho

    	
		!k-scheme factor
		r = 0.25d0*(1.d0 - this%k_)
		
		
		! x comp
		
 		!$OMP PARALLEL DO DEFAULT(none) &
		!$OMP SHARED(this,u,mu,rho,r) &
		!$OMP PRIVATE(i,j,k) &
		!$OMP PRIVATE(ip,im,ipp,imm) &
		!$OMP PRIVATE(jp,jm,jpp,jmm) &
		!$OMP PRIVATE(kp,km,kpp,kmm) &
		!$OMP PRIVATE(qp,Ap,Bp,Fp) &
		!$OMP PRIVATE(qm,Am,Bm,Fm) &
		!$OMP PRIVATE(dxx,dyy,dzz,dxxt,dyyt,dzzt) &
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
							qp = 0.5d0*(u%ux_%f_(i,j,k)+u%ux_%f_(ip,j,k))
							Ap = r*(-u%ux_%f_(im,j,k)+2.d0*u%ux_%f_(i,j,k)-u%ux_%f_(ip,j,k))
							Bp = r*(-u%ux_%f_(i,j,k)+2.d0*u%ux_%f_(ip,j,k)-u%ux_%f_(ipp,j,k))
							Fp = qp*0.5d0*(u%ux_%f_(i,j,k)+u%ux_%f_(ip,j,k)) + max(qp,0.d0)*Ap + min(qp,0.d0)*Bp
					
							qm = 0.5d0*(u%ux_%f_(i,j,k)+u%ux_%f_(im,j,k))
							Am = r*(-u%ux_%f_(imm,j,k)+2.d0*u%ux_%f_(im,j,k)-u%ux_%f_(i,j,k))
							Bm = r*(-u%ux_%f_(im,j,k)+2.d0*u%ux_%f_(i,j,k)-u%ux_%f_(ip,j,k))
							Fm = qm*0.5d0*(u%ux_%f_(i,j,k)+u%ux_%f_(im,j,k)) + max(qm,0.d0)*Am + min(qm,0.d0)*Bm
					
							this%phiX_(i,j,k) = this%phiX_(i,j,k) - (Fp-Fm)/this%ptrMesh_%dxc_(ip)
					
							!d (vu) / dy
							qp = 0.5d0*(u%uy_%f_(ip,j,k)+u%uy_%f_(i,j,k))
							Ap = r*(-u%ux_%f_(i,jm,k)+2.d0*u%ux_%f_(i,j,k)-u%ux_%f_(i,jp,k))
							Bp = r*(-u%ux_%f_(i,j,k)+2.d0*u%ux_%f_(i,jp,k)-u%ux_%f_(i,jpp,k))
							Fp = qp*0.5d0*(u%ux_%f_(i,j,k)+u%ux_%f_(i,jp,k)) + max(qp,0.d0)*Ap + min(qp,0.d0)*Bp
					
							qm = 0.5d0*(u%uy_%f_(ip,jm,k)+u%uy_%f_(i,jm,k))
							Ap = r*(-u%ux_%f_(i,jmm,k)+2.d0*u%ux_%f_(i,jm,k)-u%ux_%f_(i,j,k))
							Bp = r*(-u%ux_%f_(i,jm,k)+2.d0*u%ux_%f_(i,j,k)-u%ux_%f_(i,jp,k))
							Fm = qm*0.5d0*(u%ux_%f_(i,j,k)+u%ux_%f_(i,jm,k)) + max(qm,0.d0)*Ap + min(qm,0.d0)*Bp
					
							this%phiX_(i,j,k) = this%phiX_(i,j,k) - (Fp-Fm)/this%ptrMesh_%dyf_(j)
					
							!d (wu) / dz
							qp = 0.5d0*(u%uz_%f_(i,j,k)+u%uz_%f_(ip,j,k))
							Ap = r*(-u%ux_%f_(i,j,km)+2.d0*u%ux_%f_(i,j,k)-u%ux_%f_(i,j,kp))
							Bp = r*(-u%ux_%f_(i,j,k)+2.d0*u%ux_%f_(i,j,kp)-u%ux_%f_(i,j,kpp))
							Fp = qp*0.5d0*(u%ux_%f_(i,j,k)+u%ux_%f_(i,j,kp)) + max(qp,0.d0)*Ap + min(qp,0.d0)*Bp
					
							qm = 0.5d0*(u%uz_%f_(i,j,km)+u%uz_%f_(ip,j,km))
							Ap = r*(-u%ux_%f_(i,j,kmm)+2.d0*u%ux_%f_(i,j,km)-u%ux_%f_(i,j,k))
							Bp = r*(-u%ux_%f_(i,j,km)+2.d0*u%ux_%f_(i,j,k)-u%ux_%f_(i,j,kp))
							Fm = qm*0.5d0*(u%ux_%f_(i,j,k)+u%ux_%f_(i,j,km)) + max(qm,0.d0)*Ap + min(qm,0.d0)*Bp
					
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
		!$OMP SHARED(this,u,mu,rho,r) &
		!$OMP PRIVATE(i,j,k) &
		!$OMP PRIVATE(ip,im,ipp,imm) &
		!$OMP PRIVATE(jp,jm,jpp,jmm) &
		!$OMP PRIVATE(kp,km,kpp,kmm) &
		!$OMP PRIVATE(qp,Ap,Bp,Fp) &
		!$OMP PRIVATE(qm,Am,Bm,Fm) &
		!$OMP PRIVATE(dxx,dyy,dzz,dxxt,dyyt,dzzt) &
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
							qp = 0.5d0*(u%ux_%f_(i,j,k)+u%ux_%f_(i,jp,k))
							Ap = r*(-u%uy_%f_(im,j,k)+2.d0*u%uy_%f_(i,j,k)-u%uy_%f_(ip,j,k))
							Bp = r*(-u%uy_%f_(i,j,k)+2.d0*u%uy_%f_(ip,j,k)-u%uy_%f_(ipp,j,k))
							Fp = qp*0.5d0*(u%uy_%f_(i,j,k)+u%uy_%f_(ip,j,k)) + max(qp,0.d0)*Ap + min(qp,0.d0)*Bp
					
							qm = 0.5d0*(u%ux_%f_(im,j,k)+u%ux_%f_(im,jp,k))
							Am = r*(-u%uy_%f_(imm,j,k)+2.d0*u%uy_%f_(im,j,k)-u%uy_%f_(i,j,k))
							Bm = r*(-u%uy_%f_(im,j,k)+2.d0*u%uy_%f_(i,j,k)-u%uy_%f_(ip,j,k))
							Fm = qm*0.5d0*(u%uy_%f_(i,j,k)+u%uy_%f_(im,j,k)) + max(qm,0.d0)*Ap + min(qm,0.d0)*Bp
					
							this%phiY_(i,j,k) = this%phiY_(i,j,k) - (Fp-Fm)/this%ptrMesh_%dxf_(i)				
				
							!d (vv) / dy
							qp = 0.5d0*(u%uy_%f_(i,j,k)+u%uy_%f_(i,jp,k))
							Ap = r*(-u%uy_%f_(i,jm,k)+2.d0*u%uy_%f_(i,j,k)-u%uy_%f_(i,jp,k))
							Bp = r*(-u%uy_%f_(i,j,k)+2.d0*u%uy_%f_(i,jp,k)-u%uy_%f_(i,jpp,k))
							Fp = qp*0.5d0*(u%uy_%f_(i,j,k)+u%uy_%f_(i,jp,k)) + max(qp,0.d0)*Ap + min(qp,0.d0)*Bp
					
							qm = 0.5d0*(u%uy_%f_(i,j,k)+u%uy_%f_(i,jm,k))
							Am = r*(-u%uy_%f_(i,jmm,k)+2.d0*u%uy_%f_(i,jm,k)-u%uy_%f_(i,j,k))
							Bm = r*(-u%uy_%f_(i,jm,k)+2.d0*u%uy_%f_(i,j,k)-u%uy_%f_(i,jp,k))
							Fm = qm*0.5d0*(u%uy_%f_(i,j,k)+u%uy_%f_(i,jm,k)) + max(qm,0.d0)*Ap + min(qm,0.d0)*Bp
					
							this%phiY_(i,j,k) = this%phiY_(i,j,k) - (Fp-Fm)/this%ptrMesh_%dyc_(jp)
					
							!d (wv) / dz
							qp = 0.5d0*(u%uz_%f_(i,j,k)+u%uz_%f_(i,jp,k))
							Ap = r*(-u%uy_%f_(i,j,km)+2.d0*u%uy_%f_(i,j,k)-u%uy_%f_(i,j,kp))
							Bp = r*(-u%uy_%f_(i,j,k)+2.d0*u%uy_%f_(i,j,kp)-u%uy_%f_(i,j,kpp))
							Fp = qp*0.5d0*(u%uy_%f_(i,j,k)+u%uy_%f_(i,j,kp)) + max(qp,0.d0)*Ap + min(qp,0.d0)*Bp
					
							qm = 0.5d0*(u%uz_%f_(i,j,km)+u%uz_%f_(i,jp,km))
							Am = r*(-u%uy_%f_(i,j,kmm)+2.d0*u%uy_%f_(i,j,km)-u%uy_%f_(i,j,k))
							Bm = r*(-u%uy_%f_(i,j,km)+2.d0*u%uy_%f_(i,j,k)-u%uy_%f_(i,j,kp))
							Fm = qm*0.5d0*(u%uy_%f_(i,j,k)+u%uy_%f_(i,j,km)) + max(qm,0.d0)*Ap + min(qm,0.d0)*Bp
					
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
		!$OMP SHARED(this,u,mu,rho,r) &
		!$OMP PRIVATE(i,j,k) &
		!$OMP PRIVATE(ip,im,ipp,imm) &
		!$OMP PRIVATE(jp,jm,jpp,jmm) &
		!$OMP PRIVATE(kp,km,kpp,kmm) &
		!$OMP PRIVATE(qp,Ap,Bp,Fp) &
		!$OMP PRIVATE(qm,Am,Bm,Fm) &
		!$OMP PRIVATE(dxx,dyy,dzz,dxxt,dyyt,dzzt) &
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
							qp = 0.5d0*(u%ux_%f_(i,j,k)+u%ux_%f_(i,j,kp))
							Ap = r*(-u%uz_%f_(im,j,k)+2.d0*u%uz_%f_(i,j,k)-u%uz_%f_(ip,j,k))
							Bp = r*(-u%uz_%f_(i,j,k)+2.d0*u%uz_%f_(ip,j,k)-u%uz_%f_(ipp,j,k))
							Fp = qp*0.5d0*(u%uz_%f_(i,j,k)+u%uz_%f_(ip,j,k)) + max(qp,0.d0)*Ap + min(qp,0.d0)*Bp
					
							qm = 0.5d0*(u%ux_%f_(im,j,k)+u%ux_%f_(im,j,kp))
							Am = r*(-u%uz_%f_(imm,j,k)+2.d0*u%uz_%f_(im,j,k)-u%uz_%f_(i,j,k))
							Bm = r*(-u%uz_%f_(im,j,k)+2.d0*u%uz_%f_(i,j,k)-u%uz_%f_(ip,j,k))
							Fm = qm*0.5d0*(u%uz_%f_(i,j,k)+u%uz_%f_(im,j,k)) + max(qm,0.d0)*Ap + min(qm,0.d0)*Bp
					
							this%phiZ_(i,j,k) = this%phiZ_(i,j,k) - (Fp-Fm)/this%ptrMesh_%dxf_(i)	
					
							!d (vw) / dy
							qp = 0.5d0*(u%uy_%f_(i,j,k)+u%uy_%f_(i,j,kp))
							Ap = r*(-u%uz_%f_(i,jm,k)+2.d0*u%uz_%f_(i,j,k)-u%uz_%f_(i,jp,k))
							Bp = r*(-u%uz_%f_(i,j,k)+2.d0*u%uz_%f_(i,jp,k)-u%uz_%f_(i,jpp,k))
							Fp = qp*0.5d0*(u%uz_%f_(i,j,k)+u%uz_%f_(i,jp,k)) + max(qp,0.d0)*Ap + min(qp,0.d0)*Bp
					
							qm = 0.5d0*(u%uy_%f_(i,jm,k)+u%uy_%f_(i,jm,kp))
							Am = r*(-u%uz_%f_(i,jmm,k)+2.d0*u%uz_%f_(i,jm,k)-u%uz_%f_(i,j,k))
							Bm = r*(-u%uz_%f_(i,jm,k)+2.d0*u%uz_%f_(i,j,k)-u%uz_%f_(i,jp,k))
							Fm = qm*0.5d0*(u%uz_%f_(i,j,k)+u%uz_%f_(i,jm,k)) + max(qm,0.d0)*Ap + min(qm,0.d0)*Bp
					
							this%phiZ_(i,j,k) = this%phiZ_(i,j,k) - (Fp-Fm)/this%ptrMesh_%dyf_(j)	
					
							!d (ww) / dz
							qp = 0.5d0*(u%uz_%f_(i,j,k)+u%uz_%f_(i,j,kp))
							Ap = r*(-u%uz_%f_(i,j,km)+2.d0*u%uz_%f_(i,j,k)-u%uz_%f_(i,j,kp))
							Bp = r*(-u%uz_%f_(i,j,k)+2.d0*u%uz_%f_(i,j,kp)-u%uz_%f_(i,j,kpp))
							Fp = qp*0.5d0*(u%uz_%f_(i,j,k)+u%uz_%f_(i,j,kp)) + max(qp,0.d0)*Ap + min(qp,0.d0)*Bp
					
							qm = 0.5d0*(u%uz_%f_(i,j,k)+u%uz_%f_(i,j,km))
							Am = r*(-u%uz_%f_(i,j,kmm)+2.d0*u%uz_%f_(i,j,km)-u%uz_%f_(i,j,k))
							Bm = r*(-u%uz_%f_(i,j,km)+2.d0*u%uz_%f_(i,j,k)-u%uz_%f_(i,j,kp))
							Fm = qm*0.5d0*(u%uz_%f_(i,j,k)+u%uz_%f_(i,j,km)) + max(qm,0.d0)*Ap + min(qm,0.d0)*Bp
					
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
    subroutine addConvDiff(this,u)
    	type(momentumEqn), intent(in) :: this
    	type(vectorField), intent(inout) :: u 
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
    	type(vectorField), intent(inout) :: u
    	type(scalarField), intent(in) :: p, rho
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
    subroutine addSource(this,u,rho,st,c)
    	type(momentumEqn), intent(in) :: this
    	type(vectorField), intent(inout) :: u
    	type(scalarField), intent(in) :: rho, c
    	type(vectorField), intent(in) :: st
    	type(grid), pointer :: mesh
    	real(DP) :: dt, alpha
    	real(DP) :: invrho,r
    	integer :: i, j, k
    	
    	dt = this%ptrTime_%dt_
    	alpha = alphaRKS(this%ptrTime_)
    	
    	mesh => c%ptrMesh_
		
		!x comp
 		!$OMP PARALLEL DO DEFAULT(none) &
		!$OMP SHARED(this,mesh,u,rho,st,c) &
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
		!$OMP SHARED(this,mesh,u,rho,st,c) &
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
		!$OMP SHARED(this,mesh,u,rho,st,c) &
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
!DIR$ IF DEFINED (MG_MODE)
    subroutine makeVelocityDivFree(this,u,psi,rho)
    	type(momentumEqn), intent(in) :: this
    	type(vectorField), intent(inout) :: u 
    	type(scalarField), intent(in) :: psi, rho
    	real(DP) :: dt, alpha 
    	
    	dt = this%ptrTime_%dt_
    	alpha = alphaRKS(this%ptrTime_)
    	
		call addPressureGrad(this,u,psi,rho)	
    	!update boundaries
    	call updateBoundariesV(u)
    	
    end subroutine
    
!DIR$ ELSEIF DEFINED (FAST_MODE)

    subroutine makeVelocityDivFree(this,u,psi,rho,rho0,nl)
    	type(momentumEqn), intent(in) :: this
    	type(vectorField), intent(inout) :: u 
    	type(scalarField), intent(in) :: psi, rho
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
    	type(scalarField), intent(in) :: psi
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
!DIR$ ENDIF
!========================================================================================!

!========================================================================================!
    subroutine setEqnBounds(this,u)
    	type(momentumEqn), intent(inout) :: this
    	type(vectorField), intent(in) :: u

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
    subroutine initOldTimeFlux(this,gMesh,mesh,rt)
    	type(momentumEqn), intent(inout) :: this
    	type(grid), intent(in) :: gMesh,mesh
    	type(time), intent(in) :: rt

		if ((IS_MASTER) .AND. (rt%scheme_==s_AB2)) then
			call vectorFieldCTOR(gPhi0,'phi0',gMesh,'sx','sy','sz',1,initOpt=1,&
								 nFolder=rt%inputFold_)
		end if
		
		!decompose fluxes
		if (rt%scheme_==s_AB2) then
			call vectorFieldCTOR(phi0,'phi0',mesh,'sx','sy','sz',1,initOpt=-1)
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

