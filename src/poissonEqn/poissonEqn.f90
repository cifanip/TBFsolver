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

module poissonEqnMod

#ifdef FAST_MODE	
	use fastPoissonSolverMod
#endif
#ifdef MG_MODE
	use pcgMod
#endif
	use timeMod
	
	implicit none


	type, public :: poissonEqn
	
		!keep a pointer to grid
		type(grid), pointer :: ptrMesh_ => NULL() 
		
#ifdef FAST_MODE		
		type(fastPoissonSolver), private :: fftSolver_
#endif
#ifdef MG_MODE
		type(pcg) :: pcgs_
#endif
		
		!poisson eqn source
		type(field), private :: s_

#ifdef MG_MODE
		!poisson eqn coefficients
		type(field), private :: beta_
#endif
		
		!keep a pointer to time
		type(time), pointer :: ptrTime_ => NULL()

		
		!number of pressure levels
		integer :: nl_

		!reference density (only used in FFT based solver)
		real(DP) :: rho0_ 


	end type
	
	!pressure fields
	type(field) :: gp, p, gp0, p0
	
	
	private :: computeSource
#ifdef MG_MODE
	private :: computeBeta
#endif
	private :: info
#ifdef FAST_MODE
	private :: computeOldPressDiv
#endif
	public :: poissonEqnCTOR	
	public :: solvePoissonEqn
	public :: updatePressure


  	
contains


!========================================================================================!
#ifdef FAST_MODE
	subroutine poissonEqnCTOR(this,mesh,gMesh,gpsi,psi,rt,rhol,rhog)
		type(poissonEqn) :: this
		type(grid), intent(in), target :: mesh,gMesh
		type(field), intent(in) :: gpsi
		type(field), intent(inout) :: psi
		type(time), intent(in), target :: rt
		real(DP), intent(in) :: rhol,rhog
		integer :: n_old
#endif
#ifdef MG_MODE
	subroutine poissonEqnCTOR(this,mesh,gMesh,c,gpsi,psi,rt,rhol,rhog)
		type(poissonEqn) :: this
		type(grid), intent(in), target :: mesh,gMesh
		type(field), intent(inout) :: c
		type(field), intent(in) :: gpsi
		type(field), intent(inout) :: psi
		type(time), intent(in), target :: rt
		real(DP), intent(in) :: rhol,rhog
#endif
		
		this%ptrMesh_ => mesh
		this%ptrTime_ => rt

		!init source and beta coeff
		call fieldCTOR(this%s_,'s',mesh,'cl',psi%hd_,initOpt=-1)
		call copyBoundary(this%s_,psi)
		
#ifdef FAST_MODE
		call fastPoissonSolverCTOR(this%fftSolver_,mesh,gMesh)
#endif		
#ifdef MG_MODE
		call fieldCTOR(this%beta_,'beta',mesh,'cl',psi%hd_,initOpt=-1)
		call copyBoundary(this%beta_,c)

		call pcgCTOR(this%pcgs_,mesh,psi,this%beta_)		
#endif	

		!set reference density
		this%rho0_=min(rhol,rhog)

		!init current and old pressure field
		if (IS_MASTER) then
			call fieldCTOR(gp,'p',gMesh,'cl',halo_size=1,initOpt=4,nFolder=rt%inputFold_)
			call copyBoundary(gp,gpsi,build_htypes=.FALSE.)
		end if
		call fieldCTOR(p,'p',mesh,'cl',halo_size=1,initOpt=-1)
		call decomposeField(gp,p)
		
		if (IS_MASTER) then
			call fieldCTOR(gp0,'p0',gMesh,'cl',halo_size=1,initOpt=4,nFolder=rt%inputFold_)
			call copyBoundary(gp0,gpsi,build_htypes=.FALSE.)
		end if
		call fieldCTOR(p0,'p0',mesh,'cl',halo_size=1,initOpt=-1)
		call decomposeField(gp0,p0)
		
		!allocate old pressure field
		this%nl_=2
		call allocateOldField(psi,this%nl_)
		psi%ptrf_%ptrf_%f_=p0%f_
		
			    
	end subroutine
!========================================================================================!

!========================================================================================!
#ifdef MG_MODE
    subroutine computeSource(this,u)
        type(poissonEqn), intent(inout) :: this
        type(vfield), intent(in) :: u
        type(grid), pointer :: mesh
        integer :: nx, ny, nz
        integer :: i, j, k
        real(DP) :: dx, dy, dz
        
        
        mesh => u%ptrMesh_
        
        nx = mesh%nx_
        ny = mesh%ny_
        nz = mesh%nz_
        
 		!$OMP PARALLEL DO DEFAULT(none) &
		!$OMP SHARED(this,u,mesh) &
		!$OMP SHARED(nx,ny,nz) &
		!$OMP PRIVATE(i,j,k) &
		!$OMP PRIVATE(dx,dy,dz)
        do k = 1,nz
        	do j = 1,ny
        		do i = 1,nx
        		
        			dx = (u%ux_%f_(i,j,k)-u%ux_%f_(i-1,j,k))/mesh%dxf_(i)
        			dy = (u%uy_%f_(i,j,k)-u%uy_%f_(i,j-1,k))/mesh%dyf_(j)
        			dz = (u%uz_%f_(i,j,k)-u%uz_%f_(i,j,k-1))/mesh%dzf_(k)
        			
        			this%s_%f_(i,j,k) = dx + dy + dz
        			
        		end do
        	end do
        end do
        !$OMP END PARALLEL DO 
        
    end subroutine
#endif   
#ifdef FAST_MODE
    subroutine computeSource(this,psi,rho,u,st)
        type(poissonEqn), intent(inout) :: this
        type(field), intent(in) :: psi,rho
        type(vfield), intent(in) :: u,st
        type(grid), pointer :: mesh
        integer :: nx, ny, nz
        integer :: i, j, k
        real(DP) :: dx,dy,dz,rho0,dt,alpha,r
        real(DP) :: apx,amx,apy,amy,apz,amz
        real(DP) :: dpxp,dpxm,dpyp,dpym,dpzp,dpzm
        
        
        mesh => u%ptrMesh_
        
        nx = mesh%nx_
        ny = mesh%ny_
        nz = mesh%nz_
        
        rho0=this%rho0_
        
    	dt = this%ptrTime_%dt_
    	alpha = alphaRKS(this%ptrTime_) 
    	r = 1.d0/(alpha*dt)
        
 		!$OMP PARALLEL DO DEFAULT(none) &
		!$OMP SHARED(this,u,mesh,psi,rho,st) &
		!$OMP SHARED(nx,ny,nz,r,rho0) &
		!$OMP PRIVATE(i,j,k) &
		!$OMP PRIVATE(dx,dy,dz) &
		!$OMP PRIVATE(apx,amx,apy,amy,apz,amz) &
		!$OMP PRIVATE(dpxp,dpxm,dpyp,dpym,dpzp,dpzm)
        do k = 1,nz
        	do j = 1,ny
        		do i = 1,nx
        		
        			dx = (u%ux_%f_(i,j,k)-u%ux_%f_(i-1,j,k))/mesh%dxf_(1)
        			dy = (u%uy_%f_(i,j,k)-u%uy_%f_(i,j-1,k))/mesh%dyf_(j)
        			dz = (u%uz_%f_(i,j,k)-u%uz_%f_(i,j,k-1))/mesh%dzf_(1)
        			
        			this%s_%f_(i,j,k) = (dx + dy + dz)*r*rho0
        			
        			apx=1.d0-rho0*0.5d0*(1.d0/rho%f_(i+1,j,k)+1.d0/rho%f_(i,j,k))
        			amx=1.d0-rho0*0.5d0*(1.d0/rho%f_(i-1,j,k)+1.d0/rho%f_(i,j,k))
        			apy=1.d0-rho0*0.5d0*(1.d0/rho%f_(i,j+1,k)+1.d0/rho%f_(i,j,k))
        			amy=1.d0-rho0*0.5d0*(1.d0/rho%f_(i,j-1,k)+1.d0/rho%f_(i,j,k))
        			apz=1.d0-rho0*0.5d0*(1.d0/rho%f_(i,j,k+1)+1.d0/rho%f_(i,j,k))
        			amz=1.d0-rho0*0.5d0*(1.d0/rho%f_(i,j,k-1)+1.d0/rho%f_(i,j,k))
        			
					call computeOldPressDiv(dpxp,dpxm,dpyp,dpym,dpzp,dpzm,&
										     i,j,k,psi,st,mesh,this%nl_)
					
					dx = (apx*dpxp-amx*dpxm)/mesh%dxf_(1)
					dy = (apy*dpyp-amy*dpym)/mesh%dyf_(j)
					dz = (apz*dpzp-amz*dpzm)/mesh%dzf_(1)
					
					this%s_%f_(i,j,k)=this%s_%f_(i,j,k)+dx+dy+dz
        			

        		end do
        	end do
        end do
        !$OMP END PARALLEL DO 
        
    end subroutine
    
    subroutine computeOldPressDiv(dpxp,dpxm,dpyp,dpym,dpzp,dpzm,i,j,k,psi,st,mesh,n)
    	real(DP), intent(out) :: dpxp,dpxm,dpyp,dpym,dpzp,dpzm
    	integer, intent(in) :: i,j,k,n
    	type(field), intent(in) :: psi
    	type(vfield), intent(in) :: st
    	type(grid), intent(in) :: mesh
    	
    	select case(n)
    		case(2)
    			!x+
    			dpxp=2.d0*(psi%ptrf_%f_(i+1,j,k)-psi%ptrf_%f_(i,j,k))
    			dpxp=dpxp-(psi%ptrf_%ptrf_%f_(i+1,j,k)-psi%ptrf_%ptrf_%f_(i,j,k))
    			dpxp=dpxp/mesh%dxc_(1)
    			dpxp=dpxp-2.d0*st%ux_%ptrf_%f_(i,j,k)+st%ux_%ptrf_%ptrf_%f_(i,j,k)
    			!x-
    			dpxm=2.d0*(psi%ptrf_%f_(i,j,k)-psi%ptrf_%f_(i-1,j,k))
    			dpxm=dpxm-(psi%ptrf_%ptrf_%f_(i,j,k)-psi%ptrf_%ptrf_%f_(i-1,j,k))
    			dpxm=dpxm/mesh%dxc_(1)
    			dpxm=dpxm-2.d0*st%ux_%ptrf_%f_(i-1,j,k)+st%ux_%ptrf_%ptrf_%f_(i-1,j,k)
    			!y+
    			dpyp=2.d0*(psi%ptrf_%f_(i,j+1,k)-psi%ptrf_%f_(i,j,k))
    			dpyp=dpyp-(psi%ptrf_%ptrf_%f_(i,j+1,k)-psi%ptrf_%ptrf_%f_(i,j,k))
    			dpyp=dpyp/mesh%dyc_(j+1)
    			dpyp=dpyp-2.d0*st%uy_%ptrf_%f_(i,j,k)+st%uy_%ptrf_%ptrf_%f_(i,j,k)   
    			!y-
    			dpym=2.d0*(psi%ptrf_%f_(i,j,k)-psi%ptrf_%f_(i,j-1,k))
    			dpym=dpym-(psi%ptrf_%ptrf_%f_(i,j,k)-psi%ptrf_%ptrf_%f_(i,j-1,k))
    			dpym=dpym/mesh%dyc_(j) 	
    			dpym=dpym-2.d0*st%uy_%ptrf_%f_(i,j-1,k)+st%uy_%ptrf_%ptrf_%f_(i,j-1,k)
    			!z+
    			dpzp=2.d0*(psi%ptrf_%f_(i,j,k+1)-psi%ptrf_%f_(i,j,k))
    			dpzp=dpzp-(psi%ptrf_%ptrf_%f_(i,j,k+1)-psi%ptrf_%ptrf_%f_(i,j,k))
    			dpzp=dpzp/mesh%dzc_(1)
    			dpzp=dpzp-2.d0*st%uz_%ptrf_%f_(i,j,k)+st%uz_%ptrf_%ptrf_%f_(i,j,k)
    			!z-
    			dpzm=2.d0*(psi%ptrf_%f_(i,j,k)-psi%ptrf_%f_(i,j,k-1))
    			dpzm=dpzm-(psi%ptrf_%ptrf_%f_(i,j,k)-psi%ptrf_%ptrf_%f_(i,j,k-1))
    			dpzm=dpzm/mesh%dzc_(1)
    			dpzm=dpzm-2.d0*st%uz_%ptrf_%f_(i,j,k-1)+st%uz_%ptrf_%ptrf_%f_(i,j,k-1)
    		case default 					
    	end select

    end subroutine
#endif
!========================================================================================!

!========================================================================================!
#ifdef MG_MODE
    subroutine computeBeta(this,rho)
        type(poissonEqn), intent(inout) :: this
        type(field), intent(in) :: rho
        integer :: nx, ny, nz
        real(DP) :: dt, alpha
        integer :: i,j,k
        
    	dt = this%ptrTime_%dt_
    	alpha = alphaRKS(this%ptrTime_)
    	
        nx = this%ptrMesh_%nx_
        ny = this%ptrMesh_%ny_
        nz = this%ptrMesh_%nz_
        
 		!$OMP PARALLEL DO DEFAULT(none) &
		!$OMP SHARED(this,rho) &
		!$OMP SHARED(dt,alpha) &
		!$OMP SHARED(nx,ny,nz) &
		!$OMP PRIVATE(i,j,k)
        do k = 1,nz
        	do j = 1,ny
        		do i = 1,nx
        			this%beta_%f_(i,j,k) = (dt*alpha)/rho%f_(i,j,k)
        		end do
        	end do
        end do
        !$OMP END PARALLEL DO 
         
        call updateBoundaries(this%beta_)
        
        
    end subroutine
#endif
!========================================================================================!

!========================================================================================!
    subroutine solvePoissonEqn(this,psi,rho,u,st)
        type(poissonEqn), intent(inout) :: this
        type(field), intent(inout) :: psi
        type(field), intent(in) :: rho
        type(vfield), intent(in) :: u,st
        real(DP) :: start, finish
        
        start = MPI_Wtime()   
                   
#ifdef FAST_MODE
		call computeSource(this,psi,rho,u,st)
		call solveFPS(this%fftSolver_,psi,this%s_) 
#endif
#ifdef MG_MODE
		call computeSource(this,u)
        call computeBeta(this,rho)
        call solvePCG(this%pcgs_,this%ptrMesh_,psi,this%beta_,this%s_)
#endif

        finish = MPI_Wtime()
        
        call info(this,finish-start)
    
    end subroutine
!========================================================================================!

!========================================================================================!
    subroutine updatePressure(this,psi,p,p0,st,st0)
    	type(poissonEqn), intent(in) :: this
        type(field), intent(inout) :: p,psi,p0
        type(vfield), intent(inout) :: st,st0
        integer :: nx, ny, nz
        integer :: i,j,k
        
        nx = p%ptrMesh_%nx_
        ny = p%ptrMesh_%ny_
        nz = p%ptrMesh_%nz_
        
 		!$OMP PARALLEL DO DEFAULT(none) &
		!$OMP SHARED(p,psi) &
		!$OMP SHARED(nx,ny,nz) &
		!$OMP PRIVATE(i,j,k)
        do k = 1,nz
        	do j = 1,ny
        		do i = 1,nx
        			p%f_(i,j,k) = psi%f_(i,j,k) !+ p%f_(i,j,k)
        		end do
        	end do
        end do
        !$OMP END PARALLEL DO 

        call updateBoundaries(p)
  
        !store old pressure fields
		call storeOldField(psi,this%nl_)
		p0%f_=psi%ptrf_%ptrf_%f_
		
		!store old surface tension fields
		call storeOldFieldV(st,this%nl_)
		st0%ux_%f_=st%ux_%ptrf_%ptrf_%f_
		st0%uy_%f_=st%uy_%ptrf_%ptrf_%f_
		st0%uz_%f_=st%uz_%ptrf_%ptrf_%f_
        
    end subroutine
!========================================================================================!

!========================================================================================!
    subroutine info(this,cpuTime)
    	type(poissonEqn), intent(in) :: this
    	real(DP), intent(in) :: cpuTime
    	real(DP) :: cpuTime_max
    	type(mpiControl), pointer :: comm
    	integer :: ierror
    	
    	comm => this%ptrMesh_%ptrMPIC_
    	
    	call Mpi_Reduce(cpuTime, cpuTime_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, &
    			        comm%cartComm_, ierror)

#ifdef FAST_MODE
		if (IS_MASTER) then
			write(*,'(A,'//s_outputFormat(2:9)//')') &
				'	P Eqn:	CPU time = ', cpuTime_max
		end if
#endif
#ifdef MG_MODE	
		if (IS_MASTER) then
			write(*,'(A,'//s_intFormat(2:3)//',2(A,'//s_outputFormat(2:9)//'))') &
				'	P Eqn: n iter = ', this%pcgs_%iter_, &
				'	residual = ', this%pcgs_%res_, '	CPU time = ', cpuTime_max
		end if
#endif  
			
	end subroutine
!========================================================================================!

end module poissonEqnMod

