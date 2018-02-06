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

!DIR$ IF DEFINED (FAST_MODE)	
	use fastPoissonSolverMod
!DIR$ ELSEIF DEFINED (MG_MODE)
	use pcgMod
!DIR$ ENDIF
	use timeMod
	
	implicit none


	type, public :: poissonEqn
	
		!keep a pointer to grid
		type(grid), pointer :: ptrMesh_ => NULL() 
		
!DIR$ IF DEFINED (FAST_MODE)		
		type(fastPoissonSolver), private :: fftSolver_
!DIR$ ELSEIF DEFINED (MG_MODE)
		type(pcg) :: pcgs_
!DIR$ ENDIF
		
		!poisson eqn source
		type(scalarField), private :: s_

!DIR$ IF DEFINED (MG_MODE)		
		!poisson eqn coefficients
		type(scalarField), private :: beta_
!DIR$ ENDIF
		
		!keep a pointer to time
		type(time), pointer :: ptrTime_ => NULL()

!DIR$ IF DEFINED (FAST_MODE)		
		!number of pressure levels
		integer :: nl_
		!reference density
		real(DP) :: rho0_ 
!DIR$ ENDIF

	end type
	
	
	private :: computeSource
!DIR$ IF DEFINED (MG_MODE)
	private :: computeBeta
!DIR$ ENDIF
	private :: info
!DIR$ IF DEFINED (FAST_MODE)	
	private :: computeOldPressDiv
	private :: allocateOldPressure
	private :: storeOldPressure
!DIR$ ENDIF

	public :: poissonEqnCTOR	
	public :: solvePoissonEqn
	public :: updatePressure


  	
contains


!========================================================================================!
!DIR$ IF DEFINED (FAST_MODE)
	subroutine poissonEqnCTOR(this,mesh,gMesh,psi,rt,rhol,rhog)
		type(poissonEqn) :: this
		type(grid), intent(in), target :: mesh,gMesh
		type(scalarField), intent(inout) :: psi
		type(time), intent(in), target :: rt
		real(DP), intent(in) :: rhol,rhog
		integer :: n_old
!DIR$ ELSEIF DEFINED (MG_MODE)
	subroutine poissonEqnCTOR(this,mesh,c,psi,rt)
		type(poissonEqn) :: this
		type(grid), intent(in), target :: mesh
		type(scalarField), intent(inout) :: c
		type(scalarField), intent(inout) :: psi
		type(time), intent(in), target :: rt
!DIR$ ENDIF
		
		this%ptrMesh_ => mesh
		this%ptrTime_ => rt

		!init source and beta coeff
		call scalarFieldCTOR(this%s_,'s',mesh,'cl',psi%hd_,initOpt=-1)
		call copyBoundary(this%s_,psi)
		
!DIR$ IF DEFINED (FAST_MODE)
		call fastPoissonSolverCTOR(this%fftSolver_,mesh,gMesh)
		
		!allocate space for old time pressure
		this%nl_=2
		call allocateOldPressure(psi,this%nl_)
		
		!set reference density
		this%rho0_=min(rhol,rhog)
		
!DIR$ ELSEIF DEFINED (MG_MODE)		
		call scalarFieldCTOR(this%beta_,'beta',mesh,'cl',psi%hd_,initOpt=-1)
		call copyBoundary(this%beta_,c)

		call pcgCTOR(this%pcgs_,mesh,psi,this%beta_)		
!DIR$ ENDIF		
		
			    
	end subroutine
!========================================================================================!

!========================================================================================!
!DIR$ IF DEFINED (MG_MODE)
    subroutine computeSource(this,u)
        type(poissonEqn), intent(inout) :: this
        type(vectorField), intent(in) :: u
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
    
!DIR$ ELSEIF DEFINED (FAST_MODE)
    subroutine computeSource(this,psi,rho,u)
        type(poissonEqn), intent(inout) :: this
        type(scalarField), intent(in) :: psi,rho
        type(vectorField), intent(in) :: u
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
		!$OMP SHARED(this,u,mesh,psi,rho) &
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
										     i,j,k,psi,mesh,this%nl_)
					
					dx = (apx*dpxp-amx*dpxm)/mesh%dxf_(1)
					dy = (apy*dpyp-amy*dpym)/mesh%dyf_(j)
					dz = (apz*dpzp-amz*dpzm)/mesh%dzf_(1)
					
					this%s_%f_(i,j,k)=this%s_%f_(i,j,k)+dx+dy+dz
        			

        		end do
        	end do
        end do
        !$OMP END PARALLEL DO 
        
    end subroutine
    
    subroutine computeOldPressDiv(dpxp,dpxm,dpyp,dpym,dpzp,dpzm,i,j,k,psi,mesh,n)
    	real(DP), intent(out) :: dpxp,dpxm,dpyp,dpym,dpzp,dpzm
    	integer, intent(in) :: i,j,k,n
    	type(scalarField), intent(in) :: psi
    	type(grid), intent(in) :: mesh
    	
    	select case(n)
    		case(2)
    			!x+
    			dpxp=2.d0*(psi%ptrf_%f_(i+1,j,k)-psi%ptrf_%f_(i,j,k))
    			dpxp=dpxp-(psi%ptrf_%ptrf_%f_(i+1,j,k)-psi%ptrf_%ptrf_%f_(i,j,k))
    			dpxp=dpxp/mesh%dxc_(1)
    			!x-
    			dpxm=2.d0*(psi%ptrf_%f_(i,j,k)-psi%ptrf_%f_(i-1,j,k))
    			dpxm=dpxm-(psi%ptrf_%ptrf_%f_(i,j,k)-psi%ptrf_%ptrf_%f_(i-1,j,k))
    			dpxm=dpxm/mesh%dxc_(1)
    			!y+
    			dpyp=2.d0*(psi%ptrf_%f_(i,j+1,k)-psi%ptrf_%f_(i,j,k))
    			dpyp=dpyp-(psi%ptrf_%ptrf_%f_(i,j+1,k)-psi%ptrf_%ptrf_%f_(i,j,k))
    			dpyp=dpyp/mesh%dyc_(j+1)    
    			!y-
    			dpym=2.d0*(psi%ptrf_%f_(i,j,k)-psi%ptrf_%f_(i,j-1,k))
    			dpym=dpym-(psi%ptrf_%ptrf_%f_(i,j,k)-psi%ptrf_%ptrf_%f_(i,j-1,k))
    			dpym=dpym/mesh%dyc_(j) 	
    			!z+
    			dpzp=2.d0*(psi%ptrf_%f_(i,j,k+1)-psi%ptrf_%f_(i,j,k))
    			dpzp=dpzp-(psi%ptrf_%ptrf_%f_(i,j,k+1)-psi%ptrf_%ptrf_%f_(i,j,k))
    			dpzp=dpzp/mesh%dzc_(1)
    			!z-
    			dpzm=2.d0*(psi%ptrf_%f_(i,j,k)-psi%ptrf_%f_(i,j,k-1))
    			dpzm=dpzm-(psi%ptrf_%ptrf_%f_(i,j,k)-psi%ptrf_%ptrf_%f_(i,j,k-1))
    			dpzm=dpzm/mesh%dzc_(1)
    		case default 					
    	end select

    end subroutine
!DIR$ ENDIF
!========================================================================================!

!========================================================================================!
!DIR$ IF DEFINED (MG_MODE)
    subroutine computeBeta(this,rho)
        type(poissonEqn), intent(inout) :: this
        type(scalarField), intent(in) :: rho
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
!DIR$ ENDIF
!========================================================================================!

!========================================================================================!
    subroutine solvePoissonEqn(this,psi,rho,u)
        type(poissonEqn), intent(inout) :: this
        type(scalarField), intent(inout) :: psi
        type(scalarField), intent(in) :: rho
        type(vectorField), intent(in) :: u
        real(DP) :: start, finish
        
        start = MPI_Wtime()   
                   
!DIR$ IF DEFINED (FAST_MODE)
		call computeSource(this,psi,rho,u)
		call solveFPS(this%fftSolver_,psi,this%s_) 
!DIR$ ELSEIF DEFINED (MG_MODE)
		call computeSource(this,u)
        call computeBeta(this,rho)
        call solvePCG(this%pcgs_,this%ptrMesh_,psi,this%beta_,this%s_)
!DIR$ ENDIF

        finish = MPI_Wtime()
        
        call info(this,finish-start)
    
    end subroutine
!========================================================================================!

!========================================================================================!
!DIR$ IF DEFINED (FAST_MODE)
    subroutine updatePressure(this,p,psi)
    	type(poissonEqn), intent(in) :: this
        type(scalarField), intent(inout) :: psi
!DIR$ ELSEIF DEFINED (MG_MODE)
    subroutine updatePressure(p,psi)
        type(scalarField), intent(in) :: psi
!DIR$ ENDIF
		type(scalarField), intent(inout) :: p
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
!DIR$ IF DEFINED (FAST_MODE)
        			p%f_(i,j,k) = psi%f_(i,j,k) !+ p%f_(i,j,k)
!DIR$ ELSEIF DEFINED (MG_MODE)
					p%f_(i,j,k) = psi%f_(i,j,k) !+ p%f_(i,j,k)
!DIR$ ENDIF
        		end do
        	end do
        end do
        !$OMP END PARALLEL DO 

        call updateBoundaries(p)

!DIR$ IF DEFINED (FAST_MODE)        
        !store old pressure fields
		call storeOldPressure(psi,this%nl_)
!DIR$ ENDIF
        
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

!DIR$ IF DEFINED (FAST_MODE)
		if (IS_MASTER) then
			write(*,'(A,'//s_outputFormat(2:9)//')') &
				'	P Eqn:	CPU time = ', cpuTime_max
		end if
!DIR$ ELSEIF DEFINED (MG_MODE)		
		if (IS_MASTER) then
			write(*,'(A,'//s_intFormat(2:3)//',2(A,'//s_outputFormat(2:9)//'))') &
				'	P Eqn: n iter = ', this%pcgs_%iter_, &
				'	residual = ', this%pcgs_%res_, '	CPU time = ', cpuTime_max
		end if
!DIR$ ENDIF 
			
	end subroutine
!========================================================================================!

!========================================================================================!
!DIR$ IF DEFINED (FAST_MODE)
    recursive subroutine allocateOldPressure(q,n)
		type(scalarField), intent(inout) :: q
		integer, intent(in) :: n
		integer :: lbi,ubi,lbj,ubj,lbk,ubk

		lbi=lbound(q%f_,1)
		ubi=ubound(q%f_,1)
		lbj=lbound(q%f_,2)
		ubj=ubound(q%f_,2)
		lbk=lbound(q%f_,3)
		ubk=ubound(q%f_,3)

		if (n>0) then
			call allocatePtrf(q)
			call allocateArray(q%ptrf_%f_,lbi,ubi,lbj,ubj,lbk,ubk)
			q%ptrf_%f_=q%f_
			call allocateOldPressure(q%ptrf_,n-1)
		else
			return
		end if
        
    end subroutine

    recursive subroutine storeOldPressure(q,n)
		type(scalarField), intent(inout) :: q
		integer, intent(in) :: n

		if (n>0) then
			if (associated(q%ptrf_)) then
				call storeOldPressure(q%ptrf_,n-1)
			end if
			call assign_omp(q%ptrf_%f_,q%f_)
		else
			return
		end if
        
    end subroutine
!DIR$ ENDIF 
!========================================================================================!


end module poissonEqnMod

