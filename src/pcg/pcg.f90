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

module pcgMod
	
	use multiGridMod
	
	implicit none


	type, public :: pcg
	
		type(field), private :: d_,r_
		real(DP), allocatable, dimension(:,:,:) :: qV_,d0_
		!store metrics
		real(DP), private, allocatable, dimension(:,:) :: invx_, invy_, invz_
		
		type(multiGrid), private :: mgs_
		
		real(DP), private :: tol_
		real(DP) :: res_
		integer :: iter_
		integer, private :: maxIter_
		logical, private :: fullInfo_
		logical, private :: isSystemSingular_
		

	end type
	
	private :: computeQv
	private :: computeAlpha
	private :: updatePsi
	private :: updateDir
	private :: updateDeltaN
	private :: updateResidual
	private :: continueIterating
	private :: scalarResidual
	private :: storeMetrics
	
	public :: pcgCTOR
	public :: solvePCG

contains


!========================================================================================!
	subroutine pcgCTOR(this,mesh,psi,beta)
		type(pcg) :: this
		type(grid), intent(in) ::  mesh
		type(field), intent(inout) :: psi,beta
		type(parFile) :: pfile
		integer :: nx,ny,nz
		
		nx = mesh%nx_
		ny = mesh%ny_
		nz = mesh%nz_

		!init field
		call fieldCTOR(this%d_,'d',mesh,'cl',psi%hd_,initOpt=-1)
		call copyBoundary(this%d_,psi)
		call resetBCerrorField(this%d_)
		call fieldCTOR(this%r_,'r',mesh,'cl',psi%hd_,initOpt=-1)
		call copyBoundary(this%r_,psi)
		call resetBCerrorField(this%r_)

		!read form pfile
		call parFileCTOR(pfile,'pcg_solver','specs')
		call readParameter(pfile,this%tol_,'tolPCG')
		call readParameter(pfile,this%maxIter_,'maxIterPCG')
		call readParameter(pfile,this%fullInfo_,'fullInfoPCG')
		
		!allocate auxiliary vectors
		call allocateArray(this%qV_,1,nx,1,ny,1,nz)
		!allocate auxiliary vectors
		call allocateArray(this%d0_,1,nx,1,ny,1,nz)
		
		!build multigrid
		call multiGridCTOR(this%mgs_,mesh,this%d_,beta,this%r_)
		
		this%isSystemSingular_ = this%mgs_%smoother_%isSystemSingular_
		
		call storeMetrics(this,mesh) 
					    
	end subroutine
!========================================================================================!

!========================================================================================!
	subroutine solvePCG(this,mesh,psi,beta,q)
		type(pcg), intent(inout) :: this
		type(grid), intent(in) :: mesh
		type(field), intent(inout) :: psi,beta,q
		real(DP) :: alpha, deltaN, delta0, gamma
		integer :: nx,ny,nz
		
		nx = mesh%nx_
		ny = mesh%ny_
		nz = mesh%nz_
		
		this%iter_ = 0
		
		call computeResiduals(this%mgs_%smoother_,psi,beta,q)
		call assign_omp(this%r_%f_,this%mgs_%smoother_%resV_,1,nx,1,ny,1,nz)
		call set2zero_omp(this%d_%f_)
		call solveMG(this%mgs_,mesh,this%d_,beta,this%r_)
		call updateDeltaN(this%r_,this%d_,deltaN)
		
		do while(continueIterating(this,q))
			call computeQv(this,beta,this%d_)
			call computeAlpha(this%d_,this%qv_,deltaN,alpha)
			call updatePsi(psi,this%d_,alpha,this%isSystemSingular_)
			call updateResidual(this%r_,this%qv_,alpha)		
			call assign_omp(this%d0_,this%d_%f_,1,nx,1,ny,1,nz)
			call set2zero_omp(this%d_%f_)
			call solveMG(this%mgs_,mesh,this%d_,beta,this%r_)
			delta0 = deltaN
			call updateDeltaN(this%r_,this%d_,deltaN)
			gamma = deltaN/delta0
			call updateDir(this%d0_,this%d_,gamma)
			call assign_omp(this%d_%f_,this%d0_,1,nx,1,ny,1,nz)
			call updateBoundaries(this%d_)
		end do
		
		call updateBoundaries(psi)
					    
	end subroutine
!========================================================================================!

!========================================================================================!
	subroutine computeQv(this,beta,d)
		type(pcg), intent(inout) :: this
		type(field), intent(in) :: beta,d
		integer :: i,j,k,nx,ny,nz
		real(DP) :: aR,aL,aT,aBo,aF,aBa
		
        nx = d%nx_
        ny = d%ny_
        nz = d%nz_

		!$OMP PARALLEL DO DEFAULT(none) &
		!$OMP SHARED(this,beta,d) &
		!$OMP SHARED(nx,ny,nz) &
		!$OMP PRIVATE(i,j,k) &
		!$OMP PRIVATE(aR,aL,aT,aBo,aF,aBa)          
		do k=1,nz
			do j=1,ny
				do i=1,nx

				aR = 0.5d0*(beta%f_(i+1,j,k)+beta%f_(i,j,k))*this%invx_(1,i)
				aL = 0.5d0*(beta%f_(i,j,k)+beta%f_(i-1,j,k))*this%invx_(2,i)
				aT = 0.5d0*(beta%f_(i,j+1,k)+beta%f_(i,j,k))*this%invy_(1,j)
				aBo = 0.5d0*(beta%f_(i,j,k)+beta%f_(i,j-1,k))*this%invy_(2,j)
				aF = 0.5d0*(beta%f_(i,j,k+1)+beta%f_(i,j,k))*this%invz_(1,k)
				aBa = 0.5d0*(beta%f_(i,j,k)+beta%f_(i,j,k-1))*this%invz_(2,k)

					
				this%qv_(i,j,k) =	  aR*(d%f_(i+1,j,k)-d%f_(i,j,k))	&
									- aL*(d%f_(i,j,k)-d%f_(i-1,j,k))    &
									+ aT*(d%f_(i,j+1,k)-d%f_(i,j,k))    &
									- aBo*(d%f_(i,j,k)-d%f_(i,j-1,k))   &
									+ aF*(d%f_(i,j,k+1)-d%f_(i,j,k))    &
									- aBa*(d%f_(i,j,k)-d%f_(i,j,k-1))   
				end do
			end do
		end do
		!$OMP END PARALLEL DO
		
					    
	end subroutine
!========================================================================================!

!========================================================================================!
	subroutine computeAlpha(d,qv,deltaN,alpha)
		type(field), intent(in) :: d
		real(DP), allocatable, dimension(:,:,:), intent(in) :: qv
		real(DP), intent(in) :: deltaN
		real(DP), intent(out) :: alpha
		integer :: nx,ny,nz
		real(DP) :: denoml,denomg
		integer :: ierror, comm
		integer :: i,j,k

        nx = d%nx_
       	ny = d%ny_
        nz = d%nz_
        
        denoml = 0.d0
        
		!$OMP PARALLEL DO DEFAULT(none) &
		!$OMP SHARED(d,qv) &
		!$OMP SHARED(nx,ny,nz) &
		!$OMP PRIVATE(i,j,k) &
		!$OMP REDUCTION(+:denoml)
		do k=1,nz
			do j=1,ny
				do i=1,nx
        			denoml = denoml + d%f_(i,j,k)*qv(i,j,k)
        		end do
        	end do
        end do
        !$OMP END PARALLEL DO
        
        
        comm = d%ptrMesh_%ptrMPIC_%cartComm_

        call Mpi_Allreduce(denoml, denomg, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierror)
        
		alpha = deltaN/denomg
						    
	end subroutine
!========================================================================================!

!========================================================================================!
	subroutine updatePsi(psi,d,alpha,isSingular)
		type(field), intent(inout) :: psi
		type(field), intent(in) :: d
		real(DP), intent(in) :: alpha
		logical, intent(in) :: isSingular
		integer :: is,ie,js,je,ks,ke
		integer :: i,j,k
		
		is = lbound(psi%f_,1)
		ie = ubound(psi%f_,1)
		js = lbound(psi%f_,2)
		je = ubound(psi%f_,2)
		ks = lbound(psi%f_,3)
		ke = ubound(psi%f_,3)
		
		!$OMP PARALLEL DO DEFAULT(none) &
		!$OMP SHARED(psi,d,alpha) &
		!$OMP SHARED(is,ie,js,je,ks,ke) &
		!$OMP PRIVATE(i,j,k) 
		do k=ks,ke
			do j=js,je
				do i=is,ie
					psi%f_(i,j,k) = psi%f_(i,j,k) + alpha*d%f_(i,j,k)
				end do
			end do
		end do
		!$OMP END PARALLEL DO
		
		!correct full Neumann system 
		if (isSingular) then
			!call correctFullNeumann(psi%ptrMesh_,psi)
		end if		
		

	end subroutine
!========================================================================================!

!========================================================================================!
	subroutine updateDir(d0,d,gamma)
		real(DP), allocatable, dimension(:,:,:), intent(inout) :: d0
		type(field), intent(in) :: d
		real(DP), intent(in) :: gamma
		integer :: nx,ny,nz
		integer :: i,j,k
		
		nx = d%nx_
		ny = d%ny_
		nz = d%nz_
		
		!$OMP PARALLEL DO DEFAULT(none) &
		!$OMP SHARED(d,d0) &
		!$OMP SHARED(nx,ny,nz,gamma) &
		!$OMP PRIVATE(i,j,k)			
		do k=1,nz
			do j=1,ny
				do i=1,nx
					d0(i,j,k) = d%f_(i,j,k) + gamma*d0(i,j,k)
				end do
			end do
		end do
		!$OMP END PARALLEL DO	
		
	end subroutine
!========================================================================================!

!========================================================================================!
	subroutine updateDeltaN(va,vb,deltaN)
		type(field), intent(in) :: va,vb
		real(DP), intent(out) :: deltaN
		integer :: nx,ny,nz
		real(DP) :: deltaNl
		integer :: ierror, comm
		integer :: i,j,k

        nx = va%nx_
        ny = va%ny_
        nz = va%nz_
        
        deltaNl = 0.d0
        
 		!$OMP PARALLEL DO DEFAULT(none) &
		!$OMP SHARED(va,vb) &
		!$OMP SHARED(nx,ny,nz) &
		!$OMP PRIVATE(i,j,k) &
		!$OMP REDUCTION(+:deltaNl)
		do k=1,nz
			do j=1,ny
				do i=1,nx
        			deltaNl = deltaNl + va%f_(i,j,k)*vb%f_(i,j,k)
        		end do
        	end do
        end do
        !$OMP END PARALLEL DO       
        
        comm = va%ptrMesh_%ptrMPIC_%cartComm_
              
        call Mpi_Allreduce(deltaNl, deltaN, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierror)

					    
	end subroutine
!========================================================================================!

!========================================================================================!
	subroutine updateResidual(r,qv,alpha)
		type(field), intent(inout) :: r
		real(DP), allocatable, dimension(:,:,:), intent(in) :: qv
		real(DP), intent(in) :: alpha
		integer :: nx,ny,nz
		integer :: i,j,k

        nx = r%nx_
        ny = r%ny_
        nz = r%nz_
        
 		!$OMP PARALLEL DO DEFAULT(none) &
		!$OMP SHARED(r,qv,alpha) &
		!$OMP SHARED(nx,ny,nz) &
		!$OMP PRIVATE(i,j,k)
		do k=1,nz
			do j=1,ny
				do i=1,nx
					r%f_(i,j,k) = r%f_(i,j,k)-alpha*qv(i,j,k)
        		end do
        	end do
        end do
		!$OMP END PARALLEL DO 
					    
	end subroutine
!========================================================================================!

!========================================================================================!
    function continueIterating(this,q) RESULT(isVar)
        type(pcg), intent(inout) :: this
        type(field), intent(in) :: q
        logical :: isVar
        
        
	!check max iteration limit
	if (this%iter_ == this%maxIter_) then
		if (IS_MASTER) then 
			write(*,*) '!*************** WARNING *****************'
			write(*,*) 'EXIT PCG iterations: max iter reached'
		end if
		isVar = .FALSE.
		return
	end if
	
	call scalarResidual(this%r_,q,this%res_)
	
	!check tolerance
	if (this%res_ > this%tol_) then
	
		if (IS_MASTER) then 
			if (this%fullInfo_) then
				write(*,*) 'PCG solver: iteration ', this%iter_, ' residual = ', &
				this%res_
			end if
		end if
	   isVar = .TRUE.
	   
	else
	
		if (IS_MASTER) then
			if (this%fullInfo_) then
				write(*,*) 'Criterion met at iteration: ', this%iter_, ' residual = ', &
				this%res_
			end if
		end if
		
		isVar = .FALSE.
		
	end if
	
	this%iter_ = this%iter_ + 1
		

    end function
!========================================================================================!

!========================================================================================!
	subroutine scalarResidual(r,q,res)
		type(field), intent(in) :: r,q
		real(DP), intent(out) :: res
		integer :: nx,ny,nz
		integer :: comm, ierror
		real(DP) :: r_l,q_l,r_g,q_g

		nx = q%nx_
		ny = q%ny_
		nz = q%nz_
		
		call reduceSqrSum_omp(r%f_,1,nx,1,ny,1,nz,r_l)
		call reduceSqrSum_omp(q%f_,1,nx,1,ny,1,nz,q_l)
		
		comm = r%ptrMesh_%ptrMPIC_%cartComm_
		
		call Mpi_Allreduce(r_l, r_g, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierror)
		call Mpi_Allreduce(q_l, q_g, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierror)
		
		res = sqrt(r_g)/sqrt(q_g)
					    
	end subroutine
!========================================================================================!

!========================================================================================!
	subroutine storeMetrics(this,mesh) 
        type(pcg), intent(inout) :: this
        type(grid), intent(in) :: mesh
        integer :: nx,ny,nz
        integer :: i,j,k
		
		nx = mesh%nx_
		ny = mesh%ny_
		nz = mesh%nz_
		
		call allocateArray(this%invx_,1,2,1,nx)
		call allocateArray(this%invy_,1,2,1,ny)
		call allocateArray(this%invz_,1,2,1,nz)
		
		!x
		do i=1,nx
			this%invx_(1,i) = 1.d0/(mesh%dxc_(i+1)*mesh%dxf_(i))
			this%invx_(2,i) = 1.d0/(mesh%dxc_(i)*mesh%dxf_(i))
		end do
		
		!y
		do j=1,ny
			this%invy_(1,j) = 1.d0/(mesh%dyc_(j+1)*mesh%dyf_(j))
			this%invy_(2,j) = 1.d0/(mesh%dyc_(j)*mesh%dyf_(j))
		end do
		
		!z
		do k=1,nz
			this%invz_(1,k) = 1.d0/(mesh%dzc_(k+1)*mesh%dzf_(k))
			this%invz_(2,k) = 1.d0/(mesh%dzc_(k)*mesh%dzf_(k))
		end do


    end subroutine
!========================================================================================!



end module pcgMod

