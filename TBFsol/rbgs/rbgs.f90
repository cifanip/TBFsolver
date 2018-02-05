! ************************************************************************************** !
!    TBFsol - DNS turbulent bubbly flow solver
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

module rbgsMod
	
	use scalarFieldMod
	
	implicit none

	real(DP), parameter :: pi = 4.d0*DATAN(1.d0)

	type, public :: rbgs
	
		real(DP), private, allocatable, dimension(:,:) :: rList_, bList_ 
		real(DP), allocatable, dimension(:,:,:) :: resV_
		
		!store metrics
		real(DP), private, allocatable, dimension(:,:) :: invx_, invy_, invz_
		
		!class pointer for multi-grids
		type(rbgs), pointer :: ptrRbgs_ => NULL()
		
		type(grid), pointer, private :: ptrMesh_ => NULL()
		
		real(DP) :: res_
		
		logical :: isSystemSingular_
		

		contains
		
		final :: delete_rbgs


	end type
	
	
	private :: coarsenRbgsSolver
	private :: allocatePtrRbgs
	private :: deallocatePtrRbgs
	private :: checkSingularSystem
	private :: storeMetrics
	private :: iterateGS
	
	public :: rbgsCTOR
	public :: delete_rbgs
	public :: solveRBGS
	public :: computeResiduals
	public :: coarsenRbgsSolvers
	public :: correctFullNeumann
	

contains




!========================================================================================!
    subroutine delete_rbgs(this)
        type(rbgs), intent(inout) :: this
		
		call deallocatePtrRbgs(this) 
        
    end subroutine
!========================================================================================!

!========================================================================================!
	subroutine rbgsCTOR(this,mesh,p)
		type(rbgs) :: this
		type(grid), intent(in), target :: mesh
		type(scalarField), intent(in) :: p
		type(mpiControl), pointer :: ptrMPIC
		integer, dimension(3) :: procCoord
		integer :: nx, ny, nz
		integer :: i, j, k
		integer :: i0, j0, k0, i1, j1, k1
		integer :: nr, nb
		
		this%ptrMesh_ => mesh
		ptrMPIC => mesh%ptrMPIC_
		
		procCoord = ptrMPIC%procCoord_
		
		nx = mesh%nx_
		ny = mesh%ny_
		nz = mesh%nz_
		
		!info system
		call checkSingularSystem(this,p)
		
		!global indexes 
		i0 = mesh%i0g_
		i1 = mesh%i1g_
		j0 = mesh%j0g_
		j1 = mesh%j1g_
		k0 = mesh%k0g_
		k1 = mesh%k1g_
		
		
		!count first to allocate lists
		nr = 0
		nb = 0
		do k=k0,k1
			do j=j0,j1
				do i=i0,i1
					if ( mod(i+j+k,2) == 0 ) then
						nr = nr +1
					else
						nb = nb +1
					end if
				end do
			end do
		end do
	
		call allocateArray(this%rList_,1,3,1,nr)
		call allocateArray(this%bList_,1,3,1,nb)
		

		!set indices 
		nr = 1
		nb = 1
		do k=k0,k1
			do j=j0,j1
				do i=i0,i1
					if ( mod(i+j+k,2) == 0 ) then
						this%rList_(1,nr) = i-i0+1
						this%rList_(2,nr) = j-j0+1
						this%rList_(3,nr) = k-k0+1
						nr = nr+1
					else
						this%bList_(1,nb) = i-i0+1
						this%bList_(2,nb) = j-j0+1
						this%bList_(3,nb) = k-k0+1
						nb = nb+1
					end if
				end do
			end do
		end do
		
		
		!allocate residuals array
		call allocateArray(this%resV_,1,nx,1,ny,1,nz)
		
		!store metrics
		call storeMetrics(this)
		
		
			    
	end subroutine
!========================================================================================!

!========================================================================================!
	subroutine solveRBGS(this,p,beta,q,nIter,isToBeReset) 
        type(rbgs), intent(inout) :: this
        type(scalarField), intent(in) :: beta, q
        integer, intent(in) :: nIter
        type(scalarField), intent(inout) :: p
        logical, intent(in) :: isToBeReset
        
        !GS sweeps
        call iterateGS(this,p,beta,q,nIter,isToBeReset)

 		!update residuals
		call computeResiduals(this,p,beta,q) 
	

    end subroutine
!========================================================================================!

!========================================================================================!
    subroutine iterateGS(this,p,beta,q,nIter,isToBeReset) 
        type(rbgs), intent(inout) :: this
        type(scalarField), intent(in) :: beta, q
        integer, intent(in) :: nIter
        type(scalarField), intent(inout) :: p
        logical, intent(in) :: isToBeReset
        integer :: i, j, k, n
        real(DP) :: aR, aL, aT, aBo, aF, aBa
        integer :: iter

        
        if (isToBeReset) then
        	call initToZero(p)
        end if
        
  
        do iter = 1,nIter
        
        	!update red-list
			!$OMP PARALLEL DO DEFAULT(none) &
			!$OMP SHARED(this,p,beta,q) &
			!$OMP PRIVATE(i,j,k,n) &
			!$OMP PRIVATE(aR,aL,aT,aBo,aF,aBa)
			do n=1,size(this%rList_,2)
				i = this%rList_(1,n)
				j = this%rList_(2,n)
				k = this%rList_(3,n)
			

				aR = 0.5d0*(beta%f_(i+1,j,k)+beta%f_(i,j,k))*this%invx_(1,i)
				aL = 0.5d0*(beta%f_(i,j,k)+beta%f_(i-1,j,k))*this%invx_(2,i)
				aT = 0.5d0*(beta%f_(i,j+1,k)+beta%f_(i,j,k))*this%invy_(1,j)
				aBo = 0.5d0*(beta%f_(i,j,k)+beta%f_(i,j-1,k))*this%invy_(2,j)
				aF = 0.5d0*(beta%f_(i,j,k+1)+beta%f_(i,j,k))*this%invz_(1,k)
				aBa = 0.5d0*(beta%f_(i,j,k)+beta%f_(i,j,k-1))*this%invz_(2,k)

			
				p%f_(i,j,k) = ( &
						  	- q%f_(i,j,k) &
					      	+ (aR*p%f_(i+1,j,k)+aL*p%f_(i-1,j,k)) &
					      	+ (aT*p%f_(i,j+1,k)+aBo*p%f_(i,j-1,k)) &
					      	+ (aF*p%f_(i,j,k+1)+aBa*p%f_(i,j,k-1)) &
					      	  )/(aR+aL+aT+aBo+aF+aBa)
			end do
			!$OMP END PARALLEL DO
		
			!update boundaries
			call updateBoundaries(p)

			!update black-list
			!$OMP PARALLEL DO DEFAULT(none) &
			!$OMP SHARED(this,p,beta,q) &
			!$OMP PRIVATE(i,j,k,n) &
			!$OMP PRIVATE(aR,aL,aT,aBo,aF,aBa)
			do n=1,size(this%bList_,2)
				i = this%bList_(1,n)
				j = this%bList_(2,n)
				k = this%bList_(3,n)
			
				aR = 0.5d0*(beta%f_(i+1,j,k)+beta%f_(i,j,k))*this%invx_(1,i)
				aL = 0.5d0*(beta%f_(i,j,k)+beta%f_(i-1,j,k))*this%invx_(2,i)
				aT = 0.5d0*(beta%f_(i,j+1,k)+beta%f_(i,j,k))*this%invy_(1,j)
				aBo = 0.5d0*(beta%f_(i,j,k)+beta%f_(i,j-1,k))*this%invy_(2,j)
				aF = 0.5d0*(beta%f_(i,j,k+1)+beta%f_(i,j,k))*this%invz_(1,k)
				aBa = 0.5d0*(beta%f_(i,j,k)+beta%f_(i,j,k-1))*this%invz_(2,k)

			
				p%f_(i,j,k) = ( &
						  	- q%f_(i,j,k) &
					      	+ (aR*p%f_(i+1,j,k)+aL*p%f_(i-1,j,k)) &
					      	+ (aT*p%f_(i,j+1,k)+aBo*p%f_(i,j-1,k)) &
					      	+ (aF*p%f_(i,j,k+1)+aBa*p%f_(i,j,k-1)) &
					      	  )/(aR+aL+aT+aBo+aF+aBa)
			end do
			!$OMP END PARALLEL DO
		
			!update boundaries
			call updateBoundaries(p)   
			
			!correct full Neumann system 
			if (this%isSystemSingular_) then
				call correctFullNeumann(this%ptrMesh_,p)
			end if	
        	
        end do		
        
    end subroutine
!========================================================================================!

!========================================================================================!
    subroutine computeResiduals(this,p,beta,q)
        type(rbgs), intent(inout) :: this
        type(scalarField), intent(in) :: beta, q
        type(scalarField), intent(in) :: p
        type(mpiControl), pointer :: ptrMPIC
        integer :: i, j, k
        real(DP) :: aR, aL, aT, aBo, aF, aBa
        real(DP) :: res_l,resq_l,res_g,resq_g
        integer :: ierror
        integer :: nx, ny, nz
        
        nx = this%ptrMesh_%nx_
        ny = this%ptrMesh_%ny_
        nz = this%ptrMesh_%nz_
        
		!$OMP PARALLEL DO DEFAULT(none) &
		!$OMP SHARED(this,p,beta,q) &
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

					
					this%resV_(i,j,k) = q%f_(i,j,k) 							&
										-										&
										(										&
											  aR*(p%f_(i+1,j,k)-p%f_(i,j,k))	&
									  		- aL*(p%f_(i,j,k)-p%f_(i-1,j,k))    &
									  		+ aT*(p%f_(i,j+1,k)-p%f_(i,j,k))    &
									  		- aBo*(p%f_(i,j,k)-p%f_(i,j-1,k))   &
									  		+ aF*(p%f_(i,j,k+1)-p%f_(i,j,k))    &
									  		- aBa*(p%f_(i,j,k)-p%f_(i,j,k-1))   &	
										)
				end do
			end do
		end do
		!$OMP END PARALLEL DO
		

		call reduceSqrSum_omp(this%resV_,res_l)
		call reduceSqrSum_omp(q%f_,resq_l)
		
		ptrMPIC => this%ptrMesh_%ptrMPIC_
		call Mpi_Allreduce(res_l, res_g, 1, MPI_DOUBLE_PRECISION, MPI_SUM, ptrMPIC%cartComm_, ierror)
		call Mpi_Allreduce(resq_l, resq_g, 1, MPI_DOUBLE_PRECISION, MPI_SUM, ptrMPIC%cartComm_, ierror)
		
		this%res_ = sqrt(res_g)/(sqrt(resq_g)+tiny(1.d0))

    end subroutine
!========================================================================================!

!========================================================================================!
    recursive subroutine coarsenRbgsSolvers(this,mesh,p,n)
        type(rbgs), intent(inout) :: this
        type(grid), intent(in) :: mesh
        type(scalarField), intent(in) :: p
        integer, intent(in) :: n
		
		if (n > 1) then 
			call coarsenRbgsSolver(this,mesh,p)
			call coarsenRbgsSolvers(this%ptrRbgs_,mesh%ptrGrid_,p%ptrf_,n-1)
		else
			return
		end if
		
        
    end subroutine
!========================================================================================!

!========================================================================================!
    subroutine coarsenRbgsSolver(this,mesh,p)
        type(rbgs), intent(inout) :: this
        type(grid), intent(in) :: mesh
        type(scalarField), intent(in) :: p
		
		!allocate pointer
		call allocatePtrRbgs(this) 
		
		!init coarse solver
		 call rbgsCTOR(this%ptrRbgs_,mesh%ptrGrid_,p%ptrf_)
		
        
    end subroutine
!========================================================================================!

!========================================================================================!
	subroutine allocatePtrRbgs(this) 
		type(rbgs), intent(inout) :: this 
		integer :: err
		
		
		if (.not. associated(this%ptrRbgs_)) then
			
			allocate(this%ptrRbgs_,STAT=err)

			if (err /= 0) then
				call mpiABORT('Allocation of ptrRbgs failed ') 
			end if
		else
			call mpiABORT('Attempt to allocate an already associated ptrRbgs ')
		end if

		
	end subroutine
!========================================================================================!

!========================================================================================!
	subroutine deallocatePtrRbgs(this) 
		type(rbgs), intent(inout) :: this 

        if (associated(this%ptrRbgs_)) then
            deallocate(this%ptrRbgs_)
        end if
	
	end subroutine
!========================================================================================!

!========================================================================================!
    subroutine correctFullNeumann(mesh,p)
        type(grid), intent(in) :: mesh
        type(scalarField), intent(inout) :: p
        type(mpiControl), pointer :: ptrMPIC
        real(DP) :: pAvl, pAvg
        integer :: ierror
        integer :: is,ie,js,je,ks,ke
        integer :: nx, ny, nz
        integer :: i,j,k
        
        nx = mesh%nx_
        ny = mesh%ny_
        nz = mesh%nz_
        
		is = lbound(p%f_,1)
		ie = ubound(p%f_,1)
		js = lbound(p%f_,2)
		je = ubound(p%f_,2)
		ks = lbound(p%f_,3)
		ke = ubound(p%f_,3)
		
		call reduceSum_omp(p%f_,1,nx,1,ny,1,nz,pAvl)
		
		ptrMPIC => mesh%ptrMPIC_
		
		call Mpi_Allreduce(pAvl, pAvg, 1, MPI_DOUBLE_PRECISION, MPI_SUM, ptrMPIC%cartComm_, ierror)
		pAvg = pAvg/(mesh%nxg_*mesh%nyg_*mesh%nzg_)
		
		!make field zero-average
 		!$OMP PARALLEL DO DEFAULT(none) &
		!$OMP SHARED(p,pAvg) &
		!$OMP SHARED(is,ie,js,je,ks,ke) &
		!$OMP PRIVATE(i,j,k)
		do k=ks,ke
			do j=js,je
				do i=is,ie
					p%f_(i,j,k) = p%f_(i,j,k) - pAvg
				end do
			end do
		end do
		!$OMP END PARALLEL DO
		
		      
    end subroutine
!========================================================================================!

!========================================================================================!
	subroutine checkSingularSystem(this,f) 
		type(rbgs), intent(inout) :: this
		type(scalarField), intent(in) :: f
		type(mpiControl), pointer :: ptrMPIC
		integer :: i, ig
		integer :: ierror
		
		
		ptrMPIC => this%ptrMesh_%ptrMPIC_

		if (											  &
			(f%bLeft_%bType_   == s_fixedValue) .OR. &
			(f%bRight_%bType_  == s_fixedValue) .OR. &
			(f%bBottom_%bType_ == s_fixedValue) .OR. &
			(f%bTop_%bType_    == s_fixedValue) .OR. &
			(f%bBack_%bType_   == s_fixedValue) .OR. &
			(f%bFront_%bType_  == s_fixedValue)       &
			) then
			
			i = 1
			
		else
		
			i = 0
				
		end if
		
		call Mpi_Allreduce(i, ig, 1, MPI_INTEGER, MPI_SUM, ptrMPIC%cartComm_, ierror)
		
		if (ig == 0) then
			this%isSystemSingular_ = .TRUE.
		else
			this%isSystemSingular_ = .FALSE.
		end if
		
	end subroutine
!========================================================================================!

!========================================================================================!
	subroutine storeMetrics(this) 
        type(rbgs), intent(inout) :: this
        integer :: nx,ny,nz
        integer :: i,j,k
		
		nx = this%ptrMesh_%nx_
		ny = this%ptrMesh_%ny_
		nz = this%ptrMesh_%nz_
		
		call allocateArray(this%invx_,1,2,1,nx)
		call allocateArray(this%invy_,1,2,1,ny)
		call allocateArray(this%invz_,1,2,1,nz)
		
		!x
		do i=1,nx
			this%invx_(1,i) = 1.d0/(this%ptrMesh_%dxc_(i+1)*this%ptrMesh_%dxf_(i))
			this%invx_(2,i) = 1.d0/(this%ptrMesh_%dxc_(i)*this%ptrMesh_%dxf_(i))
		end do
		
		!y
		do j=1,ny
			this%invy_(1,j) = 1.d0/(this%ptrMesh_%dyc_(j+1)*this%ptrMesh_%dyf_(j))
			this%invy_(2,j) = 1.d0/(this%ptrMesh_%dyc_(j)*this%ptrMesh_%dyf_(j))
		end do
		
		!z
		do k=1,nz
			this%invz_(1,k) = 1.d0/(this%ptrMesh_%dzc_(k+1)*this%ptrMesh_%dzf_(k))
			this%invz_(2,k) = 1.d0/(this%ptrMesh_%dzc_(k)*this%ptrMesh_%dzf_(k))
		end do


    end subroutine
!========================================================================================!


end module rbgsMod

