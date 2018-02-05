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

module poissMatMod
	
	use scalarFieldMod
	
	implicit none

	INCLUDE 'dmumps_struc.h'

	type, public :: poissMat
	
		type (DMUMPS_STRUC), private :: mumpsPar
		type(grid), pointer, private :: ptrMesh_ => null()
		
		integer, private :: entryCount_
		
		logical, private :: isSystemSingular_
		
		!local rhs to be sent to the master 
		real(DP), allocatable, private, dimension(:,:,:) :: rhs_loc
		
		contains
		
		final :: delete_poissMat

	end type
	
	private :: allocateMat
	private :: updateSystem
	private :: scatterSolution
	private :: assembleIndexesMat
	private :: computeInternalMatrixIndexes
	private :: computeLagrangianMulIndexes
	private :: computeBoundaryMatrixIndexes
	private :: computeDiagBoundaryIndexes
	private :: computeOffDiagBoundaryIndexes
	private :: updateMatrix
	private :: computeInternalMatrixValues
	private :: computeLagrangianMulValues
	private :: computeBoundaryMatrixValues
	private :: computeDiagBoundaryValues
	private :: computeOffDiagBoundaryValues
	private :: computeRHS
	private :: gatherRHS
	private :: computeBoundarySourceRHS
	private :: setEntryIndex
	private :: setEntryValue
	private :: diagMatIndex
	private :: setOrthBoundaryIndex
	private :: setIndexBounds
	private :: setDerivIndexes
	private :: setDirSwitch
	private :: setPatchSwitch
	private :: setOffSetLoopOverBoundaries
	private :: checkSingularSystem
	
	public :: delete_poissMat
	public :: poissMatCTOR
	public :: solveMat

  	
contains

	


!========================================================================================!
    subroutine delete_poissMat(this)
        type(poissMat), intent(inout) :: this
   
        if (this%mumpsPar%JOB .ne. 0)	then
        	this%mumpsPar%JOB = -2
        	call DMUMPS(this%mumpsPar)
        end if
        
    end subroutine
!========================================================================================!

!========================================================================================!
	subroutine poissMatCTOR(this,mesh,f)
		type(poissMat) :: this
		type(grid), intent(in), target :: mesh
		type(scalarField), intent(in) :: f
		type(mpiControl), pointer :: ptrMPIC
		
		!keep a pointer to mesh
		this%ptrMesh_ => mesh
		
		!info system
		call checkSingularSystem(this,f)
		
		ptrMPIC => mesh%ptrMPIC_
		
		!set MUMPS parameters
		!communicator
		this%mumpsPar%comm = ptrMPIC%cartComm_
		!init job
		this%mumpsPar%JOB = -1
		if (this%isSystemSingular_) then
			!general symmetric
			this%mumpsPar%SYM = 1
		else
			!symmetric and positive definite matrix
			this%mumpsPar%SYM = 1
		end if
		!the master-node takes part to the solution
		this%mumpsPar%PAR = 1
		!init mumps
		call DMUMPS(this%mumpsPar)
		
		!distributed matrix
		this%mumpsPar%ICNTL(18) = 3
		!centralised solution
		this%mumpsPar%ICNTL(21) = 0
		
		!print out error messages only
		this%mumpsPar%ICNTL(1) = 1
		this%mumpsPar%ICNTL(2) = 0
		this%mumpsPar%ICNTL(3) = 0
		this%mumpsPar%ICNTL(4) = 1
		
		!allocate linear system
		call allocateMat(this,f)
		
		!compute matrix indexes
		call assembleIndexesMat(this,f) 
		
		!analysis phase
		this%mumpsPar%job = 1
		call DMUMPS(this%mumpsPar)
		    
	end subroutine
!========================================================================================!


!========================================================================================!
	subroutine allocateMat(this,f) 
		type(poissMat), intent(inout) :: this
		type(scalarField), intent(in) :: f
		integer :: nx, ny, nz
		integer :: nxg, nyg, nzg
		integer :: nInter, nBound, nEntries
		
		!local mesh size
		nx = this%ptrMesh_%nx_
		ny = this%ptrMesh_%ny_
		nz = this%ptrMesh_%nz_
		!global mesh size
		nxg = this%ptrMesh_%nxg_
		nyg = this%ptrMesh_%nyg_
		nzg = this%ptrMesh_%nzg_
		
		!number of internal matrix entries
		nInter = (nx-2)*(ny-2)*(nz-2)*4

		if ( ((nx-2) < 0 ) .OR. ((ny-2) < 0) .OR. ((nz-2) < 0 ) ) then
			call mpiABORT('Attempt to init poissMat with too coarse mesh ')
		end if
		
		!number of boundary matrix entries (off diagonal)
		
		!along x
		!bLeft: if the patch is not periodic, even if it is internal has an entry 
		!       in the lower diagonal, does not count
		if (f%bLeft_%bType_ == s_periodicBC) then 
			nBound = 2*(f%bLeft_%n1_*f%bLeft_%n2_)
		else
			nBound = f%bLeft_%n1_*f%bLeft_%n2_
		end if
		!bRight
		if ( f%bRight_%bType_ == s_parallelBC ) then
			nBound = nBound + f%bRight_%n1_*f%bRight_%n2_
		end if
		!bTop + bBottom
		nBound = nBound + 2*(f%bTop_%n1_-2)*f%bTop_%n2_
		!bFront + bBack
		nBound = nBound + 2*(f%bFront_%n1_-2)*(f%bFront_%n2_-2)
		
		!along y
		if (f%bBottom_%bType_ == s_periodicBC) then 
			nBound = nBound + 2*(f%bBottom_%n1_*f%bBottom_%n2_)
		else
			nBound = nBound + f%bBottom_%n1_*f%bBottom_%n2_
		end if
		!bTop
		if ( f%bTop_%bType_ == s_parallelBC ) then
			nBound = nBound + f%bTop_%n1_*f%bTop_%n2_
		end if
		!bLeft + bRight
		nBound = nBound + 2*(f%bRight_%n1_-2)*f%bRight_%n2_
		!bFront + bBack
		nBound = nBound + 2*(f%bFront_%n1_-2)*(f%bFront_%n2_-2)
				
		!along z
		if (f%bBack_%bType_ == s_periodicBC) then
			nBound = nBound + 2*(f%bBack_%n1_*f%bBack_%n2_)
		else
			nBound = nBound + f%bBack_%n1_*f%bBack_%n2_
		end if
		!bFront
		if ( f%bFront_%bType_ == s_parallelBC ) then
			nBound = nBound + f%bFront_%n1_*f%bFront_%n2_
		end if
		!bLeft + bRight
		nBound = nBound + 2*f%bRight_%n1_*(f%bRight_%n2_-2)
		!bTop + bBottom
		nBound = nBound + 2*(f%bTop_%n1_-2)*(f%bTop_%n2_-2)
		
		!add up central boundary entries (!note: doubled entries are summed)
		!along x
		nBound = nBound + 2*(f%bRight_%n1_*f%bRight_%n2_) 		  &
						+ 2*((f%bTop_%n1_-2)*f%bTop_%n2_)			  &
						+ 2*((f%bFront_%n1_-2)*(f%bFront_%n2_-2))		
		
		!along y				
		nBound = nBound + 2*(f%bTop_%n1_*f%bTop_%n2_) 		  	&
						+ 2*((f%bRight_%n1_-2)*f%bRight_%n2_)	 	&
						+ 2*((f%bFront_%n1_-2)*(f%bFront_%n2_-2))	
						
		!along z				
		nBound = nBound + 2*(f%bFront_%n1_*f%bFront_%n2_) 		  	&
						+ 2*(f%bRight_%n1_*(f%bRight_%n2_-2))	 		&
						+ 2*((f%bTop_%n1_-2)*(f%bTop_%n2_-2))	
		
		!total entries
		nEntries = nInter + nBound
		
		!correct full Neumann bc matrix
		if (this%isSystemSingular_) then
			nEntries = nEntries + nx*ny*nz
		end if
		
		!allocate space for mumps RHS
		if (IS_MASTER) then
			if (this%isSystemSingular_) then
				this%mumpsPar%N = nxg*nyg*nzg + 1
			else
				this%mumpsPar%N = nxg*nyg*nzg
			end if
			allocate(this%mumpsPar%rhs(this%mumpsPar%N))
			this%mumpsPar%rhs = 0.d0
		end if
		
		!allocate space for rhs_loc
		call allocateArray(this%rhs_loc,1-f%hd_,nx+f%hd_,&
						   1-f%hd_,ny+f%hd_,1-f%hd_,nz+f%hd_)
				
		!allocate mumps_vectors
		this%mumpsPar%nz_loc = nEntries
		allocate(this%mumpsPar%irn_loc(nEntries))
		allocate(this%mumpsPar%jcn_loc(nEntries))
		allocate(this%mumpsPar%A_loc(nEntries))

		
		
	end subroutine
!========================================================================================!


!========================================================================================!
	subroutine updateSystem(this,f,beta,s) 
		type(poissMat), intent(inout) :: this
		type(scalarField), intent(in) :: f, beta, s
		
		call updateMatrix(this,f,beta)
		call computeRHS(this,f,beta,s)

	end subroutine
!========================================================================================!

!========================================================================================!
    subroutine solveMat(this,f,beta,s)
        type(poissMat), intent(inout) :: this
		type(scalarField), intent(inout) :: f
		type(scalarField), intent(in) :: beta, s

        
        call updateSystem(this,f,beta,s) 
        
        !factorisation phase
        this%mumpsPar%job = 2
		call DMUMPS(this%mumpsPar)
		
        !solution phase
        this%mumpsPar%job = 3
		call DMUMPS(this%mumpsPar)
		
		!assign solution vector to solution field
		call scatterSolution(this,f)
		
		!update boundaries
		call updateBoundaries(f)
        
        
    end subroutine
!========================================================================================!

!========================================================================================!
   subroutine scatterSolution(this,f)
        type(poissMat), intent(inout) :: this
		type(scalarField), intent(inout) :: f
        integer :: nx, ny, nz
        integer :: nxg, nyg, nzg
        integer :: i, j, k, n, q, np
        integer :: i0g, j0g, k0g, i1g, j1g, k1g
        type(mpiControl), pointer :: ptrMPIC
        integer :: ierror
        integer :: tag = 0
        integer, dimension(MPI_STATUS_SIZE) :: status


        ptrMPIC => this%ptrMesh_%ptrMPIC_
        np = ptrMPIC%nProcs_
        
        nxg = this%ptrMesh_%nxg_
        nyg = this%ptrMesh_%nyg_
        nzg = this%ptrMesh_%nzg_
        
        nx = this%ptrMesh_%nx_
        ny = this%ptrMesh_%ny_
        nz = this%ptrMesh_%nz_
        
        
        !update solution in the master proc first
		if (IS_MASTER) then
		
			i0g = this%ptrMesh_%i0g_
			j0g = this%ptrMesh_%j0g_
			k0g = this%ptrMesh_%k0g_
			i1g = this%ptrMesh_%i1g_
			j1g = this%ptrMesh_%j1g_
			k1g = this%ptrMesh_%k1g_
				
			!loop over grid points
			do k=k0g,k1g
				do j=j0g,j1g
					do i=i0g,i1g
							
						call diagMatIndex(n,i,j,k,nxg,nyg)
						
						f%f_(i-i0g+1,j-j0g+1,k-k0g+1) = this%mumpsPar%rhs(n)
							
					end do
				end do
			end do		
		
		end if
		
        !send solution to slaves
		if (IS_MASTER) then
			!loop over cores
			do q=1,np-1
			
				!set global indexes proc q
				call globalIndexesFromAll(this%ptrMesh_,i0g,i1g,j0g,j1g,k0g,k1g,q,f%tp_) 
				
				!loop over grid points
				do k=k0g,k1g
					do j=j0g,j1g
						do i=i0g,i1g
							
							call diagMatIndex(n,i,j,k,nxg,nyg)
							this%rhs_loc(i-i0g+1,j-j0g+1,k-k0g+1) = this%mumpsPar%rhs(n)
							
						end do
					end do
				end do
				
				call MPI_Send(this%rhs_loc,size(this%rhs_loc),MPI_DOUBLE_PRECISION,q,tag,ptrMPIC%cartComm_,ierror)
			end do
		
		else
			call MPI_Recv(f%f_,size(f%f_),MPI_DOUBLE_PRECISION,0,tag,ptrMPIC%cartComm_,status,ierror)
		end if
        
        
    end subroutine
!========================================================================================!


!****************************************************************************************!
!************************************ Indexes *******************************************!
!****************************************************************************************!

!========================================================================================!
	subroutine assembleIndexesMat(this,f) 
		type(poissMat), intent(inout) :: this
		type(scalarField), intent(in) :: f
		
		!init entry counter
		this%entryCount_ = 0

		call computeInternalMatrixIndexes(this)
		call computeBoundaryMatrixIndexes(this,f)
		
		!correct full Neumann bc
		if (this%isSystemSingular_) then
			call computeLagrangianMulIndexes(this)
		end if
		
	end subroutine
!========================================================================================!

!========================================================================================!
	subroutine computeInternalMatrixIndexes(this) 
		type(poissMat), intent(inout) :: this
		integer :: i, j, k, n
		integer :: i0g, j0g, k0g    !global indexes start
		integer :: i1g, j1g, k1g    !global indexes end
		integer :: nxg, nyg, nzg 	!global size
		
		nxg = this%ptrMesh_%nxg_
		nyg = this%ptrMesh_%nyg_
		nzg = this%ptrMesh_%nzg_
		
		i0g = this%ptrMesh_%i0g_
		j0g = this%ptrMesh_%j0g_
		k0g = this%ptrMesh_%k0g_
		i1g = this%ptrMesh_%i1g_
		j1g = this%ptrMesh_%j1g_
		k1g = this%ptrMesh_%k1g_
		
		!nodi interni
		do k = k0g+1, k1g-1
			do j = j0g+1, j1g-1
				do i = i0g+1, i1g-1
				
					!compute diag matrix index
					call diagMatIndex(n,i,j,k,nxg,nyg)

					!matrix indexes
					!node (i+1,j,k)
					call setEntryIndex(this,n,n+1) 
					!node (i,j+1,k)
					call setEntryIndex(this,n,n+nxg) 
					!node (i,j,k+1)
					call setEntryIndex(this,n,n+nxg*nyg)
					!node (i,j,k)
					call setEntryIndex(this,n,n)
					
				end do
			end do
		end	do
		
		
	end subroutine
!========================================================================================!

!========================================================================================!
	subroutine computeLagrangianMulIndexes(this) 
		type(poissMat), intent(inout) :: this
		integer :: i, j, k, n
		integer :: i0g, j0g, k0g    !global indexes start
		integer :: i1g, j1g, k1g    !global indexes end
		integer :: nxg, nyg, nzg 	!global size
		integer :: d
		
		nxg = this%ptrMesh_%nxg_
		nyg = this%ptrMesh_%nyg_
		nzg = this%ptrMesh_%nzg_
		
		i0g = this%ptrMesh_%i0g_
		j0g = this%ptrMesh_%j0g_
		k0g = this%ptrMesh_%k0g_
		i1g = this%ptrMesh_%i1g_
		j1g = this%ptrMesh_%j1g_
		k1g = this%ptrMesh_%k1g_
		
		d = nxg*nyg*nzg + 1
		
		!nodi interni
		do k = k0g, k1g
			do j = j0g, j1g
				do i = i0g, i1g
				
					!compute diag matrix index
					call diagMatIndex(n,i,j,k,nxg,nyg)
					
					call setEntryIndex(this,n,d) 
					
				end do
			end do
		end	do
		
		
	end subroutine
!========================================================================================!

!========================================================================================!
	subroutine computeBoundaryMatrixIndexes(this,f) 
		type(poissMat), intent(inout) :: this
		type(scalarField), intent(in) :: f
		integer :: i0g, j0g, k0g    !global indexes start
		integer :: i1g, j1g, k1g    !global indexes end
		integer :: nxg, nyg, nzg 	!global size
		integer :: q
		
		!assign sizes
		nxg = this%ptrMesh_%nxg_
		nyg = this%ptrMesh_%nyg_
		nzg = this%ptrMesh_%nzg_
		
		i0g = this%ptrMesh_%i0g_
		j0g = this%ptrMesh_%j0g_
		k0g = this%ptrMesh_%k0g_
		i1g = this%ptrMesh_%i1g_
		j1g = this%ptrMesh_%j1g_
		k1g = this%ptrMesh_%k1g_
		
		!off diagonal boundary indexes
		do q=1,3
			call computeOffDiagBoundaryIndexes(this,f%bLeft_,q,i0g,i1g,j0g,j1g,k0g,k1g,nxg,nyg,nzg)
			call computeOffDiagBoundaryIndexes(this,f%bRight_,q,i0g,i1g,j0g,j1g,k0g,k1g,nxg,nyg,nzg)
			call computeOffDiagBoundaryIndexes(this,f%bBottom_,q,i0g,i1g,j0g,j1g,k0g,k1g,nxg,nyg,nzg)
			call computeOffDiagBoundaryIndexes(this,f%bTop_,q,i0g,i1g,j0g,j1g,k0g,k1g,nxg,nyg,nzg)
			call computeOffDiagBoundaryIndexes(this,f%bBack_,q,i0g,i1g,j0g,j1g,k0g,k1g,nxg,nyg,nzg)
			call computeOffDiagBoundaryIndexes(this,f%bFront_,q,i0g,i1g,j0g,j1g,k0g,k1g,nxg,nyg,nzg)
		end do
		
		!central boundary indexes
		do q=1,3
			call computeDiagBoundaryIndexes(this,f%bLeft_,q,i0g,i1g,j0g,j1g,k0g,k1g,nxg,nyg)
			call computeDiagBoundaryIndexes(this,f%bRight_,q,i0g,i1g,j0g,j1g,k0g,k1g,nxg,nyg)
			call computeDiagBoundaryIndexes(this,f%bBottom_,q,i0g,i1g,j0g,j1g,k0g,k1g,nxg,nyg)
			call computeDiagBoundaryIndexes(this,f%bTop_,q,i0g,i1g,j0g,j1g,k0g,k1g,nxg,nyg)
			call computeDiagBoundaryIndexes(this,f%bBack_,q,i0g,i1g,j0g,j1g,k0g,k1g,nxg,nyg)
			call computeDiagBoundaryIndexes(this,f%bFront_,q,i0g,i1g,j0g,j1g,k0g,k1g,nxg,nyg)
		end do

		
	end subroutine
!========================================================================================!


!========================================================================================!
	subroutine computeDiagBoundaryIndexes(this,bf,dir,i0g,i1g,j0g,j1g,k0g,k1g,nxg,nyg) 
		type(poissMat), intent(inout) :: this
		type(scalarBoundaryField), intent(in) :: bf
		integer, intent(in) :: dir
		integer, intent(in) :: i0g, j0g, k0g    !global indexes start
		integer, intent(in) :: i1g, j1g, k1g    !global indexes end
		integer, intent(in) :: nxg, nyg 	    !global size
		integer :: i, j, k, n
		integer :: di, dj, dk
		integer :: is, ie, js, je, ks, ke
		integer :: patchDir

		patchDir = bf%cartDir_
	
		!set indexes bounds 
		call setIndexBounds(patchDir,is,ie,js,je,ks,ke,i0g,i1g,j0g,j1g,k0g,k1g)
		!set off-set from mesh size
		call setOffSetLoopOverBoundaries(patchDir,dir,di,dj,dk)
		

		do k=ks+dk,ke-dk
			do j=js+dj,je-dj
				do i=is+di,ie-di
					call diagMatIndex(n,i,j,k,nxg,nyg)
					call setEntryIndex(this,n,n)
				end do
			end do
		end do	
		
	end subroutine
!========================================================================================!


	
!========================================================================================!
	subroutine computeOffDiagBoundaryIndexes(this,bf,dir,i0g,i1g,j0g,j1g,k0g,k1g,nxg,nyg,nzg) 
		type(poissMat), intent(inout) :: this
		type(scalarBoundaryField), intent(in) :: bf
		integer, intent(in) :: dir
		integer, intent(in) :: i0g, j0g, k0g    !global indexes start
		integer, intent(in) :: i1g, j1g, k1g    !global indexes end
		integer, intent(in) :: nxg, nyg, nzg 	!global size
		integer :: i, j, k, n
		integer :: is, ie, js, je, ks, ke
		integer :: d, p, di, dj, dk
		integer :: patchDir, bType
	
	
		bType = bf%bType_
		patchDir = bf%cartDir_
		
		!set index bounds 
		call setIndexBounds(patchDir,is,ie,js,je,ks,ke,i0g,i1g,j0g,j1g,k0g,k1g)
		!set off-set from mesh size
		call setOffSetLoopOverBoundaries(patchDir,dir,di,dj,dk)
		
		!set matrix entry displacement
		SELECT CASE (dir)
   		CASE (1) 
			d = 1
			p = nxg-1
		CASE (2)
			d = nxg
			p = nxg*(nyg-1)	
   		CASE (3) 
			d = nxg*nyg
			p = nxg*nyg*(nzg-1)
   		CASE DEFAULT
		END SELECT
		
	
		if ( abs(patchDir) == dir) then
	
			SELECT CASE(patchDir)
			
				CASE(-1,-2,-3)
			
					!check periodic boundary
					if (bType == s_periodicBC) then 
						do k=ks,ke
							do j=js,je
								do i = is, ie
									call diagMatIndex(n,i,j,k,nxg,nyg)
									call setEntryIndex(this,n,n+p)
								end do
							end do
						end do
					end if 
					!node +d
					do k=ks,ke
						do j=js,je
							do i=is,ie
								call diagMatIndex(n,i,j,k,nxg,nyg)
								call setEntryIndex(this,n,n+d)
							end do
						end do
					end do
			
				CASE(1,2,3)
			
					if ( bType == s_parallelBC ) then 
						!node +d
						do k=ks,ke
							do j=js,je
								do i = is, ie
									call diagMatIndex(n,i,j,k,nxg,nyg)
									call setEntryIndex(this,n,n+d)
								end do
							end do
						end do
					end if	
	
	
				CASE DEFAULT
			
			END SELECT
		
		else
				
			!node (i,j,k)+d
			do k=ks+dk,ke-dk
				do j=js+dj,je-dj
					do i = is+di, ie-di
						call diagMatIndex(n,i,j,k,nxg,nyg)
						call setEntryIndex(this,n,n+d)
					end do
				end do
			end do
			
		end if
	
	end subroutine
!========================================================================================!		
	
	

!****************************************************************************************!
!************************************* Values *******************************************!
!****************************************************************************************!

!========================================================================================!
	subroutine updateMatrix(this,f,beta) 
		type(poissMat), intent(inout) :: this
		type(scalarField), intent(in) :: f
		type(scalarField), intent(in) :: beta
		
		!init entry counter
		this%entryCount_ = 0

		call computeInternalMatrixValues(this,beta)
		call computeBoundaryMatrixValues(this,f,beta)
		
		!correct full Neumann bc
		if (this%isSystemSingular_) then
			call computeLagrangianMulValues(this)
		end if
		
	end subroutine
!========================================================================================!

!========================================================================================!
	subroutine computeInternalMatrixValues(this,beta) 
		type(poissMat), intent(inout) :: this
		type(scalarField), intent(in) :: beta
		integer :: i, j, k
		integer :: il, jl, kl
		integer :: i0g, j0g, k0g, i1g, j1g, k1g	    !global indexes 
		integer :: nxg, nyg, nzg 					!global size
		real(DP) :: aL, aR, aBo, aT, aBa, aF
		real(DP) :: V
		
		!assign sizes	
		nxg = this%ptrMesh_%nxg_
		nyg = this%ptrMesh_%nyg_
		nzg = this%ptrMesh_%nzg_
		
		i0g = this%ptrMesh_%i0g_
		j0g = this%ptrMesh_%j0g_
		k0g = this%ptrMesh_%k0g_
		i1g = this%ptrMesh_%i1g_
		j1g = this%ptrMesh_%j1g_
		k1g = this%ptrMesh_%k1g_
		
		
		!internal nodes
		do k = k0g+1, k1g-1
			do j = j0g+1, j1g-1
				do i = i0g+1, i1g-1
				
					il = i-i0g+1
					jl = j-j0g+1
					kl = k-k0g+1
					
					V = this%ptrMesh_%V_(il,jl,kl)
					
					aR = 0.5d0*(beta%f_(il+1,jl,kl)+beta%f_(il,jl,kl))/(this%ptrMesh_%dxc_(il+1)*this%ptrMesh_%dxf_(il))
					aL = 0.5d0*(beta%f_(il,jl,kl)+beta%f_(il-1,jl,kl))/(this%ptrMesh_%dxc_(il)*this%ptrMesh_%dxf_(il))
					aT = 0.5d0*(beta%f_(il,jl+1,kl)+beta%f_(il,jl,kl))/(this%ptrMesh_%dyc_(jl+1)*this%ptrMesh_%dyf_(jl))
					aBo = 0.5d0*(beta%f_(il,jl,kl)+beta%f_(il,jl-1,kl))/(this%ptrMesh_%dyc_(jl)*this%ptrMesh_%dyf_(jl))
					aF = 0.5d0*(beta%f_(il,jl,kl+1)+beta%f_(il,jl,kl))/(this%ptrMesh_%dzc_(kl+1)*this%ptrMesh_%dzf_(kl))
					aBa = 0.5d0*(beta%f_(il,jl,kl)+beta%f_(il,jl,kl-1))/(this%ptrMesh_%dzc_(kl)*this%ptrMesh_%dzf_(kl))
					
					!matrix values
					!node (i+1,j,k)
					call setEntryValue(this,-V*aR) 
					!node (i,j+1,k)
					call setEntryValue(this,-V*aT)
					!node (i,j,k+1)
					call setEntryValue(this,-V*aF)
					!node (i,j,k)
					call setEntryValue(this,V*(aR+aL+aT+aBo+aF+aba))
					
				end do

			end do

		end	do
		
		
	end subroutine
!========================================================================================!

!========================================================================================!
	subroutine computeLagrangianMulValues(this) 
		type(poissMat), intent(inout) :: this
		integer :: i, j, k
		integer :: nx, ny, nz 
		
		!assign sizes	
		nx = this%ptrMesh_%nx_
		ny = this%ptrMesh_%ny_
		nz = this%ptrMesh_%nz_
		
		!internal nodes
		do k = 1, nz
			do j = 1, ny
				do i = 1, nx
					
					call setEntryValue(this,1.d0)
					
				end do
			end do
		end	do
		
		
	end subroutine
!========================================================================================!

!========================================================================================!
	subroutine computeBoundaryMatrixValues(this,f,beta) 
		type(poissMat), intent(inout) :: this
		type(scalarField), intent(in) :: f
		type(scalarField), intent(in) :: beta
		integer :: nx, ny, nz
		integer :: q
		
		!mesh size
		nx = this%ptrMesh_%nx_
		ny = this%ptrMesh_%ny_
		nz = this%ptrMesh_%nz_
		
		!off-diagonal boundary values
		do q=1,3
			call computeOffDiagBoundaryValues(this,f%bLeft_,beta,q,nx,ny,nz)
			call computeOffDiagBoundaryValues(this,f%bRight_,beta,q,nx,ny,nz)
			call computeOffDiagBoundaryValues(this,f%bBottom_,beta,q,nx,ny,nz)
			call computeOffDiagBoundaryValues(this,f%bTop_,beta,q,nx,ny,nz)
			call computeOffDiagBoundaryValues(this,f%bBack_,beta,q,nx,ny,nz)
			call computeOffDiagBoundaryValues(this,f%bFront_,beta,q,nx,ny,nz)
		end do
		
		!diagonal boundary values
		do q=1,3
			call computeDiagBoundaryValues(this,f%bLeft_,beta,q,nx,ny,nz)
			call computeDiagBoundaryValues(this,f%bRight_,beta,q,nx,ny,nz)
			call computeDiagBoundaryValues(this,f%bBottom_,beta,q,nx,ny,nz)
			call computeDiagBoundaryValues(this,f%bTop_,beta,q,nx,ny,nz)
			call computeDiagBoundaryValues(this,f%bBack_,beta,q,nx,ny,nz)
			call computeDiagBoundaryValues(this,f%bFront_,beta,q,nx,ny,nz)
		end do
		
		
	end subroutine
!========================================================================================!

!========================================================================================!

	subroutine computeDiagBoundaryValues(this,bf,beta,dir,nx,ny,nz) 
		type(poissMat), intent(inout) :: this
		type(scalarBoundaryField), intent(in) :: bf
		type(scalarField), intent(in) :: beta
		integer, intent(in) :: dir
		integer, intent(in) :: nx, ny, nz 
		integer :: i, j, k, n
		integer :: is, ie, js, je, ks ,ke
		integer :: di, dj, dk, dim, dip, djm, djp, dkm, dkp
		integer :: patchDir, bType
		real(DP) :: ddm, ddp, vm1, vm2, vm, vp1, vp2, vp, bm, bp
		real(DP), dimension(3) :: dSwitch
		real(DP) :: mSwitch, pSwitch
		real(DP) :: V

		bType = bf%bType_
		patchDir = bf%cartDir_
		
		!set indexes bounds
		call setIndexBounds(patchDir,is,ie,js,je,ks,ke,1,nx,1,ny,1,nz)
		!set derivatives parameters
		call setDerivIndexes(dir,dip,dim,djp,djm,dkp,dkm) 
		!set off-set from mesh size
		call setOffSetLoopOverBoundaries(patchDir,dir,di,dj,dk) 
		!set orthogonal direction index to the boundary
		call setOrthBoundaryIndex(n,patchDir,nx,ny,nz)
		!set direction switch
		call setDirSwitch(dir,dSwitch)
		!set patch switch: (0,1) depending on positive or negative normal
		call setPatchSwitch(patchDir,mSwitch,pSwitch)
		
		if ( abs(patchDir) == dir) then	

			SELECT CASE (bType)
				CASE(s_fixedValue)
					do k = ks, ke
						do j= js, je
							do i = is, ie
								V = this%ptrMesh_%V_(i,j,k)
								!contribute node -1
								!contribute -n	
								ddm = this%ptrMesh_%df_(dir,n)*(bf%dhb_ - bf%dhi_)
								bm = 0.5d0*(beta%f_(i,j,k)+beta%f_(i+dim,j+djm,k+dkm))
								vm1 = mSwitch*bm/ddm
								!contribute +n	
								ddm = this%ptrMesh_%dc_(dir,n)*this%ptrMesh_%df_(dir,n)
								bm = 0.5d0*(beta%f_(i,j,k)+beta%f_(i+dim,j+djm,k+dkm))
								vm2 = -pSwitch*bm/ddm
									
								vm = vm1 + vm2
									
								!contribute node +1
								!contribute -n	
								ddp = this%ptrMesh_%dc_(dir,n+1)*this%ptrMesh_%df_(dir,n)
								bp = 0.5d0*(beta%f_(i+dip,j+djp,k+dkp)+beta%f_(i,j,k))
								vp1 = -mSwitch*bp/ddp
								!contribute +n	
								ddp = this%ptrMesh_%df_(dir,n)*(bf%dhb_ - bf%dhi_)
								bp = 0.5d0*(beta%f_(i+dip,j+dkp,k+dkp)+beta%f_(i,j,k))
								vp2 = pSwitch*bp/ddp
									
								vp = vp1 + vp2
									
								call setEntryValue(this,-V*(vm+vp))
							end do
						end do	
					end do
				CASE(s_normalGradient)
					do k = ks, ke
						do j= js, je
							do i = is, ie
								V = this%ptrMesh_%V_(i,j,k)	
								!contribute node +1
								!on patch -n
								ddp = this%ptrMesh_%dc_(dir,n+1)*this%ptrMesh_%df_(dir,n)
								bp = 0.5d0*(beta%f_(i+dip,j+djp,k+dkp)+beta%f_(i,j,k))
								vp = -mSwitch*bp/ddp
								!contribute node -1
								!on patch +n								
								ddm = this%ptrMesh_%dc_(dir,n)*this%ptrMesh_%df_(dir,n)
								bm = 0.5d0*(beta%f_(i,j,k)+beta%f_(i+dim,j+djm,k+dkm))
								vm = -pSwitch*bm/ddm
								call setEntryValue(this,-V*(vp+vm))
							end do
						end do
					end do	
				CASE DEFAULT
					do k = ks, ke
						do j= js, je
							do i = is, ie
								V = this%ptrMesh_%V_(i,j,k)
								!contribute node -1
								ddm = this%ptrMesh_%dc_(dir,n)*this%ptrMesh_%df_(dir,n)
								bm = 0.5d0*(beta%f_(i,j,k)+beta%f_(i+dim,j+djm,k+dkm))
								vm = -bm/ddm
								!contribute node +1
								ddp = this%ptrMesh_%dc_(dir,n+1)*this%ptrMesh_%df_(dir,n)
								bp = 0.5d0*(beta%f_(i+dip,j+djp,k+dkp)+beta%f_(i,j,k))
								vp = -bp/ddp
								call setEntryValue(this,-V*(vm+vp))
							end do
						end do
					end do	
			END SELECT
				
		else

			do k = ks+dk, ke-dk
				do j= js+dj, je-dj
					do i = is+di, ie-di
						V = this%ptrMesh_%V_(i,j,k)
						!contribute node -1
						ddm = dSwitch(1)*this%ptrMesh_%dxc_(i)*this%ptrMesh_%dxf_(i) +	&
							  dSwitch(2)*this%ptrMesh_%dyc_(j)*this%ptrMesh_%dyf_(j) +	&
							  dSwitch(3)*this%ptrMesh_%dzc_(k)*this%ptrMesh_%dzf_(k)
						bm = 0.5d0*(beta%f_(i,j,k)+beta%f_(i+dim,j+djm,k+dkm))
						vm = -bm/ddm
						!contribute node +1
						ddp = dSwitch(1)*this%ptrMesh_%dxc_(i+1)*this%ptrMesh_%dxf_(i) + & 
							  dSwitch(2)*this%ptrMesh_%dyc_(j+1)*this%ptrMesh_%dyf_(j) + &
							  dSwitch(3)*this%ptrMesh_%dzc_(k+1)*this%ptrMesh_%dzf_(k)
						bp = 0.5d0*(beta%f_(i+dip,j+djp,k+dkp)+beta%f_(i,j,k))
						vp = -bp/ddp
						call setEntryValue(this,-V*(vm+vp))
					end do
				end do
			end do
			
		end if
		
		
	end subroutine

!========================================================================================!


!========================================================================================!
	subroutine computeOffDiagBoundaryValues(this,bf,beta,dir,nx,ny,nz) 
		type(poissMat), intent(inout) :: this
		type(scalarBoundaryField), intent(in) :: bf
		type(scalarField), intent(in) :: beta
		integer, intent(in) :: dir
		integer, intent(in) :: nx, ny, nz 	
		integer :: i, j, k, n
		integer :: ie, is, je, js, ke, ks
		integer :: di, dj, dk, dip, djp, dkp, dim, djm, dkm
		integer :: patchDir, bType
		real(DP) :: dd, r
		real(DP), dimension(3) :: dSwitch
		real(DP) :: V

		bType = bf%bType_
		patchDir = bf%cartDir_
		
		!set indexes bounds
		call setIndexBounds(patchDir,is,ie,js,je,ks,ke,1,nx,1,ny,1,nz)
		!set derivatives parameters
		call setDerivIndexes(dir,dip,dim,djp,djm,dkp,dkm) 
		!set off-set from mesh size
		call setOffSetLoopOverBoundaries(patchDir,dir,di,dj,dk) 	
		!set orthogonal direction index to the boundary
		call setOrthBoundaryIndex(n,patchDir,nx,ny,nz)
		!set direction switch
		call setDirSwitch(dir,dSwitch)
		
	
		if ( abs(patchDir) == dir) then
	
			SELECT CASE (patchDir)
   			
   				CASE (-1,-2,-3) 
					!check periodic boundary
						if (bType == s_periodicBC) then 
							do k=ks,ke
								do j=js,je
									do i=is,ie
										V = this%ptrMesh_%V_(i,j,k)
										dd = this%ptrMesh_%dc_(dir,n)*this%ptrMesh_%df_(dir,n)
										r = 0.5d0*(beta%f_(i,j,k)+beta%f_(i+dim,j+djm,k+dkm))/dd
										call setEntryValue(this,-V*r)
									end do
								end do
							end do
						end if 
						!node +
						do k=ks,ke
							do j=js,je
								do i=is,ie
									V = this%ptrMesh_%V_(i,j,k)
									dd = this%ptrMesh_%dc_(dir,n+1)*this%ptrMesh_%df_(dir,n)
									r = 0.5d0*(beta%f_(i+dip,j+djp,k+dkp)+beta%f_(i,j,k))/dd
									call setEntryValue(this,-V*r)
								end do
							end do
						end do
					
				CASE(1,2,3)
					if ( bType == s_parallelBC ) then 
						!node +
						do k=ks,ke
							do j=js,je
								do i=is,ie
									V = this%ptrMesh_%V_(i,j,k)
									dd = this%ptrMesh_%dc_(dir,n+1)*this%ptrMesh_%df_(dir,n)
									r = 0.5d0*(beta%f_(i+dip,j+djp,k+dkp)+beta%f_(i,j,k))/dd
									call setEntryValue(this,-V*r)
								end do
							end do
						end do
					end if
				
			CASE DEFAULT
			END SELECT
				
			else
				
				do k=ks+dk,ke-dk
					do j=js+dj,je-dj
						do i = is+di, ie-di
							V = this%ptrMesh_%V_(i,j,k)
							dd = dSwitch(1)*this%ptrMesh_%dxc_(i+1)*this%ptrMesh_%dxf_(i) +	&
								 dSwitch(2)*this%ptrMesh_%dyc_(j+1)*this%ptrMesh_%dyf_(j) +	&
								 dSwitch(3)*this%ptrMesh_%dzc_(k+1)*this%ptrMesh_%dzf_(k)
							r = 0.5d0*(beta%f_(i+dip,j+djp,k+dkp)+beta%f_(i,j,k))/dd
							call setEntryValue(this,-V*r)
						end do
					end do
				end do
						
			end if
		
		
	end subroutine		

!========================================================================================!

!****************************************************************************************!
!**************************************** RHS *******************************************!
!****************************************************************************************!

!========================================================================================!
	subroutine computeRHS(this,f,beta,s) 
		type(poissMat), intent(inout) :: this
		type(scalarField), intent(in) :: f, beta, s
		integer :: i, j, k, nx, ny, nz
		real(DP) :: V
		
		!reset field to zero
		this%rhs_loc = 0.d0
		
		nx = f%nx_
		ny = f%ny_
		nz = f%nz_
		
		!source term
		do k = 1, nz
			do j = 1, ny
				do i = 1, nx
				
					V = this%ptrMesh_%V_(i,j,k)
					this%rhs_loc(i,j,k) = - V*s%f_(i,j,k)
					
				end do
			end do
		end do
		
		!boundary nodes
		call computeBoundarySourceRHS(this,f%bLeft_,beta,nx,ny,nz)
		call computeBoundarySourceRHS(this,f%bRight_,beta,nx,ny,nz)
		call computeBoundarySourceRHS(this,f%bBottom_,beta,nx,ny,nz)
		call computeBoundarySourceRHS(this,f%bTop_,beta,nx,ny,nz)
		call computeBoundarySourceRHS(this,f%bBack_,beta,nx,ny,nz)
		call computeBoundarySourceRHS(this,f%bFront_,beta,nx,ny,nz)


		!send rhs to Master
		call gatherRHS(this,f)
		
	end subroutine
!========================================================================================!


!========================================================================================!
	subroutine gatherRHS(this,f) 
		type(poissMat), intent(inout) :: this
		type(scalarField), intent(in) :: f
		integer :: i, j, k, nx, ny, nz
		integer :: nxg, nyg, nzg
		integer :: i0g, i1g, j0g, j1g, k0g, k1g
		integer :: il, jl, kl, n
		type(mpiControl), pointer :: ptrMPIC
		integer, dimension(MPI_STATUS_SIZE) :: status
		integer :: ierror 
		integer :: tag = 0
		integer :: q, np
		
		!assign sizes	
		nx = this%ptrMesh_%nx_
		ny = this%ptrMesh_%ny_
		nz = this%ptrMesh_%nz_
		
		nxg = this%ptrMesh_%nxg_
		nyg = this%ptrMesh_%nyg_
		nzg = this%ptrMesh_%nzg_
		
		ptrMPIC => this%ptrMesh_%ptrMPIC_
		np = ptrMPIC%nProcs_
		
		
		!assign master mumps rhs first
		if (IS_MASTER) then
		
			!global indexes
			i0g = this%ptrMesh_%i0g_
			j0g = this%ptrMesh_%j0g_
			k0g = this%ptrMesh_%k0g_
			i1g = this%ptrMesh_%i1g_
			j1g = this%ptrMesh_%j1g_
			k1g = this%ptrMesh_%k1g_
		
				do k=k0g,k1g
					do j=j0g,j1g
						do i=i0g,i1g
							
							call diagMatIndex(n,i,j,k,nxg,nyg)
							
							il = i-i0g+1
							jl = j-j0g+1
							kl = k-k0g+1
							
							this%mumpsPar%rhs(n) = this%rhs_loc(il,jl,kl)
							
						end do
					end do
				end do
			
		end if
		
		!gather mumps rhs 
		if (IS_MASTER) then
			
			do q=1,np-1
				 
				call MPI_Recv(this%rhs_loc,size(this%rhs_loc),MPI_DOUBLE_PRECISION,q,tag,ptrMPIC%cartComm_,status,ierror)
				 
				!set global indexes proc q
				call globalIndexesFromAll(this%ptrMesh_,i0g,i1g,j0g,j1g,k0g,k1g,q,f%tp_) 
				
				do k=k0g,k1g
					do j=j0g,j1g
						do i=i0g,i1g
							
							call diagMatIndex(n,i,j,k,nxg,nyg)
							
							il = i-i0g+1
							jl = j-j0g+1
							kl = k-k0g+1
							
							this%mumpsPar%rhs(n) = this%rhs_loc(il,jl,kl)
							
						end do
					end do
				end do
				
			end do
		else
		
			call MPI_Send(this%rhs_loc,size(this%rhs_loc),MPI_DOUBLE_PRECISION,0,tag,ptrMPIC%cartComm_,ierror)
			
		end if
		
	end subroutine
!========================================================================================!

!========================================================================================!
	subroutine computeBoundarySourceRHS(this,bf,beta,nx,ny,nz) 
		type(poissMat), intent(inout) :: this
		type(scalarBoundaryField), intent(in) :: bf
		type(scalarField), intent(in) :: beta
		integer, intent(in) :: nx, ny, nz
		integer :: i, j, k ,n
		integer :: ie, is, je, js, ke, ks
		integer :: dip, djp, dkp, dim, djm, dkm
		integer :: patchDir, bType, dir
		real(DP) :: dd, v1, v2, r, b
		real(DP) :: mSwitch, pSwitch
		real(DP) :: V
		
		bType = bf%bType_
		patchDir = bf%cartDir_	
		dir  = abs(patchDir)
		
		!set indexes bounds
		call setIndexBounds(patchDir,is,ie,js,je,ks,ke,1,nx,1,ny,1,nz)
		!set derivatives parameters
		call setDerivIndexes(dir,dip,dim,djp,djm,dkp,dkm) 	
		!set orthogonal direction index to the boundary
		call setOrthBoundaryIndex(n,patchDir,nx,ny,nz)
		!set patch switch: (0,1) depending on positive or negative normal
		call setPatchSwitch(patchDir,mSwitch,pSwitch)
		
				
		SELECT CASE(bType)
				
			CASE(s_fixedValue)
			
				do k=ks,ke
					do j=js,je
						do i=is,ie
							V = this%ptrMesh_%V_(i,j,k)
							!on patch -n
							dd = this%ptrMesh_%df_(dir,n)*(bf%dhb_ - bf%dhi_)
							b = 0.5d0*(beta%f_(i,j,k)+beta%f_(i+dim,j+djm,k+dkm))
							v1 = mSwitch*bf%bf_*b/dd
							!on patch +n
							dd = this%ptrMesh_%df_(dir,n+1)*(bf%dhb_ - bf%dhi_)
							b = 0.5d0*(beta%f_(i+dip,j+djp,k+dkp)+beta%f_(i,j,k))
							v2 = pSwitch*bf%bf_*b/dd
									
							r = V*(v1 + v2)
									
							this%rhs_loc(i,j,k) = this%rhs_loc(i,j,k) - r				
		
						end do
					end do
				end do
						
			CASE(s_normalGradient)
				do k=ks,ke
					do j=js,je
						do i=is,ie
							V = this%ptrMesh_%V_(i,j,k)
							!on patch -n
							b = 0.5d0*(beta%f_(i,j,k)+beta%f_(i+dim,j+djm,k+dkm))
							v1 = -mSwitch*bf%bf_*b
							!on patch +n
							b = 0.5d0*(beta%f_(i+dip,j+djp,k+dkp)+beta%f_(i,j,k))
							v2 = -pSwitch*bf%bf_*b								
									
							r = V*(v1 + v2)
									
							this%rhs_loc(i,j,k) = this%rhs_loc(i,j,k) - r
						end do
					end do
				end do
			CASE DEFAULT
					
		END SELECT
				
		
	end subroutine
!========================================================================================!

!****************************************************************************************!
!************************************* helpers ******************************************!
!****************************************************************************************!

!========================================================================================!
	subroutine setEntryIndex(this,r,c) 
		type(poissMat), intent(inout) :: this
		integer, intent(in) :: r, c

		this%entryCount_ = this%entryCount_+1
		this%mumpsPar%irn_loc(this%entryCount_)= r
		this%mumpsPar%jcn_loc(this%entryCount_)= c

	end subroutine
!========================================================================================!

!========================================================================================!
	subroutine setEntryValue(this,v) 
		type(poissMat), intent(inout) :: this
		real(DP), intent(in) :: v

		this%entryCount_ = this%entryCount_+1
		this%mumpsPar%A_loc(this%entryCount_)= v

	end subroutine
!========================================================================================!

!========================================================================================!
	subroutine diagMatIndex(n,i,j,k,nx,ny) 
		integer, intent(in) :: i,j,k,nx,ny
		integer, intent(out) :: n
	
		n = i + (j-1)*nx + (k-1)*nx*ny
	
	end subroutine
!========================================================================================!

!========================================================================================!
	subroutine setOrthBoundaryIndex(n,patchDir,nx,ny,nz) 
		integer, intent(in) ::  patchDir, nx, ny, nz
		integer, intent(out) :: n
		
		SELECT CASE (patchDir)
			CASE(-1)
				n=1
			CASE(1)
				n=nx
			CASE(-2)
				n=1
			CASE(2)
				n=ny
			CASE(-3)
				n=1
			CASE(3)
				n=nz
			CASE DEFAULT
		END SELECT
		
	end subroutine
!========================================================================================!


!========================================================================================!
	subroutine setIndexBounds(patchDir,is,ie,js,je,ks,ke,i0,i1,j0,j1,k0,k1) 
		integer, intent(in) ::  patchDir, i0,i1,j0,j1,k0,k1
		integer, intent(out) :: is,ie,js,je,ks,ke
		
		is = i0
		ie = i1
		js = j0
		je = j1
		ks = k0
		ke = k1
		
		SELECT CASE (patchDir)
			CASE(-1)
				is = i0
				ie = i0
			CASE(1)
				is = i1
				ie = i1 
			CASE(-2)
				js = j0
				je = j0
			CASE(2)
				js = j1
				je = j1
			CASE(-3)
				ks = k0
				ke = k0
			CASE(3)
				ks = k1
				ke = k1
   			CASE DEFAULT
		END SELECT	
		

		
	end subroutine
!========================================================================================!


!========================================================================================!
	subroutine setDerivIndexes(dir,dip,dim,djp,djm,dkp,dkm) 
		integer, intent(in) :: dir
		integer, intent(out) :: dip, dim, djp, djm, dkp, dkm

		SELECT CASE (dir)
   		CASE (1) 
   			dip = 1
			dim = -1
			djp = 0
			djm = 0
			dkp = 0
			dkm = 0
		CASE (2)
   			dip = 0
			dim = 0
			djp = 1
			djm = -1
			dkp = 0
			dkm = 0
   		CASE (3) 
   			dip = 0
			dim = 0
			djp = 0
			djm = 0
			dkp = 1
			dkm = -1
   		CASE DEFAULT
		END SELECT
		
	end subroutine
!========================================================================================!


!========================================================================================!
	subroutine setDirSwitch(dir,dSwitch) 
		integer, intent(in) :: dir
		real(DP), dimension(3), intent(out) :: dSwitch

		SELECT CASE (dir)
   		CASE (1) 
			dSwitch(1) = 1.d0
			dSwitch(2) = 0.d0
			dSwitch(3) = 0.d0
		CASE (2)
			dSwitch(1) = 0.d0
			dSwitch(2) = 1.d0
			dSwitch(3) = 0.d0
   		CASE (3) 
			dSwitch(1) = 0.d0
			dSwitch(2) = 0.d0
			dSwitch(3) = 1.d0
   		CASE DEFAULT
		END SELECT
		
	end subroutine
!========================================================================================!

!========================================================================================!
	subroutine setPatchSwitch(patchDir,mSwitch,pSwitch) 
		integer, intent(in) :: patchDir
		real(DP), intent(out) :: mSwitch, pSwitch

		SELECT CASE (patchDir)
   		CASE (-1,-2,-3) 
			mSwitch = 1.d0
			pSwitch = 0.d0
		CASE (1,2,3)
			mSwitch = 0.d0
			pSwitch = 1.d0
   		CASE DEFAULT
		END SELECT
		
	end subroutine
!========================================================================================!

!========================================================================================!
	subroutine setOffSetLoopOverBoundaries(patchDir,dir,di,dj,dk) 
		integer, intent(in) ::  patchDir, dir
		integer, intent(out) :: di, dj, dk		
		
		
		SELECT CASE (dir)
   		CASE (1) 
			SELECT CASE (patchDir)
			CASE(-1,1)
				di = 0
				dj = 0
				dk = 0
			CASE(-2,2)
				di = 1
				dj = 0
				dk = 0
			CASE(-3,3)	
				di = 1
				dj = 1
				dk = 0
   			CASE DEFAULT
			END SELECT	
		CASE (2)
			SELECT CASE (patchDir)
			CASE(-1,1)
				di = 0
				dj = 1
				dk = 0
			CASE(-2,2)
				di = 0
				dj = 0
				dk = 0
			CASE(-3,3)
				di = 1
				dj = 1
				dk = 0	
   			CASE DEFAULT
			END SELECT	
   		CASE (3) 
			SELECT CASE (patchDir)
			CASE(-1,1)
				di = 0
				dj = 0
				dk = 1
			CASE(-2,2)
				di = 1
				dj = 0
				dk = 1
			CASE(-3,3)	
				di = 0
				dj = 0
				dk = 0
   			CASE DEFAULT
			END SELECT
   		CASE DEFAULT
		END SELECT
		

		
	end subroutine
!========================================================================================!

!========================================================================================!
	subroutine checkSingularSystem(this,f) 
		type(poissMat), intent(inout) :: this
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


end module poissMatMod

