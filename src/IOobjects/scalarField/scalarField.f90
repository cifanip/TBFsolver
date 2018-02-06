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

module scalarFieldMod

	use interpolationMod
	use ompRoutinesMod
	use gridMod
	
	implicit none
	
	!boundary types
	integer, parameter :: s_fixedValue = 1
	integer, parameter :: s_normalGradient = 2
	integer, parameter :: s_periodicBC = 3
	integer, parameter :: s_parallelBC = 4
	integer, parameter :: s_calculatedBC = 5
	integer, parameter :: s_contactAngleBC = 6
	
	character(len=0), private, parameter :: s_fileDir = ''
	
	! scalarBoundaryField
	INCLUDE 'scalarBoundaryField_H.f90'
	
	! scalarField
	type, public :: scalarField
	
		character(len=:), allocatable :: fileName_
		character(len=:), allocatable :: filePath_
	
		!class pointer for multi-grids or previous time-levels
		type(scalarField), pointer :: ptrf_ => NULL()
	
		!keep a pointer to grid
		type(grid), pointer :: ptrMesh_ => NULL() 
		
		!copy of mesh size
		integer :: nx_
		integer :: ny_
		integer :: nz_
		
		!start and end-indexes (minus the halo)
		integer :: is_, js_, ks_, ie_, je_, ke_
		
		!field type
		character(len=2) :: tp_
		
		!halo dimension
		integer :: hd_
		
		!field + halo	
		real(DP), allocatable, dimension(:,:,:) :: f_
		
		!used only in MG for prolongation
		real(DP), allocatable, dimension(:,:,:) :: prol_
		
		!boundaries
		type(scalarBoundaryField) :: bRight_, bLeft_
		type(scalarBoundaryField) :: bBottom_, bTop_
		type(scalarBoundaryField) :: bBack_, bFront_
		
		!halo patches
		integer, private :: xPatch_, yPatch_, zPatch_
		integer :: xPatchEq_, yPatchEq_, zPatchEq_
		
		contains
		
		final :: delete_scalarField
		
	end type
	
	
	private :: readField
	private :: readInternalField
	private :: readDefaultField
	private :: updateHalo
	private :: buildHaloExchangeTypes
	private :: decomposeInternalScalarField 
	private :: decomposeBoundaryScalarField 
!DIR$ IF DEFINED (MG_MODE)	
	private :: updateEdges
	private :: updateCorners
	private :: coarsenField
!DIR$ ENDIF		
	private :: deallocatePtrf
	private :: checkPeriodicWrap
	private :: setIndexesBounds
	private :: writeField
	
	public :: scalarFieldCTOR
	public :: delete_scalarField
	public :: allocatePtrf
	public :: initToZero
	public :: decomposeField
	public :: updateBoundaries
!DIR$ IF DEFINED (MG_MODE)	
	public :: coarsenFields
	public :: resetBCerrorField
	public :: prolongateField
	public :: reconstructAndWriteField
!DIR$ ENDIF		
	public :: copyBoundary

!DIR$ IF DEFINED (MG_MODE)	
	interface restrictField
		module PROCEDURE restrict_scalarField
		module PROCEDURE restrict_extField_scalarField
	end interface
!DIR$ ENDIF
	
contains

!******************************** boundaryScalarField ***********************************!

INCLUDE 'scalarBoundaryField_S.f90'

!************************************** scalarField *************************************!
!========================================================================================!
    subroutine delete_scalarField(this)
        type(scalarField), intent(inout) :: this
        integer :: ierror
        
        call deallocatePtrf(this)
        
        if (associated(this%ptrMesh_)) then
        	call MPI_TYPE_FREE(this%xPatch_, ierror)
        	call MPI_TYPE_FREE(this%yPatch_, ierror)
        	call MPI_TYPE_FREE(this%zPatch_, ierror)
        	call MPI_TYPE_FREE(this%xPatchEq_, ierror)
        	call MPI_TYPE_FREE(this%yPatchEq_, ierror)
        	call MPI_TYPE_FREE(this%zPatchEq_, ierror)
        end if
        
    end subroutine
!========================================================================================!

!========================================================================================!
	subroutine scalarFieldCTOR(this,fileName,mesh,fType,hd,initOpt,nFolder)
		type(scalarField), intent(out) :: this
		type(grid), intent(in), target :: mesh
		character(len=*), intent(in) :: fileName
		character(len=2), intent(in) :: fType
		integer, intent(in) :: hd
		integer, intent(in) :: initOpt
		integer, intent(in), optional :: nFolder
		integer :: opt_nFolder
		type(mpiControl), pointer :: ptrMPIC
		integer :: nx, ny, nz
		
		!optional argument
		if (present(nFolder)) then
			opt_nFolder = nFolder
		else
			opt_nFolder = 0
		end if
		
		this%fileName_ = fileName
		this%filePath_ = s_fileDir//'/'//fileName
		
		this%ptrMesh_ => mesh
		ptrMPIC => mesh%ptrMPIC_
		
		nx = mesh%nx_
		ny = mesh%ny_
		nz = mesh%nz_

		this%nx_ = nx
		this%ny_ = ny
		this%nz_ = nz
		
		!field type 
		this%tp_ = fType
		
		!halo dimension
		this%hd_ = hd
		
		if (hd > mesh%hd_) then
			call mpiABORT('Attempt to assign halo larger than the mesh halo ')
		end if

		!indexes 
		call setIndexesBounds(this,this%tp_)

		call allocateArray(this%f_,this%is_-hd,this%ie_+hd,	&
						   this%js_-hd,this%je_+hd,			&
						   this%ks_-hd,this%ke_+hd)

		!init options
		SELECT CASE (initOpt)
   			CASE (-1)
				!init to 0 whole field (halo included)
				call initToZero(this)
   			CASE (0)
				!init to 0 whole field (halo included)
				call initToZero(this)
				!default zero-grad BC
				call readDefaultField(this)
   			CASE (1)
				!init to 0 whole field (halo included)
				call initToZero(this)
				!read internal field + boundary patches
				call readField(this,nFolder)
   			CASE DEFAULT
   				call mpiABORT('Invalid init field option ')
		END SELECT
		
	end subroutine
!========================================================================================!

!========================================================================================!
    subroutine readField(this,nFolder)
        type(scalarField), intent(inout) :: this
        integer, intent(in) :: nFolder
        character(len=10) :: dirName
        
        write(dirName,s_intFormat) nFolder
        
        open(UNIT=s_IOunitNumber,FILE=adjustl(trim(dirName)//this%filePath_), &
        	 form='UNFORMATTED',ACCESS='STREAM',STATUS='OLD',ACTION='READ')
        
        !read internal field
		call readInternalField(this)

		!read boundaries
		call scalarBoundaryFieldCTOR(this%bLeft_,this,1,.TRUE.)
		call scalarBoundaryFieldCTOR(this%bRight_,this,2,.TRUE.)
		call scalarBoundaryFieldCTOR(this%bBottom_,this,3,.TRUE.)
		call scalarBoundaryFieldCTOR(this%bTop_,this,4,.TRUE.)
		call scalarBoundaryFieldCTOR(this%bBack_,this,5,.TRUE.)
		call scalarBoundaryFieldCTOR(this%bFront_,this,6,.TRUE.)
		
		!consistency check periodic patches - wrap around option
		call checkPeriodicWrap(this)
		
		close(s_IOunitNumber)
        
    end subroutine
!========================================================================================!

!========================================================================================!
    subroutine readDefaultField(this)
        type(scalarField), intent(inout) :: this
        
		!init default boundaries
		call scalarBoundaryFieldCTOR(this%bLeft_,this,1,.FALSE.)
		call scalarBoundaryFieldCTOR(this%bRight_,this,2,.FALSE.)
		call scalarBoundaryFieldCTOR(this%bBottom_,this,3,.FALSE.)
		call scalarBoundaryFieldCTOR(this%bTop_,this,4,.FALSE.)
		call scalarBoundaryFieldCTOR(this%bBack_,this,5,.FALSE.)
		call scalarBoundaryFieldCTOR(this%bFront_,this,6,.FALSE.) 
      
    end subroutine
!========================================================================================!

!========================================================================================!
    subroutine checkPeriodicWrap(this)
        type(scalarField), intent(inout) :: this
        type(mpiControl), pointer :: ptrMPIC
        logical, dimension(3) :: wrap
        
		ptrMPIC => this%ptrMesh_%ptrMPIC_
		wrap = ptrMPIC%wrapAround_
		
		if (wrap(1) .eqv. .TRUE.) then
			if ((this%bLeft_%bType_ /= s_periodicBC) &
			   .OR. (this%bRight_%bType_ /= s_periodicBC)) then
				call mpiABORT('Boundary Patch in x direction not periodic but MPI wrap = true ')  
			end if
		end if
		
		if (wrap(2) .eqv. .TRUE.) then
			if ((this%bTop_%bType_ /= s_periodicBC) &
			   .OR. (this%bBottom_%bType_ /= s_periodicBC)) then
				call mpiABORT('Boundary Patch in y direction not periodic but MPI wrap = true ')  
			end if
		end if
		
		if (wrap(3) .eqv. .TRUE.) then
			if ((this%bBack_%bType_ /= s_periodicBC) &
			   .OR. (this%bFront_%bType_ /= s_periodicBC)) then
				call mpiABORT('Boundary Patch in z direction not periodic but MPI wrap = true ')  
			end if
		end if
        
    end subroutine
!========================================================================================!

!========================================================================================!
    subroutine readInternalField(this)
        type(scalarField), intent(inout) :: this
        real(DP), allocatable, dimension(:,:,:) :: tmp
        integer :: n, size
        
        call allocateArray(tmp,this%is_,this%ie_,	&
        					   this%js_,this%je_,	&
        					   this%ks_,this%ke_)
        
        size = (this%ie_-this%is_+1)*(this%je_-this%js_+1)*(this%ke_-this%ks_+1)

		read(s_IOunitNumber) n
		if (n /= size) then
			call mpiABORT('Dimension of field '//this%fileName_//' not equal to mesh size')
		end if
		
		read(s_IOunitNumber) tmp(:,:,:)
		
		call assign_omp(this%f_,tmp,this%is_,this%ie_,this%js_,this%je_,this%ks_,this%ke_)			
        
    end subroutine
!========================================================================================!

!========================================================================================!
    subroutine writeField(this,nFolder)
        type(scalarField), intent(inout) :: this
        integer, intent(in) :: nFolder
        real(DP), allocatable, dimension(:,:,:) :: tmp
        integer :: size
        character(len=20) :: dirName
        integer :: ios
        
        
        if (IS_MASTER) then
        
        	call allocateArray(tmp,this%is_,this%ie_,	&
        					       this%js_,this%je_,	&
        					       this%ks_,this%ke_)
        	
        	call assign_omp(tmp,this%f_,this%is_,this%ie_,this%js_,this%je_,this%ks_,this%ke_)				   
        
        	write(dirName,s_intFormat) nFolder
        
       		size = (this%ie_-this%is_+1)*(this%je_-this%js_+1)*(this%ke_-this%ks_+1)
        
        	open(UNIT=s_IOunitNumber,FILE=adjustl(trim(dirName)//this%filePath_), &
        	 		form='UNFORMATTED',ACCESS='STREAM',STATUS='REPLACE',ACTION='WRITE')
        
        		!write internal field
        		write(s_IOunitNumber) size
        		
        		write(s_IOunitNumber) tmp
			
				!write boundary field
				call writeBoundaryField(this%bLeft_)
				call writeBoundaryField(this%bRight_)
				call writeBoundaryField(this%bBottom_)
				call writeBoundaryField(this%bTop_)
				call writeBoundaryField(this%bBack_)
				call writeBoundaryField(this%bFront_)
			
			close(unit=s_IOunitNumber,iostat=ios,STATUS='KEEP')
		
		end if
        
    end subroutine
!========================================================================================!


!========================================================================================!
    subroutine initToZero(this)
        type(scalarField), intent(inout) :: this
        
        call set2zero_omp(this%f_)
        
    end subroutine
!========================================================================================!

!========================================================================================!
    subroutine decomposeField(this,f,dealloc_gfield)
        type(scalarField), intent(inout) :: this
    	type(scalarField), intent(inout) :: f
    	logical, intent(in), optional :: dealloc_gfield
    	logical :: opt_dealloc_gfield
    	integer :: ierror
    	
		if (present(dealloc_gfield)) then
			opt_dealloc_gfield = dealloc_gfield
		else
			opt_dealloc_gfield = .TRUE.
		end if
    	
    	!check types
		if (IS_MASTER) then
    		if (this%tp_ /= f%tp_) then
    			call mpiABORT('Attempt to decompose fields of different types ')
    		end if
    	end if
    	
    	!check halo dim
		if (IS_MASTER) then
    		if (this%hd_ /= f%hd_) then
    			call mpiABORT('Attempt to decompose fields with different halo dim ')
    		end if
    	end if

		!decompose internal field
		call decomposeInternalScalarField(this,f)
		
		!build halo exchange types
		call buildHaloExchangeTypes(f)

		!decompose boundary field
		call decomposeBoundaryScalarField(this,f)
		
		call MPI_BARRIER(f%ptrMesh_%ptrMPIC_%cartComm_,ierror)
!DIR$ IF DEFINED (MEM_SAVE)
		if ((IS_MASTER).AND.(opt_dealloc_gfield)) then
			call deallocateArray(this%f_)
		end if
!DIR$ ENDIF
        
    end subroutine
!========================================================================================!


!========================================================================================!
   subroutine decomposeInternalScalarField(this,lf)
        type(scalarField), intent(in) :: this
        type(scalarField), intent(inout) :: lf    	
        integer :: ierror
        integer, dimension(MPI_STATUS_SIZE) :: status
        integer :: i0, i1, j0, j1, k0, k1, n
        integer :: i, j, k
        integer :: il, jl, kl
        integer :: tag = 0
        type(mpiControl), pointer :: ptrMPIC

        
		ptrMPIC => lf%ptrMesh_%ptrMPIC_
		
		!decompose domain and send to slaves
		if (IS_MASTER) then
			do n=1,ptrMPIC%nProcs_-1
			
				!set global indexes proc n
				call globalIndexesFromAll(lf%ptrMesh_,i0,i1,j0,j1,k0,k1,n,lf%tp_) 
		
				il = 0
				jl = 0
				kl = 0		
				do k=k0,k1
					do j=j0,j1
						do i=i0,i1
							lf%f_(lf%is_ + il, lf%js_ + jl, lf%ks_ + kl) = this%f_(i,j,k)
							il = il +1
						end do
						il = 0
						jl = jl +1
					end do
					il = 0
					jl = 0
					kl = kl + 1
				end do
				
				!note: only the internal field is updated
				call MPI_Send(lf%f_,size(lf%f_),MPI_DOUBLE_PRECISION,n,tag,ptrMPIC%cartComm_,ierror)
			end do
		
		else
			call MPI_Recv(lf%f_,size(lf%f_),MPI_DOUBLE_PRECISION,0,tag,ptrMPIC%cartComm_,status,ierror)
		end if
        
		!decompose domain and assign to master
		if (IS_MASTER) then

			!set global indexes proc 0
			call globalIndexesFromAll(lf%ptrMesh_,i0,i1,j0,j1,k0,k1,0,lf%tp_) 

			il = 0
			jl = 0
			kl = 0		
			do k=k0,k1
				do j=j0,j1
					do i=i0,i1
						lf%f_(lf%is_ + il, lf%js_ + jl, lf%ks_ + kl) = this%f_(i,j,k)
						il = il + 1
					end do
					il = 0
					jl = jl + 1
				end do
				il = 0
				jl = 0
				kl = kl + 1
			end do
			
		end if
        
    end subroutine
!========================================================================================!


!========================================================================================!
   subroutine decomposeBoundaryScalarField(this,lf)
        type(scalarField), intent(in) :: this
        type(scalarField), intent(inout) :: lf    	
        type(mpiControl), pointer :: ptrMPIC
        integer, dimension(3) :: npax
        integer :: npx, npy, npz
        
        ptrMPIC => lf%ptrMesh_%ptrMPIC_

		npax = ptrMPIC%nProcsAxis_
		npx = npax(1)
		npy = npax(2)
		npz = npax(3)		
        
        !decompose boundary patches
		call decomposeBoundaryPatch(this%bLeft_,lf%bLeft_,ptrMPIC,-1,0)
		call decomposeBoundaryPatch(this%bRight_,lf%bRight_,ptrMPIC,1,npx-1)
		call decomposeBoundaryPatch(this%bBottom_,lf%bBottom_,ptrMPIC,-2,0)
		call decomposeBoundaryPatch(this%bTop_,lf%bTop_,ptrMPIC,2,npy-1)
		call decomposeBoundaryPatch(this%bBack_,lf%bBack_,ptrMPIC,-3,0)
		call decomposeBoundaryPatch(this%bFront_,lf%bFront_,ptrMPIC,3,npz-1)

		!update boundaries
		call updateBoundaries(lf)

    end subroutine
!========================================================================================!


!========================================================================================!
   subroutine buildHaloExchangeTypes(this)
        type(scalarField), intent(inout) :: this
        integer :: ierror
        integer :: sx, sy, sz, hd
        
        sx = (this%ie_-this%is_+1)
        sy = (this%je_-this%js_+1)
        sz = (this%ke_-this%ks_+1)
        
        hd = this%hd_
        
        
        !**************************** HALO EXCHANGE *************************!
		!*********************   left - right boundary  *********************!
		!memory layout patch along \hat x
        call MPI_TYPE_vector( (sy+2*hd)*(sz+2*hd), hd , (2*hd+sx) , &
        					  MPI_DOUBLE_PRECISION , this%xPatch_,ierror)
		call MPI_Type_commit(this%xPatch_,ierror)
		
		!*********************  bottom - top boundary  *********************!	
		!memory layout patch along \hat y
        call MPI_TYPE_vector( sz+2*hd, hd*(sx+2*hd) , (2*hd+sx)*(2*hd+sy) , &
        					  MPI_DOUBLE_PRECISION , this%yPatch_,ierror)
		call MPI_Type_commit(this%yPatch_,ierror)
		
		!*********************  back - front boundary  *********************!
        !memory layout patch along \hat z
        call MPI_TYPE_contiguous( hd*(sy+2*hd)*(sx +2*hd) ,MPI_DOUBLE_PRECISION , this%zPatch_,ierror)
		call MPI_Type_commit(this%zPatch_,ierror)
		
		
		!********************* EQUATE PERIODIC BOUNDARY TYPE ****************!
		!*********************   left - right boundary  *********************!
		!memory layout patch along \hat x
        call MPI_TYPE_vector( (sy+2*hd)*(sz+2*hd), 1 , (2*hd+sx) , &
        					  MPI_DOUBLE_PRECISION , this%xPatchEq_,ierror)
		call MPI_Type_commit(this%xPatchEq_,ierror)
		
		!*********************  bottom - top boundary  *********************!	
		!memory layout patch along \hat y
        call MPI_TYPE_vector( sz+2*hd, (sx+2*hd) , (2*hd+sx)*(2*hd+sy) , &
        					  MPI_DOUBLE_PRECISION , this%yPatchEq_,ierror)
		call MPI_Type_commit(this%yPatchEq_,ierror)
		
		!*********************  back - front boundary  *********************!
        !memory layout patch along \hat z
        call MPI_TYPE_contiguous( (sy+2*hd)*(sx +2*hd) ,MPI_DOUBLE_PRECISION , this%zPatchEq_,ierror)
		call MPI_Type_commit(this%zPatchEq_,ierror)


    end subroutine
!========================================================================================!

!========================================================================================!
   subroutine updateHalo(this)
        type(scalarField), intent(inout) :: this
        type(mpiControl), pointer :: ptrMPIC
        integer, dimension(MPI_STATUS_SIZE,4) :: status
        integer, dimension(4) :: requests
        integer :: ierror
        integer :: is, js, ks, ie, je, ke, hd, tag1, tag2
        
        
        ptrMPIC => this%ptrMesh_%ptrMPIC_
        
		is = this%is_
		ie = this%ie_
		js = this%js_
		je = this%je_
		ks = this%ks_
		ke = this%ke_
		
		hd = this%hd_
		
		tag1=0
		tag2=1

        ! **************************  communication -  halo  ************************** !

        
        SELECT CASE (this%tp_)
        	
        	CASE('cl')
        	
        	!******************************* right-left *********************************!
  			!recv from right
  			call MPI_IRECV(this%f_(ie+1,js-hd,ks-hd),1,this%xPatch_, ptrMPIC%rightNe_, &
  					       tag1, ptrMPIC%cartComm_, requests(1), ierror)
			!recv from left     		
  			call MPI_IRECV(this%f_(is-hd,js-hd,ks-hd),1,this%xPatch_, ptrMPIC%leftNe_, &
  					       tag2, ptrMPIC%cartComm_, requests(2), ierror)         
        	!send to left
        	call MPI_ISSEND(this%f_(is,js-hd,ks-hd), 1, this%xPatch_, ptrMPIC%leftNe_, &
        			   	   tag1, ptrMPIC%cartComm_, requests(4), ierror) 		   	         	
        	!send to right
        	call MPI_ISSEND(this%f_(ie-(hd-1),js-hd,ks-hd), 1, this%xPatch_, ptrMPIC%rightNe_, &
        			   	   tag2, ptrMPIC%cartComm_, requests(3), ierror) 
   	   
        	call MPI_WAITALL(4, requests, status, ierror)     

 	               	
        	!******************************* top-bottom *********************************!
  			!recv from top
  			call MPI_IRECV(this%f_(is-hd,je+1,ks-hd),1,this%yPatch_, ptrMPIC%topNe_, &
  					       tag1, ptrMPIC%cartComm_, requests(1), ierror)	       		   	   
			!recv from bottom   		
  			call MPI_IRECV(this%f_(is-hd,js-hd,ks-hd),1,this%yPatch_, ptrMPIC%bottomNe_, &
  					       tag2, ptrMPIC%cartComm_, requests(2), ierror)        	
        	!send to bottom
        	call MPI_ISSEND(this%f_(is-hd,js,ks-hd), 1, this%yPatch_, ptrMPIC%bottomNe_, &
        			   	   tag1, ptrMPIC%cartComm_, requests(4), ierror) 
        	!send to top
        	call MPI_ISSEND(this%f_(is-hd,je-(hd-1),ks-hd), 1, this%yPatch_, ptrMPIC%topNe_, &
        			   	   tag2, ptrMPIC%cartComm_, requests(3), ierror) 
        			   	   
        	call MPI_WAITALL(4, requests, status, ierror) 

        	        	
        	!******************************* front-back *********************************!
  			!recv from front
  			call MPI_IRECV(this%f_(is-hd,js-hd,ke+1),1,this%zPatch_, ptrMPIC%frontNe_, &
  					       tag1, ptrMPIC%cartComm_, requests(1), ierror)
			!recv from back   		
  			call MPI_IRECV(this%f_(is-hd,js-hd,ks-hd),1,this%zPatch_, ptrMPIC%backNe_, &
  					       tag2, ptrMPIC%cartComm_, requests(2), ierror) 
        	!send to back
        	call MPI_ISSEND(this%f_(is-hd,js-hd,ks), 1, this%zPatch_, ptrMPIC%backNe_, &
        			   	   tag1, ptrMPIC%cartComm_, requests(4), ierror)       	
        	!send to front
        	call MPI_ISSEND(this%f_(is-hd,js-hd,ke-(hd-1)), 1, this%zPatch_, ptrMPIC%frontNe_, &
        			   	   tag2, ptrMPIC%cartComm_, requests(3), ierror) 
        			   	   
        	call MPI_WAITALL(4, requests, status, ierror)      	

        				  		
        	CASE('sx')
        	
        	!******************************* right-left *********************************!
  			!recv from right
  			call MPI_IRECV(this%f_(ie+1,js-hd,ks-hd),1,this%xPatch_, ptrMPIC%rightNe_, &
  					       tag1, ptrMPIC%cartComm_, requests(1), ierror)
			!recv from left     		
  			call MPI_IRECV(this%f_(is-hd,js-hd,ks-hd),1,this%xPatch_, ptrMPIC%leftNe_, &
  					       tag2, ptrMPIC%cartComm_, requests(2), ierror)   
        	!send to left
        	call MPI_ISSEND(this%f_(is+1,js-hd,ks-hd), 1, this%xPatch_, ptrMPIC%leftNe_, &
        			   	   tag1, ptrMPIC%cartComm_, requests(4), ierror)      	
        	!send to right
        	call MPI_ISSEND(this%f_(ie-hd,js-hd,ks-hd), 1, this%xPatch_, ptrMPIC%rightNe_, &
        			   	   tag2, ptrMPIC%cartComm_, requests(3), ierror) 
        			   	   
        	call MPI_WAITALL(4, requests, status, ierror)   
        	       	      	
        	!******************************* top-bottom *********************************!
  			!recv from top
  			call MPI_IRECV(this%f_(is-hd,je+1,ks-hd),1,this%yPatch_, ptrMPIC%topNe_, &
  					       tag1, ptrMPIC%cartComm_, requests(1), ierror)
			!recv from bottom   		
  			call MPI_IRECV(this%f_(is-hd,js-hd,ks-hd),1,this%yPatch_, ptrMPIC%bottomNe_, &
  					       tag2, ptrMPIC%cartComm_, requests(2), ierror)   
        	!send to bottom
        	call MPI_ISSEND(this%f_(is-hd,js,ks-hd), 1, this%yPatch_, ptrMPIC%bottomNe_, &
        			   	   tag1, ptrMPIC%cartComm_, requests(4), ierror)      	
        	!send to top
        	call MPI_ISSEND(this%f_(is-hd,je-(hd-1),ks-hd), 1, this%yPatch_, ptrMPIC%topNe_, &
        			   	   tag2, ptrMPIC%cartComm_, requests(3), ierror) 
        			   	   
        	call MPI_WAITALL(4, requests, status, ierror) 
        	       	      	
        	!******************************* front-back *********************************!
  			!recv from front
  			call MPI_IRECV(this%f_(is-hd,js-hd,ke+1),1,this%zPatch_, ptrMPIC%frontNe_, &
  					       tag1, ptrMPIC%cartComm_, requests(1), ierror)
			!recv from back   		
  			call MPI_IRECV(this%f_(is-hd,js-hd,ks-hd),1,this%zPatch_, ptrMPIC%backNe_, &
  					       tag2, ptrMPIC%cartComm_, requests(2), ierror)  
        	!send to back
        	call MPI_ISSEND(this%f_(is-hd,js-hd,ks), 1, this%zPatch_, ptrMPIC%backNe_, &
        			   	   tag1, ptrMPIC%cartComm_, requests(4), ierror)       	
        	!send to front
        	call MPI_ISSEND(this%f_(is-hd,js-hd,ke-(hd-1)), 1, this%zPatch_, ptrMPIC%frontNe_, &
        			   	   tag2, ptrMPIC%cartComm_, requests(3), ierror) 
        			   	   
        	call MPI_WAITALL(4, requests, status, ierror)


        	CASE('sy')
        	
        	!******************************* right-left *********************************!
  			!recv from right
  			call MPI_IRECV(this%f_(ie+1,js-hd,ks-hd),1,this%xPatch_, ptrMPIC%rightNe_, &
  					       tag1, ptrMPIC%cartComm_, requests(1), ierror)
			!recv from left     		
  			call MPI_IRECV(this%f_(is-hd,js-hd,ks-hd),1,this%xPatch_, ptrMPIC%leftNe_, &
  					       tag2, ptrMPIC%cartComm_, requests(2), ierror) 
        	!send to left
        	call MPI_ISSEND(this%f_(is,js-hd,ks-hd), 1, this%xPatch_, ptrMPIC%leftNe_, &
        			   	   tag1, ptrMPIC%cartComm_, requests(4), ierror)        	
        	!send to right
        	call MPI_ISSEND(this%f_(ie-(hd-1),js-hd,ks-hd), 1, this%xPatch_, ptrMPIC%rightNe_, &
        			   	   tag2, ptrMPIC%cartComm_, requests(3), ierror) 
        			   	   
        	call MPI_WAITALL(4, requests, status, ierror)    	
        	               	
        	!******************************* top-bottom *********************************!
  			!recv from top
  			call MPI_IRECV(this%f_(is-hd,je+1,ks-hd),1,this%yPatch_, ptrMPIC%topNe_, &
  					       tag1, ptrMPIC%cartComm_, requests(1), ierror)
			!recv from bottom   		
  			call MPI_IRECV(this%f_(is-hd,js-hd,ks-hd),1,this%yPatch_, ptrMPIC%bottomNe_, &
  					       tag2, ptrMPIC%cartComm_, requests(2), ierror)  
        	!send to bottom
        	call MPI_ISSEND(this%f_(is-hd,js+1,ks-hd), 1, this%yPatch_, ptrMPIC%bottomNe_, &
        			   	   tag1, ptrMPIC%cartComm_, requests(4), ierror)       	
        	!send to top
        	call MPI_ISSEND(this%f_(is-hd,je-hd,ks-hd), 1, this%yPatch_, ptrMPIC%topNe_, &
        			   	   tag2, ptrMPIC%cartComm_, requests(3), ierror) 
        			   	   
        	call MPI_WAITALL(4, requests, status, ierror) 
        	      	
        	!******************************* front-back *********************************!
  			!recv from front
  			call MPI_IRECV(this%f_(is-hd,js-hd,ke+1),1,this%zPatch_, ptrMPIC%frontNe_, &
  					       tag1, ptrMPIC%cartComm_, requests(1), ierror)
			!recv from back   		
  			call MPI_IRECV(this%f_(is-hd,js-hd,ks-hd),1,this%zPatch_, ptrMPIC%backNe_, &
  					       tag2, ptrMPIC%cartComm_, requests(2), ierror)
        	!send to back
        	call MPI_ISSEND(this%f_(is-hd,js-hd,ks), 1, this%zPatch_, ptrMPIC%backNe_, &
        			   	   tag1, ptrMPIC%cartComm_, requests(4), ierror)        	
        	!send to front
        	call MPI_ISSEND(this%f_(is-hd,js-hd,ke-(hd-1)), 1, this%zPatch_, ptrMPIC%frontNe_, &
        			   	   tag2, ptrMPIC%cartComm_, requests(3), ierror)  
        			   	   
        	call MPI_WAITALL(4, requests, status, ierror) 
        	
        	 			  		
        	CASE('sz')
        	
        	!******************************* right-left *********************************!
  			!recv from right
  			call MPI_IRECV(this%f_(ie+1,js-hd,ks-hd),1,this%xPatch_, ptrMPIC%rightNe_, &
  					       tag1, ptrMPIC%cartComm_, requests(1), ierror)
			!recv from left     		
  			call MPI_IRECV(this%f_(is-hd,js-hd,ks-hd),1,this%xPatch_, ptrMPIC%leftNe_, &
  					       tag2, ptrMPIC%cartComm_, requests(2), ierror) 
        	!send to left
        	call MPI_ISSEND(this%f_(is,js-hd,ks-hd), 1, this%xPatch_, ptrMPIC%leftNe_, &
        			   	   tag1, ptrMPIC%cartComm_, requests(4), ierror)        	
        	!send to right
        	call MPI_ISSEND(this%f_(ie-(hd-1),js-hd,ks-hd), 1, this%xPatch_, ptrMPIC%rightNe_, &
        			   	   tag2, ptrMPIC%cartComm_, requests(3), ierror) 
        			   	   
        	call MPI_WAITALL(4, requests, status, ierror)     	
        	             	
        	!******************************* top-bottom *********************************!
  			!recv from top
  			call MPI_IRECV(this%f_(is-hd,je+1,ks-hd),1,this%yPatch_, ptrMPIC%topNe_, &
  					       tag1, ptrMPIC%cartComm_, requests(1), ierror)
			!recv from bottom   		
  			call MPI_IRECV(this%f_(is-hd,js-hd,ks-hd),1,this%yPatch_, ptrMPIC%bottomNe_, &
  					       tag2, ptrMPIC%cartComm_, requests(2), ierror) 
        	!send to bottom
        	call MPI_ISSEND(this%f_(is-hd,js,ks-hd), 1, this%yPatch_, ptrMPIC%bottomNe_, &
        			   	   tag1, ptrMPIC%cartComm_, requests(4), ierror)       	
        	!send to top
        	call MPI_ISSEND(this%f_(is-hd,je-(hd-1),ks-hd), 1, this%yPatch_, ptrMPIC%topNe_, &
        			   	   tag2, ptrMPIC%cartComm_, requests(3), ierror)  
        			   	   
        	call MPI_WAITALL(4, requests, status, ierror) 
        	
        	!******************************* front-back *********************************!
  			!recv from front
  			call MPI_IRECV(this%f_(is-hd,js-hd,ke+1),1,this%zPatch_, ptrMPIC%frontNe_, &
  					       tag1, ptrMPIC%cartComm_, requests(1), ierror)
			!recv from back   		
  			call MPI_IRECV(this%f_(is-hd,js-hd,ks-hd),1,this%zPatch_, ptrMPIC%backNe_, &
  					       tag2, ptrMPIC%cartComm_, requests(2), ierror)  
        	!send to back
        	call MPI_ISSEND(this%f_(is-hd,js-hd,ks+1), 1, this%zPatch_, ptrMPIC%backNe_, &
        			   	   tag1, ptrMPIC%cartComm_, requests(4), ierror)       	
        	!send to front
        	call MPI_ISSEND(this%f_(is-hd,js-hd,ke-hd), 1, this%zPatch_, ptrMPIC%frontNe_, &
        			   	   tag2, ptrMPIC%cartComm_, requests(3), ierror) 
        			   	   
        	call MPI_WAITALL(4, requests, status, ierror)  


        	CASE DEFAULT
        END SELECT
     

    end subroutine
!========================================================================================!

!========================================================================================!
   subroutine updateBoundaries(this)
        type(scalarField), intent(inout) :: this

		!update boundary fields
		call updateBoundaryField(this%bLeft_,this)
		call updateBoundaryField(this%bRight_,this)
		call updateBoundaryField(this%bBottom_,this)
		call updateBoundaryField(this%bTop_,this)
		call updateBoundaryField(this%bBack_,this)
		call updateBoundaryField(this%bFront_,this)

!DIR$ IF DEFINED (MG_MODE)		
		if (this%tp_ == 'cl') then
			!update external edges
			call updateEdges(this)
		
			!update external corners
			call updateCorners(this)
		end if
!DIR$ ENDIF
		
		!update halo
		call updateHalo(this)

    end subroutine
!========================================================================================!

!DIR$ IF DEFINED (MG_MODE)
!========================================================================================!
   subroutine updateEdges(this)
        type(scalarField), intent(inout) :: this
        integer :: i, nx, ny, nz
        type(grid), pointer :: ptrMesh
        real(DP) :: gx, gy, gz
        
        ptrMesh => this%ptrMesh_
	
		!check the 12 edges
		nx = this%nx_
		ny = this%ny_
		nz = this%nz_
        
        ! ALONG Y
        ! left-back
        if ( (this%bLeft_%isExternal_) .AND. (this%bBack_%isExternal_) ) then
           
           do i = 1,ny
           		gx = (this%f_(2,i,1)-this%f_(0,i,1))/(ptrMesh%dxc_(2)+ptrMesh%dxc_(1))
           		gz = (this%f_(1,i,2)-this%f_(1,i,0))/(ptrMesh%dzc_(2)+ptrMesh%dzc_(1))
           		this%f_(0,i,0)= this%f_(1,i,1) + gx*(ptrMesh%xc_(0)-ptrMesh%xc_(1)) + &
           						gz*(ptrMesh%zc_(0)-ptrMesh%zc_(1))
           end do
        
        end if
		
        ! right-back
        if ( (this%bRight_%isExternal_) .AND. (this%bBack_%isExternal_) ) then
           
           do i = 1,ny
           		gx = (this%f_(nx+1,i,1)-this%f_(nx-1,i,1))/(ptrMesh%dxc_(nx+1)+ptrMesh%dxc_(nx))
           		gz = (this%f_(nx,i,2)-this%f_(nx,i,0))/(ptrMesh%dzc_(2)+ptrMesh%dzc_(1))
           		this%f_(nx+1,i,0)= this%f_(nx,i,1) + gx*(ptrMesh%xc_(nx+1)-ptrMesh%xc_(nx)) + &
           						gz*(ptrMesh%zc_(0)-ptrMesh%zc_(1))
           end do
        
        end if
        
        ! right-front
        if ( (this%bRight_%isExternal_) .AND. (this%bFront_%isExternal_) ) then
           
           do i = 1,ny
           		gx = (this%f_(nx+1,i,nz)-this%f_(nx-1,i,nz))/(ptrMesh%dxc_(nx+1)+ptrMesh%dxc_(nx))
           		gz = (this%f_(nx,i,nz+1)-this%f_(nx,i,nz-1))/(ptrMesh%dzc_(nz+1)+ptrMesh%dzc_(nz))
           		this%f_(nx+1,i,nz+1)= this%f_(nx,i,nz) + gx*(ptrMesh%xc_(nx+1)-ptrMesh%xc_(nx)) + &
           						gz*(ptrMesh%zc_(nz+1)-ptrMesh%zc_(nz))
           end do
        
        end if
        
        ! left-front
        if( (this%bLeft_%isExternal_) .AND. (this%bFront_%isExternal_) ) then
           
           do i = 1,ny
           		gx = (this%f_(2,i,nz)-this%f_(0,i,nz))/(ptrMesh%dxc_(2)+ptrMesh%dxc_(1))
           		gz = (this%f_(1,i,nz+1)-this%f_(1,i,nz-1))/(ptrMesh%dzc_(nz+1)+ptrMesh%dzc_(nz))
           		this%f_(0,i,nz+1)= this%f_(1,i,nz) + gx*(ptrMesh%xc_(0)-ptrMesh%xc_(1)) + &
           						gz*(ptrMesh%zc_(nz+1)-ptrMesh%zc_(nz))
           end do
        
        end if
        
        ! ALONG X
        ! back-top
        if( (this%bBack_%isExternal_) .AND. (this%bTop_%isExternal_) ) then
           
           do i = 1,nx
           		gy = (this%f_(i,ny+1,1)-this%f_(i,ny-1,1))/(ptrMesh%dyc_(ny+1)+ptrMesh%dyc_(ny))
           		gz = (this%f_(i,ny,2)-this%f_(i,ny,0))/(ptrMesh%dzc_(2)+ptrMesh%dzc_(1))
           		this%f_(i,ny+1,0)= this%f_(i,ny,1) + gy*(ptrMesh%yc_(ny+1)-ptrMesh%yc_(ny)) + &
           						gz*(ptrMesh%zc_(0)-ptrMesh%zc_(1))
           end do
        
        end if

        ! back-bottom
        if( (this%bBack_%isExternal_) .AND. (this%bBottom_%isExternal_) ) then
           
           do i = 1,nx
           		gy = (this%f_(i,2,1)-this%f_(i,0,1))/(ptrMesh%dyc_(2)+ptrMesh%dyc_(1))
           		gz = (this%f_(i,1,2)-this%f_(i,1,0))/(ptrMesh%dzc_(2)+ptrMesh%dzc_(1))
           		this%f_(i,0,0)= this%f_(i,1,1) + gy*(ptrMesh%yc_(0)-ptrMesh%yc_(1)) + &
           						gz*(ptrMesh%zc_(0)-ptrMesh%zc_(1))
           end do
        
        end if
        
        ! front-bottom
        if( (this%bFront_%isExternal_) .AND. (this%bBottom_%isExternal_) ) then
           
           do i = 1,nx
           		gy = (this%f_(i,2,nz)-this%f_(i,0,nz))/(ptrMesh%dyc_(2)+ptrMesh%dyc_(1))
           		gz = (this%f_(i,1,nz+1)-this%f_(i,1,nz-1))/(ptrMesh%dzc_(nz+1)+ptrMesh%dzc_(nz))
           		this%f_(i,0,nz+1)= this%f_(i,1,nz) + gy*(ptrMesh%yc_(0)-ptrMesh%yc_(1)) + &
           						gz*(ptrMesh%zc_(nz+1)-ptrMesh%zc_(nz))
           end do
        
        end if
        
        ! front-top
        if( (this%bFront_%isExternal_) .AND. (this%bTop_%isExternal_) ) then
           
           do i = 1,nx
           		gy = (this%f_(i,ny+1,nz)-this%f_(i,ny-1,nz))/(ptrMesh%dyc_(ny+1)+ptrMesh%dyc_(ny))
           		gz = (this%f_(i,ny,nz+1)-this%f_(i,ny,nz-1))/(ptrMesh%dzc_(nz+1)+ptrMesh%dzc_(nz))
           		this%f_(i,ny+1,nz+1)= this%f_(i,ny,nz) + gy*(ptrMesh%yc_(ny+1)-ptrMesh%yc_(ny)) + &
           						gz*(ptrMesh%zc_(nz+1)-ptrMesh%zc_(nz))
           end do
        
        end if
        
        ! ALONG z
        ! top-left
        if( (this%bTop_%isExternal_) .AND. (this%bLeft_%isExternal_) ) then
           
           do i = 1,nz
           		gx = (this%f_(2,ny,i)-this%f_(0,ny,i))/(ptrMesh%dxc_(2)+ptrMesh%dxc_(1))
           		gy = (this%f_(1,ny+1,i)-this%f_(1,ny-1,i))/(ptrMesh%dyc_(ny+1)+ptrMesh%dyc_(ny))
           		this%f_(0,ny+1,i)= this%f_(1,ny,i) + gx*(ptrMesh%xc_(0)-ptrMesh%xc_(1)) + &
           						gy*(ptrMesh%yc_(ny+1)-ptrMesh%yc_(ny))
           end do
        
        end if
        
        ! top-right
        if( (this%bTop_%isExternal_) .AND. (this%bRight_%isExternal_) ) then
           
           do i = 1,nz
           		gx = (this%f_(nx+1,ny,i)-this%f_(nx-1,ny,i))/(ptrMesh%dxc_(nx+1)+ptrMesh%dxc_(nx))
           		gy = (this%f_(nx,ny+1,i)-this%f_(nx,ny-1,i))/(ptrMesh%dyc_(ny+1)+ptrMesh%dyc_(ny))
           		this%f_(nx+1,ny+1,i)= this%f_(nx,ny,i) + gx*(ptrMesh%xc_(nx+1)-ptrMesh%xc_(nx)) + &
           						gy*(ptrMesh%yc_(ny+1)-ptrMesh%yc_(ny))
           end do
        
        end if
        
        ! bottom-right
        if( (this%bBottom_%isExternal_) .AND. (this%bRight_%isExternal_) ) then
           
           do i = 1,nz
           		gx = (this%f_(nx+1,1,i)-this%f_(nx-1,1,i))/(ptrMesh%dxc_(nx+1)+ptrMesh%dxc_(nx))
           		gy = (this%f_(1,2,i)-this%f_(1,0,i))/(ptrMesh%dyc_(2)+ptrMesh%dyc_(1))
           		this%f_(nx+1,0,i)= this%f_(nx,1,i) + gx*(ptrMesh%xc_(nx+1)-ptrMesh%xc_(nx)) + &
           						gy*(ptrMesh%yc_(0)-ptrMesh%yc_(1))
           end do
        
        end if
        
        ! bottom-left
        if( (this%bBottom_%isExternal_) .AND. (this%bLeft_%isExternal_) ) then
           
           do i = 1,nz
           		gx = (this%f_(2,1,i)-this%f_(0,1,i))/(ptrMesh%dxc_(2)+ptrMesh%dxc_(1))
           		gy = (this%f_(1,2,i)-this%f_(1,0,i))/(ptrMesh%dyc_(2)+ptrMesh%dyc_(1))
           		this%f_(0,0,i)= this%f_(1,1,i) + gx*(ptrMesh%xc_(0)-ptrMesh%xc_(1)) + &
           						gy*(ptrMesh%yc_(0)-ptrMesh%yc_(1))
           end do

        
        end if		

		
    end subroutine
!========================================================================================!

!========================================================================================!
   subroutine updateCorners(this)
        type(scalarField), intent(inout) :: this
        integer :: nx, ny, nz
        type(grid), pointer :: ptrMesh
        real(DP) :: gx, gy, gz

        ptrMesh => this%ptrMesh_
	
		!check 8 external corners
		nx = this%nx_
		ny = this%ny_
		nz = this%nz_  
		
		! top - back - left
		if( (this%bTop_%isExternal_) &
			.AND. (this%bBack_%isExternal_) &
			.AND. (this%bLeft_%isExternal_)) then      
			
			gx = (this%f_(2,ny,1)-this%f_(0,ny,1))/(ptrMesh%dxc_(2)+ptrMesh%dxc_(1))
			gy = (this%f_(1,ny+1,1)-this%f_(1,ny-1,1))/(ptrMesh%dyc_(ny+1)+ptrMesh%dyc_(ny))
			gz = (this%f_(1,ny,2)-this%f_(1,ny,0))/(ptrMesh%dzc_(2)+ptrMesh%dzc_(1))
           	this%f_(0,ny+1,0)= this%f_(1,ny,1) + gx*(ptrMesh%xc_(0)-ptrMesh%xc_(1)) + &
           										 gy*(ptrMesh%yc_(ny+1)-ptrMesh%yc_(ny)) + &
           										 gz*(ptrMesh%zc_(0)-ptrMesh%zc_(1))
			
		end if
        

		! top - back - right
		if( (this%bTop_%isExternal_) &
			.AND. (this%bBack_%isExternal_) &
			.AND. (this%bRight_%isExternal_)) then      
			
			gx = (this%f_(nx+1,ny,1)-this%f_(nx-1,ny,1))/(ptrMesh%dxc_(nx+1)+ptrMesh%dxc_(nx))
			gy = (this%f_(nx,ny+1,1)-this%f_(nx,ny-1,1))/(ptrMesh%dyc_(ny+1)+ptrMesh%dyc_(ny))
			gz = (this%f_(nx,ny,2)-this%f_(nx,ny,0))/(ptrMesh%dzc_(2)+ptrMesh%dzc_(1))
           	this%f_(nx+1,ny+1,0)= this%f_(nx,ny,1) + gx*(ptrMesh%xc_(nx+1)-ptrMesh%xc_(nx)) + &
           										 gy*(ptrMesh%yc_(ny+1)-ptrMesh%yc_(ny)) + &
           										 gz*(ptrMesh%zc_(0)-ptrMesh%zc_(1))
			
		end if
		

		! top - front - right
		if( (this%bTop_%isExternal_) &
			.AND. (this%bFront_%isExternal_) &
			.AND. (this%bRight_%isExternal_)) then      
			
			gx = (this%f_(nx+1,ny,nz)-this%f_(nx-1,ny,nz))/(ptrMesh%dxc_(nx+1)+ptrMesh%dxc_(nx))
			gy = (this%f_(nx,ny+1,nz)-this%f_(nx,ny-1,nz))/(ptrMesh%dyc_(ny+1)+ptrMesh%dyc_(ny))
			gz = (this%f_(nx,ny,nz+1)-this%f_(nx,ny,nz-1))/(ptrMesh%dzc_(nz+1)+ptrMesh%dzc_(nz))
           	this%f_(nx+1,ny+1,nz+1)= this%f_(nx,ny,nz) + gx*(ptrMesh%xc_(nx+1)-ptrMesh%xc_(nx)) + &
           										 gy*(ptrMesh%yc_(ny+1)-ptrMesh%yc_(ny)) + &
           										 gz*(ptrMesh%zc_(nz+1)-ptrMesh%zc_(nz))
			
		end if
		

		! top - front - left
		if( (this%bTop_%isExternal_) &
			.AND. (this%bFront_%isExternal_) &
			.AND. (this%bLeft_%isExternal_)) then      
			
			gx = (this%f_(2,ny,nz)-this%f_(0,ny,nz))/(ptrMesh%dxc_(2)+ptrMesh%dxc_(1))
			gy = (this%f_(1,ny+1,nz)-this%f_(1,ny-1,nz))/(ptrMesh%dyc_(ny+1)+ptrMesh%dyc_(ny))
			gz = (this%f_(1,ny,nz+1)-this%f_(1,ny,nz-1))/(ptrMesh%dzc_(nz+1)+ptrMesh%dzc_(nz))
           	this%f_(0,ny+1,nz+1)= this%f_(1,ny,nz) + gx*(ptrMesh%xc_(0)-ptrMesh%xc_(1)) + &
           										 gy*(ptrMesh%yc_(ny+1)-ptrMesh%yc_(ny)) + &
           										 gz*(ptrMesh%zc_(nz+1)-ptrMesh%zc_(nz))
			
		end if
		
	
		! bottom - back - left
		if( (this%bBottom_%isExternal_) &
			.AND. (this%bBack_%isExternal_) &
			.AND. (this%bLeft_%isExternal_)) then      
			
			gx = (this%f_(2,1,1)-this%f_(0,1,1))/(ptrMesh%dxc_(2)+ptrMesh%dxc_(1))
			gy = (this%f_(1,2,1)-this%f_(1,0,1))/(ptrMesh%dyc_(2)+ptrMesh%dyc_(1))
			gz = (this%f_(1,1,2)-this%f_(1,1,0))/(ptrMesh%dzc_(2)+ptrMesh%dzc_(1))
           	this%f_(0,0,0)= this%f_(1,1,1) + gx*(ptrMesh%xc_(0)-ptrMesh%xc_(1)) + &
           										 gy*(ptrMesh%yc_(0)-ptrMesh%yc_(1)) + &
           										 gz*(ptrMesh%zc_(0)-ptrMesh%zc_(1))
			
		end if
		
		
		! bottom - back - right
		if( (this%bBottom_%isExternal_) &
			.AND. (this%bBack_%isExternal_) &
			.AND. (this%bRight_%isExternal_)) then      
			
			gx = (this%f_(nx+1,1,1)-this%f_(nx-1,1,1))/(ptrMesh%dxc_(nx+1)+ptrMesh%dxc_(nx))
			gy = (this%f_(nx,2,1)-this%f_(nx,0,1))/(ptrMesh%dyc_(2)+ptrMesh%dyc_(1))
			gz = (this%f_(nx,1,2)-this%f_(nx,1,0))/(ptrMesh%dzc_(2)+ptrMesh%dzc_(1))
           	this%f_(nx+1,0,0)= this%f_(nx,1,1) + gx*(ptrMesh%xc_(nx+1)-ptrMesh%xc_(nx)) + &
           										 gy*(ptrMesh%yc_(0)-ptrMesh%yc_(1)) + &
           										 gz*(ptrMesh%zc_(0)-ptrMesh%zc_(1))
			
		end if
		
		
		! bottom - front - right
		if( (this%bBottom_%isExternal_) &
			.AND. (this%bFront_%isExternal_) &
			.AND. (this%bRight_%isExternal_)) then      
			
			gx = (this%f_(nx+1,1,nz)-this%f_(nx-1,1,nz))/(ptrMesh%dxc_(nx+1)+ptrMesh%dxc_(nx))
			gy = (this%f_(nx,2,nz)-this%f_(nx,0,nz))/(ptrMesh%dyc_(2)+ptrMesh%dyc_(1))
			gz = (this%f_(nx,1,nz+1)-this%f_(nx,1,nz-1))/(ptrMesh%dzc_(nz+1)+ptrMesh%dzc_(nz))
           	this%f_(nx+1,0,nz+1)= this%f_(nx,1,nz) + gx*(ptrMesh%xc_(nx+1)-ptrMesh%xc_(nx)) + &
           										 gy*(ptrMesh%yc_(0)-ptrMesh%yc_(1)) + &
           										 gz*(ptrMesh%zc_(nz+1)-ptrMesh%zc_(nz))
			
		end if
		
		
		! bottom - front - left
		if( (this%bBottom_%isExternal_) &
			.AND. (this%bFront_%isExternal_) &
			.AND. (this%bLeft_%isExternal_)) then      
			
			gx = (this%f_(2,1,nz)-this%f_(0,1,nz))/(ptrMesh%dxc_(2)+ptrMesh%dxc_(1))
			gy = (this%f_(1,2,nz)-this%f_(1,0,nz))/(ptrMesh%dyc_(2)+ptrMesh%dyc_(1))
			gz = (this%f_(1,1,nz+1)-this%f_(1,1,nz-1))/(ptrMesh%dzc_(nz+1)+ptrMesh%dzc_(nz))
           	this%f_(0,0,nz+1)= this%f_(1,1,nz) + gx*(ptrMesh%xc_(0)-ptrMesh%xc_(1)) + &
           										 gy*(ptrMesh%yc_(0)-ptrMesh%yc_(1)) + &
           										 gz*(ptrMesh%zc_(nz+1)-ptrMesh%zc_(nz))
			
		end if


        
    end subroutine
!========================================================================================!

!========================================================================================!
   subroutine coarsenField(this)
        type(scalarField), intent(inout) :: this	!fine field at the current grid level
        integer :: nx, ny, nz
		
		!allocate coarse field
		call allocatePtrf(this)
		
		!allocate prol-field
		call allocateArray(this%prol_,1,this%nx_,1,this%ny_,1,this%nz_)
		
		!check for pointer association to coarse grid
		if (.not. associated(this%ptrMesh_%ptrGrid_)) then
			call mpiABORT('Attempt to coarsen field before coarsening grid ')
		end if
		
		!init base class ioFile
		this%ptrf_%fileName_ = this%fileName_//'_c'
		this%ptrf_%filePath_ = s_fileDir//'/'//this%fileName_
		
		!associate coarse mesh
		this%ptrf_%ptrMesh_ => this%ptrMesh_%ptrGrid_
		
		!dimensions
		nx = this%ptrf_%ptrMesh_%nx_
		ny = this%ptrf_%ptrMesh_%ny_
		nz = this%ptrf_%ptrMesh_%nz_
		this%ptrf_%nx_ = nx
		this%ptrf_%ny_ = ny
		this%ptrf_%nz_ = nz
		
		!field type
		this%ptrf_%tp_ = this%tp_
		
		!halo dim
		this%ptrf_%hd_ = this%hd_
		
		!coarse indexes
		call setIndexesBounds(this%ptrf_,this%tp_)
		
		call allocateArray(this%ptrf_%f_,											&
							this%ptrf_%is_-this%hd_,this%ptrf_%ie_+this%hd_,&
							this%ptrf_%js_-this%hd_,this%ptrf_%je_+this%hd_,&
							this%ptrf_%ks_-this%hd_,this%ptrf_%ke_+this%hd_)
							
		call initToZero(this%ptrf_)
		
		!build halo exchange types	
		call buildHaloExchangeTypes(this%ptrf_)		
		
		!coarsen boundary fields
		call coarsenBoundary(this%bLeft_,this%ptrf_%bLeft_,		&
										 this%ptrf_%ptrMesh_,		&
										 this%ptrf_%tp_)
										 
		call coarsenBoundary(this%bRight_,this%ptrf_%bRight_,		&
										 this%ptrf_%ptrMesh_,		&
										 this%ptrf_%tp_)
										 
		call coarsenBoundary(this%bBottom_,this%ptrf_%bBottom_,	&
										 this%ptrf_%ptrMesh_,		&
										 this%ptrf_%tp_)
										 
		call coarsenBoundary(this%bTop_,this%ptrf_%bTop_,			&
										 this%ptrf_%ptrMesh_,	&
										 this%ptrf_%tp_)
										 
		call coarsenBoundary(this%bBack_,this%ptrf_%bBack_,		&
										 this%ptrf_%ptrMesh_,		&
										 this%ptrf_%tp_)
										 
		call coarsenBoundary(this%bFront_,this%ptrf_%bFront_,		&
										 this%ptrf_%ptrMesh_,		&
										 this%ptrf_%tp_)
		

		
    end subroutine
!========================================================================================!

!========================================================================================!
    recursive subroutine coarsenFields(this,n)
        type(scalarField), intent(inout) :: this
        integer, intent(in) :: n
		
		if (n > 1) then 
			call coarsenField(this)
			call coarsenFields(this%ptrf_,n-1)
		else
			return
		end if
		
        
    end subroutine
!========================================================================================!

!========================================================================================!
    recursive subroutine resetBCerrorField(this)
        type(scalarField), intent(inout) :: this
		
		if ( associated(this%ptrf_) ) then
			
			!bLeft
			if ( (this%ptrf_%bLeft_%bType_ == s_fixedValue)     &
				.OR. 								      			  &
			   (this%ptrf_%bLeft_%bType_ == s_normalGradient) ) then
				
				this%ptrf_%bLeft_%bf_ = 0.d0
			end if
			
			!bRight
			if ( (this%ptrf_%bRight_%bType_ == s_fixedValue)     &
				.OR. 								      			   &
			   (this%ptrf_%bRight_%bType_ == s_normalGradient) ) then
				
				this%ptrf_%bRight_%bf_ = 0.d0
			end if
			
			!bBottom
			if ( (this%ptrf_%bBottom_%bType_ == s_fixedValue)     &
				.OR. 								      			    &
			   (this%ptrf_%bBottom_%bType_ == s_normalGradient) ) then
				
				this%ptrf_%bBottom_%bf_ = 0.d0
			end if
			
			!bTop
			if ( (this%ptrf_%bTop_%bType_ == s_fixedValue)     &
				.OR. 								      			 &
			   (this%ptrf_%bTop_%bType_ == s_normalGradient) ) then
				
				this%ptrf_%bTop_%bf_ = 0.d0
			end if
			
			!bBack
			if ( (this%ptrf_%bBack_%bType_ == s_fixedValue)     &
				.OR. 								      			  &
			   (this%ptrf_%bBack_%bType_ == s_normalGradient) ) then
				
				this%ptrf_%bBack_%bf_ = 0.d0
			end if
			
			!bFront
			if ( (this%ptrf_%bFront_%bType_ == s_fixedValue)     &
				.OR. 								     			   &
			   (this%ptrf_%bFront_%bType_ == s_normalGradient) ) then
				
				this%ptrf_%bFront_%bf_ = 0.d0
			end if
			
			call resetBCerrorField(this%ptrf_)
			
		else
			return
		end if

    end subroutine
!========================================================================================!

!========================================================================================!
	subroutine restrict_scalarField(this) 
		type(scalarField), intent(inout) :: this
		integer :: nxf, nyf, nzf, nxc, nyc, nzc

		!fine mesh size
		nxf = this%nx_
		nyf = this%ny_
		nzf = this%nz_
			
		!coarse mesh size
		if (associated(this%ptrf_)) then
			nxc = this%ptrf_%nx_
			nyc = this%ptrf_%ny_
			nzc = this%ptrf_%nz_
		else
			call mpiABORT('Attempt to restrict on not-allocated scalarField pointer ')
		end if
		
		call restrictionOp3D(this%f_(1:nxf,1:nyf,1:nzf),				&
							this%ptrf_%f_(1:nxc,1:nyc,1:nzc),		&
							this%ptrMesh_%xc_(1:nxf),					&
							this%ptrMesh_%yc_(1:nyf),					&
							this%ptrMesh_%zc_(1:nzf),					&
							this%ptrMesh_%ptrGrid_%xc_(1:nxc),			&
							this%ptrMesh_%ptrGrid_%yc_(1:nyc),			&
							this%ptrMesh_%ptrGrid_%zc_(1:nzc))		
							
		call updateBoundaries(this%ptrf_)
	
	end subroutine
!========================================================================================!

!========================================================================================!
	subroutine restrict_extField_scalarField(this,f) 
		type(scalarField), intent(inout) :: this
		real(DP), allocatable, dimension(:,:,:) :: f
		integer :: nxf, nyf, nzf, nxc, nyc, nzc

		!fine mesh size
		nxf = this%nx_
		nyf = this%ny_
		nzf = this%nz_
			
		!coarse mesh size
		if (associated(this%ptrf_)) then
			nxc = this%ptrf_%nx_
			nyc = this%ptrf_%ny_
			nzc = this%ptrf_%nz_
		else
			call mpiABORT('Attempt to restrict on not-allocated scalarField pointer ')
		end if
		
		call restrictionOp3D(f,											&
							this%ptrf_%f_(1:nxc,1:nyc,1:nzc),		&
							this%ptrMesh_%xc_(1:nxf),					&
							this%ptrMesh_%yc_(1:nyf),					&
							this%ptrMesh_%zc_(1:nzf),					&
							this%ptrMesh_%ptrGrid_%xc_(1:nxc),			&
							this%ptrMesh_%ptrGrid_%yc_(1:nyc),			&
							this%ptrMesh_%ptrGrid_%zc_(1:nzc))		
							
		call updateBoundaries(this%ptrf_)
	
	end subroutine
!========================================================================================!

!========================================================================================!
	subroutine prolongateField(this,x) 
		type(scalarField), intent(inout) :: this
		real(DP), allocatable, dimension(:,:,:), intent(inout) :: x
		integer :: nxf, nyf, nzf

		!fine mesh size
		nxf = this%nx_
		nyf = this%ny_
		nzf = this%nz_
		
		call prolongationOp3D(x,									&
							  this%ptrf_%f_,					&
							  this%ptrMesh_%xc_(1:nxf),				&
							  this%ptrMesh_%yc_(1:nyf),				&
							  this%ptrMesh_%zc_(1:nzf),				&
							  this%ptrMesh_%ptrGrid_%xc_,			&
							  this%ptrMesh_%ptrGrid_%yc_,			&
							  this%ptrMesh_%ptrGrid_%zc_)	
			
	end subroutine
!========================================================================================!

!DIR$ ENDIF

!========================================================================================!
	subroutine allocatePtrf(this) 
		type(scalarField), intent(inout) :: this
		integer :: err
		
		
		if (.not. associated(this%ptrf_)) then
			
			allocate(this%ptrf_,STAT=err)

			if (err /= 0) then
				call mpiABORT('Allocation of ptrScField failed ') 
			end if
		else
			call mpiABORT('Attempt to allocate an associated ptrScField ')
		end if

		
	end subroutine
!========================================================================================!

!========================================================================================!
	subroutine deallocatePtrf(this) 
		type(scalarField), intent(inout) :: this

        if (associated(this%ptrf_)) then
            deallocate(this%ptrf_)
        end if
	
	end subroutine
!========================================================================================!

!========================================================================================!
	subroutine reconstructAndWriteField(this,gf,nFolder,write_to_file)
		type(scalarField), intent(inout) :: this
		type(scalarField), intent(inout) :: gf
		integer, intent(in) :: nFolder
		logical, intent(in), optional :: write_to_file
		logical :: opt_write_to_file
        integer :: ierror
        integer, dimension(MPI_STATUS_SIZE) :: status
        integer :: i0, i1, j0, j1, k0, k1, n
        integer :: i, j, k
        integer :: il, jl, kl
        integer :: tag = 0
        type(mpiControl), pointer :: ptrMPIC
        

		if (present(write_to_file)) then
			opt_write_to_file = write_to_file
		else
			opt_write_to_file = .TRUE.
		end if        

        
        ptrMPIC => this%ptrMesh_%ptrMPIC_
        
!DIR$ IF DEFINED (MEM_SAVE)
		if (IS_MASTER) then
			call reAllocateArray(gf%f_,gf%is_-this%hd_,gf%ie_+this%hd_,	&
						   	     gf%js_-this%hd_,gf%je_+this%hd_,		&
						   	     gf%ks_-this%hd_,gf%ke_+this%hd_)
			call initToZero(gf)
		end if
!DIR$ ENDIF
        
        
		!reconstruct domain on the master proc first
		if (IS_MASTER) then

			!set global indexes proc 0
			call globalIndexesFromAll(this%ptrMesh_,i0,i1,j0,j1,k0,k1,0,this%tp_) 
			
			il = 0
			jl = 0
			kl = 0	
			do k=k0,k1
				do j=j0,j1
					do i=i0,i1
						gf%f_(i,j,k) = this%f_(this%is_ + il, this%js_ + jl, this%ks_ + kl)
						il = il +1
					end do
					il = 0
					jl = jl +1
				end do
				il = 0
				jl = 0
				kl = kl + 1
			end do
			
		end if	
        
		!reconstruct domain and send to master
		if (IS_MASTER) then
			do n=1,ptrMPIC%nProcs_-1
			
				call MPI_Recv(this%f_,size(this%f_),MPI_DOUBLE_PRECISION,n,tag,ptrMPIC%cartComm_,status,ierror)
			
				!set global indexes proc n
				call globalIndexesFromAll(this%ptrMesh_,i0,i1,j0,j1,k0,k1,n,this%tp_) 
				
				il = 0
				jl = 0
				kl = 0	
				do k=k0,k1
					do j=j0,j1
						do i=i0,i1
							gf%f_(i,j,k) = this%f_(this%is_ + il, this%js_ + jl, this%ks_ + kl)
							il = il +1
						end do
						il = 0
						jl = jl +1
					end do
					il = 0
					jl = 0
					kl = kl + 1
				end do
				
			end do
		
		else
		
			!note: only the internal field is being reconstructed 
			!(no need to update halo since the original global field is only a placeholder for IN-OUT operations)
			call MPI_Send(this%f_,size(this%f_),MPI_DOUBLE_PRECISION,0,tag,ptrMPIC%cartComm_,ierror)
			
		end if
        
        
        !reset field on the master proc previously overwritten for send-recv MPI
		if (IS_MASTER) then

			!set global indexes proc 0
			call globalIndexesFromAll(this%ptrMesh_,i0,i1,j0,j1,k0,k1,0,this%tp_) 
			
			il = 0
			jl = 0
			kl = 0			
			do k=k0,k1
				do j=j0,j1
					do i=i0,i1
						this%f_(this%is_ + il, this%js_ + jl, this%ks_ + kl) = gf%f_(i,j,k) 
						il = il +1
					end do
					il = 0
					jl = jl +1
				end do
				il = 0
				jl = 0
				kl = kl + 1
			end do
			
		end if
		
		!this avoid a temporary copy
		call updateBoundaries(this)
		
		if (opt_write_to_file) then
			call writeField(gf,nFolder)
		end if
		
!DIR$ IF DEFINED (MEM_SAVE)
		if (IS_MASTER) then
			call deallocateArray(gf%f_)
		end if
!DIR$ ENDIF
	

	end subroutine
!========================================================================================!

!========================================================================================!
	subroutine setIndexesBounds(this,tp) 
		type(scalarField), intent(inout) :: this
		character(len=2), intent(in) :: tp
		integer :: nx, ny, nz
		
		nx = this%ptrMesh_%nx_
		ny = this%ptrMesh_%ny_
		nz = this%ptrMesh_%nz_
		
		SELECT CASE (tp)
			CASE('cl')
				this%is_ = 1
				this%ie_ = nx
				this%js_ = 1
				this%je_ = ny
				this%ks_ = 1
				this%ke_ = nz
			CASE('sx')
				this%is_ = 0
				this%ie_ = nx
				this%js_ = 1
				this%je_ = ny
				this%ks_ = 1
				this%ke_ = nz
			CASE('sy')
				this%is_ = 1
				this%ie_ = nx
				this%js_ = 0
				this%je_ = ny
				this%ks_ = 1
				this%ke_ = nz
			CASE('sz')
				this%is_ = 1
				this%ie_ = nx
				this%js_ = 1
				this%je_ = ny
				this%ks_ = 0
				this%ke_ = nz
   			CASE DEFAULT
   				call mpiABORT('Invalid scalarField type ')
			END SELECT
		
			
	end subroutine
!========================================================================================!

!========================================================================================!
	subroutine copyBoundary(cpf,f)
		type(scalarField), intent(inout) :: cpf
		type(scalarField), intent(in) :: f
		
		call buildHaloExchangeTypes(cpf)
		
		cpf%bRight_ = f%bRight_
		cpf%bLeft_ = f%bLeft_
		cpf%bBottom_ = f%bBottom_
		cpf%bTop_ = f%bTop_
		cpf%bBack_ = f%bBack_
		cpf%bFront_ = f%bFront_
		
		
	end subroutine
!========================================================================================!





	
end module scalarFieldMod


	



