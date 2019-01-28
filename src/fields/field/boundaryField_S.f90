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

!========================================================================================!
	subroutine boundaryFieldCTOR(this,f,bNumber,readBC,bType,bValue)
		type(boundaryField) :: this
		type(field), intent(in) :: f
		integer, intent(in) :: bNumber
		logical, intent(in) :: readBC
		integer, intent(in), optional :: bType
		real(DP), intent(in), optional :: bValue

		if (readBC) then
			call readBoundaryField(this,f,bNumber,bType,bValue)
		else
			call initDeafaultBoundaryField(this,f,bNumber)
		end if

	end subroutine
!========================================================================================!

!========================================================================================!
    subroutine initDeafaultBoundaryField(this,f,bNumber)
        type(boundaryField), intent(inout) :: this
        type(field), intent(in) :: f
        integer, intent(in) :: bNumber
        type(mpiControl), pointer :: ptrMPIC
        logical, dimension(3) :: wrap
        	
		ptrMPIC => f%ptrMesh_%ptrMPIC_
		wrap = ptrMPIC%wrapAround_	
		
		
		SELECT CASE (bNumber)
   			CASE (1)
  				if (wrap(1)) then
					call initBoundary(this,f,bNumber,s_periodicBC)
					this%bf_=0.d0
				else
					call initBoundary(this,f,bNumber,s_normalGradient)
					this%bf_=0.d0			
				end if 				
   			CASE (2) 
  				if (wrap(1)) then
					call initBoundary(this,f,bNumber,s_periodicBC)
					this%bf_=0.d0
				else
					call initBoundary(this,f,bNumber,s_normalGradient)
					this%bf_=0.d0			
				end if 
   			CASE (3) 
  				if (wrap(2)) then
					call initBoundary(this,f,bNumber,s_periodicBC)
					this%bf_=0.d0
				else
					call initBoundary(this,f,bNumber,s_normalGradient)
					this%bf_=0.d0			
				end if 
   			CASE (4)
  				if (wrap(2)) then
					call initBoundary(this,f,bNumber,s_periodicBC)
					this%bf_=0.d0
				else
					call initBoundary(this,f,bNumber,s_normalGradient)
					this%bf_=0.d0			
				end if  
   			CASE (5) 
  				if (wrap(3)) then
					call initBoundary(this,f,bNumber,s_periodicBC)
					this%bf_=0.d0
				else
					call initBoundary(this,f,bNumber,s_normalGradient)
					this%bf_=0.d0			
				end if 
   			CASE (6) 
  				if (wrap(3)) then
					call initBoundary(this,f,bNumber,s_periodicBC)
					this%bf_=0.d0
				else
					call initBoundary(this,f,bNumber,s_normalGradient)
					this%bf_=0.d0			
				end if 
   			CASE DEFAULT
		END SELECT
			
		!set metrics
		call setMetrics(this,f%ptrMesh_	,f%tp_)
			        
    end subroutine
!========================================================================================!

!========================================================================================!
    subroutine readBoundaryField(this,f,bNumber,bType,bValue)
        type(boundaryField), intent(inout) :: this
        type(field), intent(in) :: f
        integer, intent(in) :: bNumber,bType
        real(DP), intent(in) :: bValue
        integer :: bn
			
		if ( (bType == s_fixedValue)     &
			 .OR. 					  	 &
			 (bType == s_normalGradient) &
			 .OR.						 &
			 (bType == s_periodicBC) 	 &
			 .OR.						 &
			 (bType == s_calculatedBC) 	 &
			 .OR.						 &
			 (bType == s_contactAngleBC) &
		   ) then
		   call initBoundary(this,f,bNumber,bType)
		   call initField(this,f,bValue)
		else			
			call mpiABORT('Invalid boundary type for field '//f%fileName_)
		end if
			
		!set metrics
		call setMetrics(this,f%ptrMesh_	,f%tp_)
			        
    end subroutine
!========================================================================================!

!========================================================================================!
    subroutine setMetrics(this,mesh,tp)
        type(boundaryField), intent(inout) :: this
		type(grid), intent(in) :: mesh
		character(len=2), intent(in) :: tp
		integer :: nx, ny, nz
		
		
		nx = mesh%nx_
		ny = mesh%ny_
		nz = mesh%nz_
		
		!default case = 'cl'
		SELECT CASE (this%cartDir_)
   			CASE (-1) !left boundary
				this%dhi_ = mesh%dxc_(1)
				this%dhb_ = mesh%xf_(0) - mesh%xc_(0)
   			CASE (1) !right boundary
				this%dhi_ = mesh%dxc_(nx+1)
				this%dhb_ = mesh%xc_(nx+1) - mesh%xf_(nx)
   			CASE (-2) !bottom boundary
				this%dhi_ = mesh%dyc_(1)
				this%dhb_ = mesh%yf_(0) - mesh%yc_(0)
   			CASE (2) !top boundary
				this%dhi_ = mesh%dyc_(ny+1)
				this%dhb_ = mesh%yc_(ny+1) - mesh%yf_(ny)
   			CASE (-3) !back boundary
				this%dhi_ = mesh%dzc_(1)
				this%dhb_ = mesh%zf_(0) - mesh%zc_(0)
   			CASE (3) !front boundary
				this%dhi_ = mesh%dzc_(nz+1)
				this%dhb_ = mesh%zc_(nz+1) - mesh%zf_(nz)
   			CASE DEFAULT
		END SELECT
		
		SELECT CASE (tp)
			
			CASE('sx')
			
				SELECT CASE (this%cartDir_)
   					CASE (-1) !left boundary
						this%dhi_ = mesh%dxf_(0) + mesh%dxf_(1)
						this%dhb_ = mesh%dxf_(0)
   					CASE (1) !right boundary
						this%dhi_ = mesh%dxf_(nx) + mesh%dxf_(nx+1)
						this%dhb_ = mesh%dxf_(nx+1)
   					CASE DEFAULT
				END SELECT
			
			CASE('sy')
			
				SELECT CASE (this%cartDir_)
   					CASE (-2) !bottom boundary
						this%dhi_ = mesh%dyf_(0) + mesh%dyf_(1)
						this%dhb_ = mesh%dyf_(0)
   					CASE (2) !top boundary
						this%dhi_ = mesh%dyf_(ny) + mesh%dyf_(ny+1)
						this%dhb_ = mesh%dyf_(ny+1)
   					CASE DEFAULT
				END SELECT
				
			CASE('sz')
			
				SELECT CASE (this%cartDir_)
   					CASE (-3) !back boundary
						this%dhi_ = mesh%dzf_(0) + mesh%dzf_(1)
						this%dhb_ = mesh%dzf_(0)
   					CASE (3) !front boundary
						this%dhi_ = mesh%dzf_(nz) + mesh%dzf_(nz+1)
						this%dhb_ = mesh%dzf_(nz+1)
   					CASE DEFAULT
				END SELECT
			
		
		    CASE DEFAULT
		END SELECT
		   	        
    end subroutine
!========================================================================================!

!========================================================================================!
    subroutine initBoundary(this,f,bNumber,bType)
        type(boundaryField), intent(inout) :: this
        type(field), intent(in) :: f
        integer, intent(in) :: bNumber, bType

		this%bNumber_ = bNumber
		this%bType_ = bType
		
		!set internal-external
		if (bType == s_fixedValue) then
			this%isExternal_ = .TRUE.
		else if (bType == s_normalGradient) then 
			this%isExternal_ = .TRUE.
		else if (bType == s_contactAngleBC) then 	
			this%isExternal_ = .TRUE.
		else 
			this%isExternal_ = .FALSE.
		end if
		
		!set dimensions (mesh size boundaries)
		SELECT CASE (bNumber)
   		CASE (1) !left boundary
   			this%n1_ = f%ny_
   			this%n2_ = f%nz_
			this%cartDir_ = -1
   			this%dir1_ = 2
   			this%dir2_ = 3
   		CASE (2) !right boundary
   			this%n1_ = f%ny_
   			this%n2_ = f%nz_
			this%cartDir_ = 1
   			this%dir1_ = 2
   			this%dir2_ = 3
   		CASE (3) !bottom boundary
   			this%n1_ = f%nx_
   			this%n2_ = f%nz_
   			this%cartDir_ = -2
   			this%dir1_ = 1
   			this%dir2_ = 3
   		CASE (4) !top
   			this%n1_ = f%nx_
   			this%n2_ = f%nz_
   			this%cartDir_ = 2
   			this%dir1_ = 1
   			this%dir2_ = 3
   		CASE (5) !back boundary
   			this%n1_ = f%nx_
   			this%n2_ = f%ny_
   			this%cartDir_ = -3
   			this%dir1_ = 1
   			this%dir2_ = 2
   		CASE (6) !front boundary
   			this%n1_ = f%nx_
   			this%n2_ = f%ny_
   			this%cartDir_ = 3
   			this%dir1_ = 1
   			this%dir2_ = 2
   		CASE DEFAULT
		END SELECT
		
        
    end subroutine
!========================================================================================!

!========================================================================================!
    subroutine initField(this,f,bValue)
        type(boundaryField), intent(inout) :: this
        type(field), intent(in) :: f
        real(DP), intent(in) :: bValue
        type(mpiControl), pointer :: ptrMPIC
        logical, dimension(3) :: wrap


		!init field
		SELECT CASE (this%bType_)
   		CASE (s_fixedValue,s_normalGradient,s_calculatedBC,s_contactAngleBC) 
			call initBoundaryValues(this,bValue)
		CASE (s_periodicBC) 
			ptrMPIC => f%ptrMesh_%ptrMPIC_
			wrap = ptrMPIC%wrapAround_
			call initPeriodicBoundary(this,wrap)
   		CASE DEFAULT
		END SELECT
		
        
    end subroutine
!========================================================================================!

!========================================================================================!
    subroutine initBoundaryValues(this,bValue)
        type(boundaryField), intent(inout) :: this
        real(DP), intent(in) :: bValue
        	
		!init boundary value
		this%bf_ = bValue
				        
    end subroutine
!========================================================================================!

!========================================================================================!
    subroutine initPeriodicBoundary(this,wrap)
        type(boundaryField), intent(inout) :: this
        logical, dimension(3) :: wrap
        real(DP) :: x
        character(len=2) :: str
        
        !default value (set to zero)
		this%bf_ = 0.d0
		
		!check consistency with decomposeDict
		if (wrap(abs(this%cartDir_)) .eqv. .FALSE.) then
			write(str,'(I2)') this%bNumber_
			call mpiABORT('Boundary Patch ' // trim(str) // ' initialised as periodic '&
			//'but MPI wrapAround is not ')
		end if
			        
    end subroutine
!========================================================================================!

!========================================================================================!
    subroutine decomposeBoundaryPatch(this,p,ptrMPIC,patchDir,patchCoordIndex)
        type(boundaryField), intent(in) :: this
        type(boundaryField), intent(inout) :: p
        type(mpiControl), pointer, intent(in) :: ptrMPIC
        integer, intent(in) :: patchDir, patchCoordIndex
        integer, dimension(3) :: pCoord, np
        integer :: pCoordN, n
        integer :: n1, n2
        

		if (IS_MASTER) then
			!check proc-by-proc
			do n=1,ptrMPIC%nProcs_-1
			
				pCoordN = ptrMPIC%gCoords_(abs(patchDir),n)
				
				if ( pCoordN == patchCoordIndex ) then
					call sendFieldFromMaster(this,n,ptrMPIC)
				else
					call sendParallelPatchFromMaster(this,n,ptrMPIC)
				end if
			
			end do
			
		else
				!receive
				pCoord = ptrMPIC%procCoord_
				pCoordN = pCoord(abs(patchDir))
				
				if ( pCoordN == patchCoordIndex ) then
					p%cartDir_ = patchDir
					call receiveFieldFromMaster(p,ptrMPIC)
				else
					!init proc patch
					p%cartDir_ = patchDir
					call receiveParallelPatchFromMaster(p,ptrMPIC)
					
				end if
		end if
		
		
        !decompose master patch
        if (IS_MASTER) then
        
        	pCoordN = ptrMPIC%gCoords_(abs(patchDir),0)
        	
        	if ( pCoordN == patchCoordIndex ) then

        		p%bType_ = this%bType_
        		p%bNumber_ = this%bNumber_
        		p%isExternal_ = this%isExternal_
        		p%cartDir_ = this%cartDir_
        		p%dir1_ = this%dir1_
				p%dir2_ = this%dir2_
				p%dhi_ = this%dhi_
				p%dhb_ = this%dhb_
        	
        		np = ptrMPIC%nProcsAxis_
        
        		n1 = this%n1_ / np(this%dir1_)
        		n2 = this%n2_ / np(this%dir2_)
        	
        		p%n1_ = n1
        		p%n2_ = n2 
        	
				!set patch value
        		p%bf_ = this%bf_

        	else
        		!init Proc patch
        		p%bType_ = s_parallelBC
        		p%bNumber_ = this%bNumber_
        		p%isExternal_ = .FALSE.
        		p%cartDir_ = this%cartDir_
        		p%dir1_ = this%dir1_
				p%dir2_ = this%dir2_
				p%dhi_ = 0d0
				p%dhb_ = 0d0
        	
        		np = ptrMPIC%nProcsAxis_
        
        		n1 = this%n1_ / np(this%dir1_)
        		n2 = this%n2_ / np(this%dir2_)
        	
        		p%n1_ = n1
        		p%n2_ = n2 
        		
        		p%bf_ = 0d0
        		
        	end if
        	
        end if !end if IS_MASTER
        
        	        
    end subroutine
!========================================================================================!

!========================================================================================!
    subroutine sendFieldFromMaster(this,n,ptrMPIC)
        type(boundaryField), intent(in) :: this
        type(mpiControl), pointer, intent(in) :: ptrMPIC
        integer, intent(in) :: n
        character(len=:), allocatable :: buff
        integer, dimension(3) :: nProcsAxis
		integer :: ierror, sizeBuff
		integer :: tag = 0
		integer :: position
		integer :: n1, n2
		
		!reset position for packing
		position = 0
        
        nProcsAxis = ptrMPIC%nProcsAxis_
        
        n1 = this%n1_ / nProcsAxis(this%dir1_)
        n2 = this%n2_ / nProcsAxis(this%dir2_)
        
        !allocate buffer
        sizeBuff = 6*integer_size + 1*logical_size + 3*realDP_size
        call allocateArray(buff,sizeBuff)
        !send buffer size
        call MPI_Send(sizeBuff,1,MPI_INTEGER,n,tag,ptrMPIC%cartComm_,ierror)
        	
        !pack data
   		call MPI_PACK(this%bType_, 1, MPI_INTEGER, buff, sizeBuff, position, ptrMPIC%cartComm_, ierror)
   		call MPI_PACK(n1, 1, MPI_INTEGER, buff, sizeBuff, position, ptrMPIC%cartComm_, ierror)  
   		call MPI_PACK(n2, 1, MPI_INTEGER, buff, sizeBuff, position, ptrMPIC%cartComm_, ierror) 
   		call MPI_PACK(this%bNumber_, 1, MPI_INTEGER, buff, sizeBuff, position, ptrMPIC%cartComm_, ierror) 
   		call MPI_PACK(this%isExternal_, 1, MPI_LOGICAL, buff, sizeBuff, position, ptrMPIC%cartComm_, ierror) 
   		call MPI_PACK(this%dir1_, 1, MPI_INTEGER, buff, sizeBuff, position, ptrMPIC%cartComm_, ierror)
   		call MPI_PACK(this%dir2_, 1, MPI_INTEGER, buff, sizeBuff, position, ptrMPIC%cartComm_, ierror)
   		call MPI_PACK(this%dhi_, 1, MPI_DOUBLE_PRECISION, buff, sizeBuff, position, ptrMPIC%cartComm_, ierror)
   		call MPI_PACK(this%dhb_, 1, MPI_DOUBLE_PRECISION, buff, sizeBuff, position, ptrMPIC%cartComm_, ierror)
   		call MPI_PACK(this%bf_, 1, MPI_DOUBLE_PRECISION, buff, sizeBuff, position, ptrMPIC%cartComm_, ierror)

   		!send data
   		call MPI_Send(buff,sizeBuff,MPI_CHARACTER,n,tag,ptrMPIC%cartComm_,ierror)
   			        	        
    end subroutine
!========================================================================================!

!========================================================================================!
    subroutine receiveFieldFromMaster(p,ptrMPIC)
        type(boundaryField), intent(inout) :: p
        type(mpiControl), pointer, intent(in) :: ptrMPIC
        character(len=:), allocatable :: buff
		integer :: ierror, sizeBuff
		integer, dimension(MPI_STATUS_SIZE) :: status
		integer :: tag = 0
		integer :: position
		
		!reset position for unpacking
		position = 0

        !receive buffer size
        call MPI_Recv(sizeBuff,1,MPI_INTEGER,0,tag,ptrMPIC%cartComm_,status,ierror)
        !allocate buffer
        call allocateArray(buff,sizeBuff)
   		!receive data
   		call MPI_Recv(buff,sizeBuff,MPI_CHARACTER,0,tag,ptrMPIC%cartComm_,status,ierror)
        !unpack data
   		call MPI_UNPACK(buff, sizeBuff, position, p%bType_, 1, MPI_INTEGER, ptrMPIC%cartComm_, ierror)
   		call MPI_UNPACK(buff, sizeBuff, position, p%n1_, 1, MPI_INTEGER, ptrMPIC%cartComm_, ierror)   
   		call MPI_UNPACK(buff, sizeBuff, position, p%n2_, 1, MPI_INTEGER, ptrMPIC%cartComm_, ierror) 
   		call MPI_UNPACK(buff, sizeBuff, position, p%bNumber_, 1, MPI_INTEGER, ptrMPIC%cartComm_, ierror) 
   		call MPI_UNPACK(buff, sizeBuff, position, p%isExternal_, 1, MPI_LOGICAL, ptrMPIC%cartComm_, ierror)
   		call MPI_UNPACK(buff, sizeBuff, position, p%dir1_, 1, MPI_INTEGER, ptrMPIC%cartComm_, ierror)
   		call MPI_UNPACK(buff, sizeBuff, position, p%dir2_, 1, MPI_INTEGER, ptrMPIC%cartComm_, ierror)    
   		call MPI_UNPACK(buff, sizeBuff, position, p%dhi_, 1, MPI_DOUBLE_PRECISION, ptrMPIC%cartComm_, ierror)
   		call MPI_UNPACK(buff, sizeBuff, position, p%dhb_, 1, MPI_DOUBLE_PRECISION, ptrMPIC%cartComm_, ierror) 
   		call MPI_UNPACK(buff, sizeBuff, position, p%bf_, 1, MPI_DOUBLE_PRECISION, ptrMPIC%cartComm_, ierror)
        	        
    end subroutine
!========================================================================================!


!========================================================================================!
    subroutine sendParallelPatchFromMaster(this,n,ptrMPIC)
        type(boundaryField), intent(in) :: this
        type(mpiControl), pointer, intent(in) :: ptrMPIC
        integer, intent(in) :: n
        integer, dimension(5) :: buff
        integer, dimension(3) :: nProcsAxis
		integer :: tag = 0
		integer :: ierror
		integer :: n1, n2
        
        nProcsAxis = ptrMPIC%nProcsAxis_
        
        n1 = this%n1_ / nProcsAxis(this%dir1_)
        n2 = this%n2_ / nProcsAxis(this%dir2_)
        
		buff(1) = n1
		buff(2) = n2
		buff(3) = this%bNumber_
		buff(4) = this%dir1_
		buff(5) = this%dir2_

   		!send data
   		call MPI_Send(buff,size(buff),MPI_INTEGER,n,tag,ptrMPIC%cartComm_,ierror)
   			
        	        
    end subroutine
!========================================================================================!

!========================================================================================!
    subroutine receiveParallelPatchFromMaster(p,ptrMPIC)
        type(boundaryField), intent(inout) :: p
        type(mpiControl), pointer, intent(in) :: ptrMPIC
		integer, dimension(MPI_STATUS_SIZE) :: status
		integer, dimension(5) :: buff
		integer :: tag = 0
		integer :: ierror
		
		p%bType_ = s_parallelBC
		p%isExternal_ = .FALSE.
		p%bf_ = 0d0

        !receive data
        call MPI_Recv(buff,size(buff),MPI_INTEGER,0,tag,ptrMPIC%cartComm_,status,ierror)
        
        p%n1_ = 	 buff(1)
        p%n2_ = 	 buff(2)
        p%bNumber_ = buff(3)
        p%dir1_ = 	 buff(4)
        p%dir2_ = 	 buff(5)
        
        p%dhi_ = 0d0
        p%dhb_ = 0d0
        
        	        
    end subroutine
!========================================================================================!

!========================================================================================!
    subroutine updateBoundaryField(this,f)
        type(boundaryField), intent(in) :: this
        type(field), intent(inout) :: f
		
		SELECT CASE (this%bType_)
   		CASE (s_fixedValue) 
			call updateFixedValueBoundary(this,f)
   		CASE (s_normalGradient) 
			call updateNormalGradientBoundary(this,f)
		CASE (s_calculatedBC)
			call updateCalculatedBoundary(this,f)
   		CASE (s_contactAngleBC)
   			call updateContanctAngleBoundary(this,f)
   		CASE DEFAULT
		END SELECT		
        	        
    end subroutine
!========================================================================================!

!========================================================================================!
    subroutine updateFixedValueBoundary(this,f)
        type(boundaryField), intent(in) :: this
        type(field), intent(inout) :: f
		real(DP) :: dhb, dhi
        integer :: is, js, ks, ie, je, ke
        
		is = f%is_
		ie = f%ie_
		js = f%js_
		je = f%je_
		ks = f%ks_
		ke = f%ke_
		
		dhb = this%dhb_
		dhi = this%dhi_
		
		
		SELECT CASE (f%tp_)
			
			CASE ('cl')
		
				SELECT CASE (this%cartDir_)
					CASE (-1) !left boundary
						f%f_(is-1,js:je,ks:ke) = (f%f_(is,js:je,ks:ke)*dhb-this%bf_*dhi)/(dhb-dhi)
					CASE (1) !right boundary
						f%f_(ie+1,js:je,ks:ke) = (f%f_(ie,js:je,ks:ke)*dhb-this%bf_*dhi)/(dhb-dhi)
					CASE (-2) !bottom boundary
						f%f_(is:ie,js-1,ks:ke) = (f%f_(is:ie,js,ks:ke)*dhb-this%bf_*dhi)/(dhb-dhi)
					CASE (2) !top boundary
						f%f_(is:ie,je+1,ks:ke) = (f%f_(is:ie,je,ks:ke)*dhb-this%bf_*dhi)/(dhb-dhi)
					CASE (-3) !back boundary
						f%f_(is:ie,js:je,ks-1) = (f%f_(is:ie,js:je,ks)*dhb-this%bf_*dhi)/(dhb-dhi)
					CASE (3) !front boundary
						f%f_(is:ie,js:je,ke+1) = (f%f_(is:ie,js:je,ke)*dhb-this%bf_*dhi)/(dhb-dhi)		
   					CASE DEFAULT
				END SELECT
				
			CASE('sx')
			
				SELECT CASE (this%cartDir_)
					CASE (-1) !left boundary
						f%f_(is-1,js:je,ks:ke) = (f%f_(is+1,js:je,ks:ke)*dhb-this%bf_*dhi)/(dhb-dhi)
					CASE (1) !right boundary
						f%f_(ie+1,js:je,ks:ke) = (f%f_(ie-1,js:je,ks:ke)*dhb-this%bf_*dhi)/(dhb-dhi)
					CASE (-2) !bottom boundary
						f%f_(is:ie,js-1,ks:ke) = (f%f_(is:ie,js,ks:ke)*dhb-this%bf_*dhi)/(dhb-dhi)
					CASE (2) !top boundary
						f%f_(is:ie,je+1,ks:ke) = (f%f_(is:ie,je,ks:ke)*dhb-this%bf_*dhi)/(dhb-dhi)
					CASE (-3) !back boundary
						f%f_(is:ie,js:je,ks-1) = (f%f_(is:ie,js:je,ks)*dhb-this%bf_*dhi)/(dhb-dhi)
					CASE (3) !front boundary
						f%f_(is:ie,js:je,ke+1) = (f%f_(is:ie,js:je,ke)*dhb-this%bf_*dhi)/(dhb-dhi)		
   					CASE DEFAULT
				END SELECT
			
			CASE('sy')
			
				SELECT CASE (this%cartDir_)
					CASE (-1) !left boundary
						f%f_(is-1,js:je,ks:ke) = (f%f_(is,js:je,ks:ke)*dhb-this%bf_*dhi)/(dhb-dhi)
					CASE (1) !right boundary
						f%f_(ie+1,js:je,ks:ke) = (f%f_(ie,js:je,ks:ke)*dhb-this%bf_*dhi)/(dhb-dhi)
					CASE (-2) !bottom boundary
						f%f_(is:ie,js-1,ks:ke) = (f%f_(is:ie,js+1,ks:ke)*dhb-this%bf_*dhi)/(dhb-dhi)
					CASE (2) !top boundary
						f%f_(is:ie,je+1,ks:ke) = (f%f_(is:ie,je-1,ks:ke)*dhb-this%bf_*dhi)/(dhb-dhi)
					CASE (-3) !back boundary
						f%f_(is:ie,js:je,ks-1) = (f%f_(is:ie,js:je,ks)*dhb-this%bf_*dhi)/(dhb-dhi)
					CASE (3) !front boundary
						f%f_(is:ie,js:je,ke+1) = (f%f_(is:ie,js:je,ke)*dhb-this%bf_*dhi)/(dhb-dhi)		
   					CASE DEFAULT
				END SELECT
			
			CASE('sz')
			
				SELECT CASE (this%cartDir_)
					CASE (-1) !left boundary
						f%f_(is-1,js:je,ks:ke) = (f%f_(is,js:je,ks:ke)*dhb-this%bf_*dhi)/(dhb-dhi)
					CASE (1) !right boundary
						f%f_(ie+1,js:je,ks:ke) = (f%f_(ie,js:je,ks:ke)*dhb-this%bf_*dhi)/(dhb-dhi)
					CASE (-2) !bottom boundary
						f%f_(is:ie,js-1,ks:ke) = (f%f_(is:ie,js,ks:ke)*dhb-this%bf_*dhi)/(dhb-dhi)
					CASE (2) !top boundary
						f%f_(is:ie,je+1,ks:ke) = (f%f_(is:ie,je,ks:ke)*dhb-this%bf_*dhi)/(dhb-dhi)
					CASE (-3) !back boundary
						f%f_(is:ie,js:je,ks-1) = (f%f_(is:ie,js:je,ks+1)*dhb-this%bf_*dhi)/(dhb-dhi)
					CASE (3) !front boundary
						f%f_(is:ie,js:je,ke+1) = (f%f_(is:ie,js:je,ke-1)*dhb-this%bf_*dhi)/(dhb-dhi)		
   					CASE DEFAULT
				END SELECT
			
			CASE DEFAULT
		
		END SELECT
		
   			
        	        
    end subroutine
!========================================================================================!

!========================================================================================!
    subroutine updateNormalGradientBoundary(this,f)
        type(boundaryField), intent(in) :: this
        type(field), intent(inout) :: f
		real(DP) :: dhi
        integer :: is, js, ks, ie, je, ke
        
		is = f%is_
		ie = f%ie_
		js = f%js_
		je = f%je_
		ks = f%ks_
		ke = f%ke_
		
		dhi = this%dhi_
		
		SELECT CASE (f%tp_)
			
			CASE ('cl')
		
				SELECT CASE (this%cartDir_)
					CASE (-1) !left boundary
						f%f_(is-1,js:je,ks:ke) = f%f_(is,js:je,ks:ke) + this%bf_*dhi		
					CASE (1) !right boundary
						f%f_(ie+1,js:je,ks:ke) = f%f_(ie,js:je,ks:ke) + this%bf_*dhi	
					CASE (-2) !bottom boundary
						f%f_(is:ie,js-1,ks:ke) = f%f_(is:ie,js,ks:ke) + this%bf_*dhi	
					CASE (2) !top boundary
						f%f_(is:ie,je+1,ks:ke) = f%f_(is:ie,je,ks:ke) + this%bf_*dhi	
					CASE (-3) !back boundary
						f%f_(is:ie,js:je,ks-1) = f%f_(is:ie,js:je,ks) + this%bf_*dhi		
					CASE (3) !front boundary
						f%f_(is:ie,js:je,ke+1) = f%f_(is:ie,js:je,ke) + this%bf_*dhi	
   					CASE DEFAULT
				END SELECT
		
			CASE ('sx')
			
				SELECT CASE (this%cartDir_)
					CASE (-1) !left boundary
						f%f_(is-1,js:je,ks:ke) = f%f_(is+1,js:je,ks:ke) + this%bf_*dhi		
					CASE (1) !right boundary
						f%f_(ie+1,js:je,ks:ke) = f%f_(ie-1,js:je,ks:ke) + this%bf_*dhi
					CASE (-2) !bottom boundary
						f%f_(is:ie,js-1,ks:ke) = f%f_(is:ie,js,ks:ke) + this%bf_*dhi	
					CASE (2) !top boundary
						f%f_(is:ie,je+1,ks:ke) = f%f_(is:ie,je,ks:ke) + this%bf_*dhi	
					CASE (-3) !back boundary
						f%f_(is:ie,js:je,ks-1) = f%f_(is:ie,js:je,ks) + this%bf_*dhi		
					CASE (3) !front boundary
						f%f_(is:ie,js:je,ke+1) = f%f_(is:ie,js:je,ke) + this%bf_*dhi	
   					CASE DEFAULT
				END SELECT
			
			CASE ('sy')
			
				SELECT CASE (this%cartDir_)
					CASE (-1) !left boundary
						f%f_(is-1,js:je,ks:ke) = f%f_(is,js:je,ks:ke) + this%bf_*dhi		
					CASE (1) !right boundary
						f%f_(ie+1,js:je,ks:ke) = f%f_(ie,js:je,ks:ke) + this%bf_*dhi	
					CASE (-2) !bottom boundary
						f%f_(is:ie,js-1,ks:ke) = f%f_(is:ie,js+1,ks:ke) + this%bf_*dhi	
					CASE (2) !top boundary
						f%f_(is:ie,je+1,ks:ke) = f%f_(is:ie,je-1,ks:ke) + this%bf_*dhi	
					CASE (-3) !back boundary
						f%f_(is:ie,js:je,ks-1) = f%f_(is:ie,js:je,ks) + this%bf_*dhi		
					CASE (3) !front boundary
						f%f_(is:ie,js:je,ke+1) = f%f_(is:ie,js:je,ke) + this%bf_*dhi	
   					CASE DEFAULT
				END SELECT
			
			CASE ('sz')
			
				SELECT CASE (this%cartDir_)
					CASE (-1) !left boundary
						f%f_(is-1,js:je,ks:ke) = f%f_(is,js:je,ks:ke) + this%bf_*dhi		
					CASE (1) !right boundary
						f%f_(ie+1,js:je,ks:ke) = f%f_(ie,js:je,ks:ke) + this%bf_*dhi	
					CASE (-2) !bottom boundary
						f%f_(is:ie,js-1,ks:ke) = f%f_(is:ie,js,ks:ke) + this%bf_*dhi	
					CASE (2) !top boundary
						f%f_(is:ie,je+1,ks:ke) = f%f_(is:ie,je,ks:ke) + this%bf_*dhi	
					CASE (-3) !back boundary
						f%f_(is:ie,js:je,ks-1) = f%f_(is:ie,js:je,ks+1) + this%bf_*dhi		
					CASE (3) !front boundary
						f%f_(is:ie,js:je,ke+1) = f%f_(is:ie,js:je,ke-1) + this%bf_*dhi	
   					CASE DEFAULT
				END SELECT
			
			CASE DEFAULT
			
		END SELECT
   			
        	        
    end subroutine
!========================================================================================!

!========================================================================================!
    subroutine updateCalculatedBoundary(this,f)
        type(boundaryField), intent(in) :: this
        type(field), intent(inout) :: f
		real(DP) :: q0,q1,q2,q3,b1,b2,b3
		type(grid), pointer :: mesh
        integer :: is, js, ks, ie, je, ke
        
		is = f%is_
		ie = f%ie_
		js = f%js_
		je = f%je_
		ks = f%ks_
		ke = f%ke_
		
		mesh => f%ptrMesh_
			
		
		SELECT CASE (f%tp_)
			
			CASE ('cl')
		
				SELECT CASE (this%cartDir_)
					CASE (-1) !left boundary
						q0=mesh%xc_(is-1)
						q1=mesh%xc_(is)
						q2=mesh%xc_(is+1)
						q3=mesh%xc_(is+2)
						b1=(q0-q2)*(q0-q3)/((q1-q2)*(q1-q3))
						b2=(q0-q1)*(q0-q3)/((q2-q1)*(q2-q3))
						b3=(q0-q1)*(q0-q2)/((q3-q1)*(q3-q2))	

						f%f_(is-1,js:je,ks:ke) = b1*f%f_(is,js:je,ks:ke)+	&
												 b2*f%f_(is+1,js:je,ks:ke)+ &
												 b3*f%f_(is+2,js:je,ks:ke)
		

					CASE (1) !right boundary
						q0=mesh%xc_(ie+1)
						q1=mesh%xc_(ie)
						q2=mesh%xc_(ie-1)
						q3=mesh%xc_(ie-2)
						b1=(q0-q2)*(q0-q3)/((q1-q2)*(q1-q3))
						b2=(q0-q1)*(q0-q3)/((q2-q1)*(q2-q3))
						b3=(q0-q1)*(q0-q2)/((q3-q1)*(q3-q2))
				
						f%f_(ie+1,js:je,ks:ke) = b1*f%f_(ie,js:je,ks:ke)+   &
												 b2*f%f_(ie-1,js:je,ks:ke)+ &
												 b3*f%f_(ie-2,js:je,ks:ke)


					CASE (-2) !bottom boundary
						q0=mesh%yc_(js-1)
						q1=mesh%yc_(js)
						q2=mesh%yc_(js+1)
						q3=mesh%yc_(js+2)
						b1=(q0-q2)*(q0-q3)/((q1-q2)*(q1-q3))
						b2=(q0-q1)*(q0-q3)/((q2-q1)*(q2-q3))
						b3=(q0-q1)*(q0-q2)/((q3-q1)*(q3-q2))
				
						f%f_(is:ie,js-1,ks:ke) = b1*f%f_(is:ie,js,ks:ke)+	&
										 		 b2*f%f_(is:ie,js+1,ks:ke)+ &
										 		 b3*f%f_(is:ie,js+2,ks:ke)


					CASE (2) !top boundary
						q0=mesh%yc_(je+1)
						q1=mesh%yc_(je)
						q2=mesh%yc_(je-1)
						q3=mesh%yc_(je-2)
						b1=(q0-q2)*(q0-q3)/((q1-q2)*(q1-q3))
						b2=(q0-q1)*(q0-q3)/((q2-q1)*(q2-q3))
						b3=(q0-q1)*(q0-q2)/((q3-q1)*(q3-q2))
				
						f%f_(is:ie,je+1,ks:ke) = b1*f%f_(is:ie,je,ks:ke)+	&
										         b2*f%f_(is:ie,je-1,ks:ke)+ &
										         b3*f%f_(is:ie,je-2,ks:ke)

					CASE (-3) !back boundary
						q0=mesh%zc_(ks-1)
						q1=mesh%zc_(ks)
						q2=mesh%zc_(ks+1)
						q3=mesh%zc_(ks+2)
						b1=(q0-q2)*(q0-q3)/((q1-q2)*(q1-q3))
						b2=(q0-q1)*(q0-q3)/((q2-q1)*(q2-q3))
						b3=(q0-q1)*(q0-q2)/((q3-q1)*(q3-q2))	
				
						f%f_(is:ie,js:je,ks-1) = b1*f%f_(is:ie,js:je,ks)+	&
											     b2*f%f_(is:ie,js:je,ks+1)+ &
											     b3*f%f_(is:ie,js:je,ks+2)
	
					CASE (3) !front boundary
						q0=mesh%zc_(ke+1)
						q1=mesh%zc_(ke)
						q2=mesh%zc_(ke-1)
						q3=mesh%zc_(ke-2)
						b1=(q0-q2)*(q0-q3)/((q1-q2)*(q1-q3))
						b2=(q0-q1)*(q0-q3)/((q2-q1)*(q2-q3))
						b3=(q0-q1)*(q0-q2)/((q3-q1)*(q3-q2))

						f%f_(is:ie,js:je,ke+1) = b1*f%f_(is:ie,js:je,ke)+	&
												 b2*f%f_(is:ie,js:je,ke-1)+ &
												 b3*f%f_(is:ie,js:je,ke-2)

   					CASE DEFAULT
				END SELECT
			
			CASE DEFAULT
				call mpiABORT('Calculated BC for staggered fields not supported ')
		END SELECT
   			
        	        
    end subroutine
!========================================================================================!

!========================================================================================!
    subroutine updateContanctAngleBoundary(this,f)
        type(boundaryField), intent(in) :: this
        type(field), intent(inout) :: f
		type(grid), pointer :: mesh
        integer :: is, js, ks, ie, je, ke
        
		is = f%is_
		ie = f%ie_
		js = f%js_
		je = f%je_
		ks = f%ks_
		ke = f%ke_
		
		mesh => f%ptrMesh_ 
		
		SELECT CASE (f%tp_)
			
			CASE ('cl')
		
				SELECT CASE (this%cartDir_)
					CASE (-1) !left boundary
						f%f_(is-1,js:je,ks:ke) = 0.d0
					CASE (1) !right boundary
						f%f_(ie+1,js:je,ks:ke) = 0.d0
					CASE (-2) !bottom boundary
						f%f_(is:ie,js-1,ks:ke) = 0.d0
					CASE (2) !top boundary
						f%f_(is:ie,je+1,ks:ke) = 0.d0
					CASE (-3) !back boundary
						f%f_(is:ie,js:je,ks-1) = 0.d0	
					CASE (3) !front boundary
						f%f_(is:ie,js:je,ke+1) = 0.d0
   					CASE DEFAULT
				END SELECT
			
			CASE DEFAULT
				call mpiABORT('Contact angle BC for staggered fields not supported ')
		END SELECT
   			
        	        
    end subroutine
!========================================================================================!

!========================================================================================!
   subroutine coarsenBoundary(this,fc,mc,tp)
        type(boundaryField), intent(in) :: this		!fine field
		type(boundaryField), intent(inout) :: fc		!coarse field
		type(grid), intent(in) :: mc						!coarse mesh
		character(len=2), intent(in) :: tp					
		
		fc%bType_ = this%bType_
		fc%n1_ = this%n1_/2
		fc%n2_ = this%n2_/2
		fc%bNumber_ = this%bNumber_
		fc%isExternal_ = this%isExternal_ 
		fc%cartDir_ = this%cartDir_
		fc%dir1_ = this%dir1_
		fc%dir2_ = this%dir2_

   		fc%bf_ = this%bf_
      	
      	!set coarse metrics
      	call setMetrics(fc,mc,tp) 

    end subroutine
!========================================================================================!


subroutine writeBoundaryField(this)
	type(boundaryField), intent(IN) :: this

	write(s_IOunitNumber) this%bNumber_
	write(s_IOunitNumber) this%bType_
	write(s_IOunitNumber) this%bf_
		
end subroutine






