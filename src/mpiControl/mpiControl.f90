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

module mpiControlMod
	
	use allocateArraysMod
	use parFileMod
	
	implicit none
	

	type, public :: mpiControl
	
		!cartesian communicator
		integer :: cartComm_
		
		integer, dimension(3) :: nProcsAxis_
		integer, dimension(3) :: procCoord_
		logical, dimension(3) :: wrapAround_
		
		integer, allocatable, dimension(:,:) :: gCoords_ 
		
		integer :: rank_
		integer :: nProcs_
		
		!mesh size
		integer, private :: nx_, ny_, nz_
		
		!neighbours procs
		integer :: leftNe_, rightNe_, backNe_, frontNe_, bottomNe_, topNe_
		
		!halo dim
		integer, private :: haloDim_
		
		contains	
			 

	end type
	
	private :: readMeshpfile
	private :: readDecomposepfile
	private :: setNeighboursProc
	private :: setGlobalCoord
#ifdef FAST_MODE
	private :: checkPeriodicDirFFT
	private :: checkPencilDecomp
#endif
	
	public :: mpiControlCTOR
  	
	
contains


!========================================================================================!
	subroutine mpiControlCTOR(this) 
		type(mpiControl), intent(out) :: this
		integer :: ierr
		
        !set number of procs
        call MPI_COMM_SIZE(MPI_COMM_WORLD, this%nProcs_, ierr)
		
		!read meshpfile
		call readMeshpfile(this)
		
		!read decomposepfile
		call readDecomposepfile(this)
		
		!create cartesian communicator
		call MPI_CART_CREATE(MPI_COMM_WORLD,   &
    					3,                	   &
    					this%nProcsAxis_,  	   &
        				this%wrapAround_, 	   &
        				1,                     &
        				this%cartComm_,		   &
        				ierr)
        
        !set proc rank
        call MPI_COMM_RANK(this%cartComm_, this%rank_, ierr)
        
		!set proc coordinates
		call MPI_CART_COORDS(this%cartComm_, this%rank_, 3, this%procCoord_, ierr)

		!set proc neighbours
		call setNeighboursProc(this)
		
		!set global cart coordinates to All procs
		call setGlobalCoord(this)

#ifdef FAST_MODE
		call checkPeriodicDirFFT(this)
		call checkPencilDecomp(this)
#endif
		
		        		
		  
	end subroutine
!========================================================================================!


!========================================================================================!
	subroutine readDecomposepfile(this)
		type(mpiControl), intent(inout) :: this
		type(parFile) :: pfile
		
		call parFileCTOR(pfile,'decompose','specs')
		
		!set number of procs per axis
		if (IS_PAR) then
			call readParameter(pfile,this%nProcsAxis_(1),'px')
			call readParameter(pfile,this%nProcsAxis_(2),'py')
			call readParameter(pfile,this%nProcsAxis_(3),'pz') 
		else
			this%nProcsAxis_(1) = 1
			this%nProcsAxis_(2) = 1
			this%nProcsAxis_(3) = 1
		end if
		
		if ( (this%nProcs_ /= this%nProcsAxis_(1)*this%nProcsAxis_(2)*this%nProcsAxis_(3)) &
		    .AND. (IS_PAR)) then
			call mpiABORT('Number of read in procs not equal to MPI_SIZE')
		end if
		
		!check even division mesh size - number of procs 
		if ( mod(this%nx_,this%nProcsAxis_(1)) /= 0) then
			call mpiABORT('Number of procs in x direction does not divide evenly nx')
		end if
		if ( mod(this%ny_,this%nProcsAxis_(2)) /= 0) then
			call mpiABORT('Number of procs in y direction does not divide evenly ny')
		end if
		if ( mod(this%nz_,this%nProcsAxis_(3)) /= 0) then
			call mpiABORT('Number of procs in z direction does not divide evenly nz')
		end if
		
		!set wrap around option 
		call readParameter(pfile,this%wrapAround_(1),'wrap_x')
		call readParameter(pfile,this%wrapAround_(2),'wrap_y')
		call readParameter(pfile,this%wrapAround_(3),'wrap_z') 

	end subroutine
!========================================================================================!

!========================================================================================!
	subroutine readMeshpfile(this)
		type(mpiControl), intent(inout) :: this
		type(parFile) :: pfile
		
		 call parFileCTOR(pfile,'mesh','specs')
		
		call readParameter(pfile,this%nx_,'nx')
		call readParameter(pfile,this%ny_,'ny')
		call readParameter(pfile,this%nz_,'nz')


	end subroutine
!========================================================================================!

!========================================================================================!
	subroutine setNeighboursProc(this)
		type(mpiControl), intent(inout) :: this
		integer :: ierror
		
		!along x
		call MPI_CART_SHIFT(this%cartComm_, 0,  1, this%leftNe_, this%rightNe_, ierror)

		!along y
		call MPI_CART_SHIFT(this%cartComm_, 1,  1, this%bottomNe_, this%topNe_, ierror)

		!along z
		call MPI_CART_SHIFT(this%cartComm_, 2,  1, this%backNe_, this%frontNe_, ierror)


	end subroutine
!========================================================================================!

!========================================================================================!
	subroutine setGlobalCoord(this)
		type(mpiControl), intent(inout) :: this
		integer :: ierror
		
		!gather global proc coords
		call allocateArray(this%gCoords_,1,3,0,this%nProcs_-1)
		
		call MPI_ALLGather(this%procCoord_,3,MPI_INTEGER,this%gCoords_,3,MPI_INTEGER,this%cartComm_,ierror)

	end subroutine
!========================================================================================!

#ifdef FAST_MODE
!========================================================================================!
	subroutine checkPeriodicDirFFT(this)
		type(mpiControl), intent(in) :: this
		
		if (.NOT.(this%wrapAround_(1))) then
			call mpiAbort('Direction x is not periodic')
		end if
		
		if (.NOT.(this%wrapAround_(3))) then
			call mpiAbort('Direction z is not periodic')
		end if

	end subroutine
!========================================================================================!

!========================================================================================!
	subroutine checkPencilDecomp(this)
		type(mpiControl), intent(in) :: this
		
		if (this%nProcsAxis_(2)>1) then
			call mpiAbort('n procs along pencil > 1')
		end if

	end subroutine
!========================================================================================!
#endif	
	
end module mpiControlMod


