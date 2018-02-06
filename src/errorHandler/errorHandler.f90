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

module errorHandlerMod

	use initMpiGvarMod

	implicit none
	
	public :: mpiABORT
	
contains

!========================================================================================!
	subroutine mpiABORT(msg,opt) 
		character(len=*), intent(in) :: msg
		integer, optional :: opt
		integer :: errorCode, ierror

		if (present(opt)) then
			errorCode = opt
		else
			errorCode = -1
		end if
		
		!if (IS_MASTER) then
			write(*,*) '/***************************************************************/'
			write(*,*) 'mpiABORT CALLED with ERROR CODE: ', errorCode
			write(*,*) 'ERROR MESSAGE: ', msg
			write(*,*) '/***************************************************************/'
		!end if
		
		call MPI_ABORT(MPI_COMM_WORLD,errorcode,ierror)
		
	end subroutine
!========================================================================================!

	
end module errorHandlerMod


	



