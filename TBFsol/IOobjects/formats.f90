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

module formatsMod
	
	implicit none
	
	integer, protected :: s_nFilesToWrite = 0
	
	integer, parameter :: s_IOunitNumber = 50
	character(len=4), parameter :: s_intFormat = '(I6)'
	character(len=3), parameter :: s_charFormat = '(A)'
	character(len=4), parameter :: s_logicalFormat = '(L7)'
	character(len=11), parameter :: s_doubleFormat = '(ES25.15E3)'
	character(len=10), parameter :: s_outputFormat = '(ES11.4E2)'
	
end module formatsMod
	


	



