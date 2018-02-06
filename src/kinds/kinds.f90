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

module kindsMod
	implicit none
	
	!kinds
	integer, parameter :: SP = KIND(0.e0)
	integer, parameter :: DP = KIND(0.d0)
	
	!size
	integer, parameter :: realSP_size = SIZEOF(0.e0)
	integer, parameter :: realDP_size = SIZEOF(0.d0)
	integer, parameter :: integer_size = SIZEOF(0)
	integer, parameter :: logical_size = SIZEOF(.TRUE.)

end module kindsMod
