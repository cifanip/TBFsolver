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

	! scalarBoundaryField
	type, public :: scalarBoundaryField
		
		integer :: bType_
		integer :: n1_, n2_
		integer :: bNumber_
		logical :: isExternal_
		integer :: cartDir_ !x: +-1; y: +-2; z: +-3
		integer :: dir1_, dir2_ !cart directions along n1 and n2
		
		real(DP) :: bf_
		
		real(DP) :: dhi_
		real(DP) :: dhb_
		
	end type
	
	private :: initDeafaultBoundaryField
	private :: readBoundaryField
	private :: setMetrics
	private :: initBoundary
	private :: initField
	private :: initBoundaryValues
	private :: initPeriodicBoundary
	private :: decomposeBoundaryPatch
	private :: sendFieldFromMaster
	private :: receiveFieldFromMaster
	private :: sendParallelPatchFromMaster
	private :: receiveParallelPatchFromMaster
	private :: updateBoundaryField
	private :: updateFixedValueBoundary
	private :: updateNormalGradientBoundary
	private :: updateCalculatedBoundary
	private :: updateContanctAngleBoundary
	private :: coarsenBoundary
	private :: writeBoundaryField
	
	public :: scalarBoundaryFieldCTOR
  