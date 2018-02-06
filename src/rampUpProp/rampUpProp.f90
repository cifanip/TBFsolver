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

module rampUpPropMod

	use dictionaryMod
	
	implicit none
	
	type, public :: rampUpProp
	
		real(DP), private :: Tr_
		real(DP), private :: startValue_,endValue_
		logical, private :: isRampUp_


		contains	
	end type
	
	private :: reSet
	private :: rampUp
	
	public :: rampUpPropCTOR
	public :: updateProp

	
contains

!========================================================================================!
	subroutine rampUpPropCTOR(this,p,startValue,endValue)
		type(rampUpProp), intent(out) :: this
		real(DP), intent(inout) :: p
		real(DP) :: startValue,endValue
		type(dictionary) :: dict
		
		call dictionaryCTOR(dict,'timeControl','specs')
		call readParameter(dict,this%Tr_,'Tr')
		call readParameter(dict,this%isRampUp_,'isRampUp')
		
		this%startValue_ = startValue
		this%endValue_ = endValue
		
		call reSet(this,p)


	end subroutine
!========================================================================================!

!========================================================================================!
	subroutine reSet(this,p)
		type(rampUpProp), intent(in) :: this
		real(DP), intent(inout) :: p
		
		if (this%isRampUp_) then
			p = this%startValue_
		end if

	end subroutine
!========================================================================================!

!========================================================================================!
	function rampUp(this,t) result(isRamp)
		type(rampUpProp), intent(in) :: this
		real(DP), intent(in) :: t
		logical :: isRamp
		
		if (t <= this%Tr_) then
			isRamp = .TRUE.
		else
			isRamp = .FALSE.
		end if

	end function
!========================================================================================!

!========================================================================================!
	subroutine updateProp(this,t,p)
		type(rampUpProp), intent(in) :: this
		real(DP), intent(inout) :: p
		real(DP), intent(in) :: t
		
		if (rampUp(this,t) .AND. this%isRampUp_) then
			p=this%startValue_+( (this%endValue_-this%startValue_)/this%Tr_  )*t
		else
			p = this%endValue_
		end if

	end subroutine
!========================================================================================!


end module rampUpPropMod


	



