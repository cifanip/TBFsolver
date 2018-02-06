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

module dictionaryMod

	use errorHandlerMod
	use kindsMod
	use formatsMod
	
	implicit none
	
	type, public :: dictionary
	
		character(len=:), allocatable :: fileName_
		character(len=:), allocatable :: filePath_
		
	end type
	
	!generic read parameter routines
	interface readParameter
		module PROCEDURE readParInt
		module PROCEDURE readParBool
		module PROCEDURE readParReal
	end interface
	
	private :: frstnb
	
	public :: dictionaryCTOR


contains

!========================================================================================!
	subroutine dictionaryCTOR(this,fileName,fileDir)
		type(dictionary), intent(out) :: this
		character(len=*), intent(in) :: fileName
		character(len=*), intent(in) :: fileDir 
		
		this%fileName_ = fileName
		this%filePath_ = fileDir//'/'//fileName
		
	end subroutine
!========================================================================================!

!========================================================================================!
	subroutine readParInt(this,x,name) 
		type(dictionary), intent(in) :: this
		integer, intent(out) :: x
		character(len=*), intent(in) :: name
		character(len=100) :: whole
		integer :: i1, i2, ios
		
		open(UNIT=s_IOunitNumber,FILE=this%filePath_,STATUS='OLD',ACTION='READ')
		do 
			read(s_IOunitNumber,s_charFormat,IOSTAT=ios) whole
			
				if (ios /= 0) then
					call mpiABORT('Parameter '//name//' not found in file '//this%fileName_)
				else
					!parse whole
					i1 = index(whole,' ')
						if (whole(1:i1-1) == name) then
							i2 = frstnb(whole(i1+1:),name)
							read(whole(i1+i2:),s_intFormat) x
							exit
						end if
				end if
		end do
		close(s_IOunitNumber)
	
	end subroutine
!========================================================================================!

!========================================================================================!
	subroutine readParBool(this,x,name) 
		type(dictionary), intent(in) :: this
		logical, intent(out) :: x
		character(len=*), intent(in) :: name
		character(len=100) :: whole
		integer :: i1, i2, ios
		
		open(UNIT=s_IOunitNumber,FILE=this%filePath_,STATUS='OLD',ACTION='READ')
		do 
			read(s_IOunitNumber,s_charFormat,IOSTAT=ios) whole
			
				if (ios /= 0) then
					call mpiABORT('Parameter '//name//' not found in file '//this%fileName_)
				else
					!parse whole
					i1 = index(whole,' ')
						if (whole(1:i1-1) == name) then
							i2 = frstnb(whole(i1+1:),name)
							read(whole(i1+i2:),s_logicalFormat) x
							exit
						end if
				end if
		end do
		close(s_IOunitNumber)
	
	end subroutine
!========================================================================================!

!========================================================================================!
	subroutine readParReal(this,x,name) 
		type(dictionary), intent(in) :: this
		real(DP), intent(out) :: x
		character(len=*), intent(in) :: name
		character(len=100) :: whole
		integer :: i1, i2, ios
		
		open(UNIT=s_IOunitNumber,FILE=this%filePath_,STATUS='OLD',ACTION='READ')
		do 
			read(s_IOunitNumber,s_charFormat,IOSTAT=ios) whole
			
				if (ios /= 0) then
					call mpiABORT('Parameter '//name//' not found in file '//this%fileName_)
				else
					!parse whole
					i1 = index(whole,' ')
						if (whole(1:i1-1) == name) then
							i2 = frstnb(whole(i1+1:),name)
							read(whole(i1+i2:),s_doubleFormat) x
							exit
						end if
				end if
		end do
		close(s_IOunitNumber)
	
	end subroutine
!========================================================================================!

!========================================================================================!
	function frstnb(str,par) result(ind)
		character(len=*), intent(in) :: str, par
		integer :: ind, i
		
		do i=1,len(str)
			if (str(i:i) /= ' ') then
				ind = i
				return 
			end if
		end do
		
		call mpiABORT('Invalid entry while parsing parameter: ' // par)
	
	end function
!========================================================================================!

	
end module dictionaryMod


	



