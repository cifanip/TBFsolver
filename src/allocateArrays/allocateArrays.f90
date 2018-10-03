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

module allocateArraysMod

	use kindsMod
	use errorHandlerMod

	implicit none
	
	interface allocateArray
#ifdef FAST_MODE		
		module PROCEDURE allocateArray1D_complex
		module PROCEDURE allocateArray3D_complex
#endif	
    	module PROCEDURE allocateArray1D_double
    	module PROCEDURE allocateArray2D_double
    	module PROCEDURE allocateArray3D_double
    	module PROCEDURE allocateArray4D_double
    	
    	module PROCEDURE allocateArray1D_integer
    	module PROCEDURE allocateArray2D_integer
    	module PROCEDURE allocateArray3D_integer
    	
    	module PROCEDURE allocateArray1D_logical
    	module PROCEDURE allocateArray2D_logical
    	module PROCEDURE allocateArray3D_logical
    	
    	module PROCEDURE allocateArray1D_char
    	
    	module PROCEDURE allocatePointer1D_integer
    	module PROCEDURE allocatePointer1D_double
    	module PROCEDURE allocatePointer3D_double
  	end interface
  	
	interface deallocateArray
#ifdef FAST_MODE
		module PROCEDURE deallocateArray1D_complex
		module PROCEDURE deallocateArray3D_complex
#endif 	
    	module PROCEDURE deallocateArray1D_double
    	module PROCEDURE deallocateArray2D_double
    	module PROCEDURE deallocateArray3D_double
    	module PROCEDURE deallocateArray4D_double
    	
    	module PROCEDURE deallocateArray1D_integer
    	module PROCEDURE deallocateArray2D_integer
    	module PROCEDURE deallocateArray3D_integer
    	
    	module PROCEDURE deallocateArray1D_logical
    	module PROCEDURE deallocateArray3D_logical
    	
    	module PROCEDURE deallocateArray1D_char
    	
    	module PROCEDURE deallocatePointer1D_integer
    	module PROCEDURE deallocatePointer1D_double
    	module PROCEDURE deallocatePointer3D_double
  	end interface
  	
	interface reAllocateArray
    	module PROCEDURE reAllocateArray1D_double
    	module PROCEDURE reAllocateArray2D_double
    	module PROCEDURE reAllocateArray3D_double
    	
    	module PROCEDURE reAllocateArray1D_integer
    	module PROCEDURE reAllocateArray2D_integer
    	module PROCEDURE reAllocateArray3D_integer
    	
    	module PROCEDURE reAllocateArray1D_logical
    	module PROCEDURE reAllocateArray3D_logical
    	
    	module PROCEDURE reAllocateArray1D_char
    	
    	module PROCEDURE reAllocatePointer1D_integer
    	module PROCEDURE reAllocatePointer1D_double
    	module PROCEDURE reAllocatePointer3D_double
  	end interface
  	
	interface assignArray
		module PROCEDURE assignArray1D_integer
		module PROCEDURE assignArray3D_integer
		module PROCEDURE assignArray3D_double
		module PROCEDURE assignArray3D_logical
  	end interface
	
contains

! 										allocate 
!========================================================================================!
#ifdef FAST_MODE
	subroutine allocateArray1D_complex(v,st,en) 
		integer, intent(IN) :: st, en
		complex(DP), allocatable, dimension(:), intent(INOUT) :: v
		integer :: err
		
		if (.not. allocated(v)) then
			allocate(v(st:en),STAT=err)
			
			if (err /= 0) then
				call mpiABORT('Allocation of v in allocateArray1D_complex failed ') 
			end if
			
		else
			call mpiABORT('Attempt to allocate allocated array ') 
		end if
			
	end subroutine
	
	subroutine allocateArray3D_complex(v,st1,en1,st2,en2,st3,en3) 
		integer, intent(IN) :: st1, en1, st2, en2, st3, en3
		complex(DP), allocatable, dimension(:,:,:), intent(INOUT) :: v
		integer :: err

		if (.not. allocated(v)) then
			allocate(v(st1:en1,st2:en2,st3:en3),STAT=err)
			
			if (err /= 0) then
				call mpiABORT('Allocation of v in allocateArray3D_double failed ')
			end if
			
		else
			call mpiABORT('Attempt to allocate allocated array ') 
		end if
	
		
	end subroutine
	
#endif
!========================================================================================!

!========================================================================================!
	subroutine allocateArray1D_double(v,st,en) 
		integer, intent(IN) :: st, en
		real(DP), allocatable, dimension(:), intent(INOUT) :: v
		integer :: err
		
		if (.not. allocated(v)) then
			allocate(v(st:en),STAT=err)
			
			if (err /= 0) then
				call mpiABORT('Allocation of v in allocateArray1D_double failed ') 
			end if
			
		else
			call mpiABORT('Attempt to allocate allocated array ') 
		end if
		
	end subroutine
!========================================================================================!

!========================================================================================!	
	subroutine allocateArray2D_double(v,st1,en1,st2,en2) 
		integer, intent(IN) :: st1, en1, st2, en2
		real(DP), allocatable, dimension(:,:), intent(INOUT) :: v
		integer :: err

		if (.not. allocated(v)) then
			allocate(v(st1:en1,st2:en2),STAT=err)
			
			if (err /= 0) then
				call mpiABORT('Allocation of v in allocateArray2D_double failed ')
			end if
			
		else
			call mpiABORT('Attempt to allocate allocated array ') 
		end if
		
		
	end subroutine
!========================================================================================!

!========================================================================================!	
	subroutine allocateArray3D_double(v,st1,en1,st2,en2,st3,en3) 
		integer, intent(IN) :: st1, en1, st2, en2, st3, en3
		real(DP), allocatable, dimension(:,:,:), intent(INOUT) :: v
		integer :: err

		if (.not. allocated(v)) then
			allocate(v(st1:en1,st2:en2,st3:en3),STAT=err)
			
			if (err /= 0) then
				call mpiABORT('Allocation of v in allocateArray3D_double failed ')
			end if
			
		else
			call mpiABORT('Attempt to allocate allocated array ') 
		end if
	
		
	end subroutine
!========================================================================================!

!========================================================================================!	
	subroutine allocateArray4D_double(v,st1,en1,st2,en2,st3,en3,st4,en4) 
		integer, intent(IN) :: st1, en1, st2, en2, st3, en3, st4, en4
		real(DP), allocatable, dimension(:,:,:,:), intent(INOUT) :: v
		integer :: err

		if (.not. allocated(v)) then
			allocate(v(st1:en1,st2:en2,st3:en3,st4:en4),STAT=err)
			
			if (err /= 0) then
				call mpiABORT('Allocation of v in allocateArray4D_double failed ')
			end if
			
		else
			call mpiABORT('Attempt to allocate allocated array ') 
		end if
		
		
	end subroutine
!========================================================================================!

!========================================================================================!
	subroutine allocateArray1D_integer(v,st,en) 
		integer, intent(IN) :: st, en
		integer, allocatable, dimension(:), intent(INOUT) :: v
		integer :: err

		if (.not. allocated(v)) then
			allocate(v(st:en),STAT=err)
			
			if (err /= 0) then
				call mpiABORT('Allocation of v in allocateArray1D_integer failed ') 
			end if
			
		else
			call mpiABORT('Attempt to allocate allocated array ') 
		end if
		
	end subroutine
!========================================================================================!

!========================================================================================!	
	subroutine allocateArray2D_integer(v,st1,en1,st2,en2) 
		integer, intent(IN) :: st1, en1, st2, en2
		integer, allocatable, dimension(:,:), intent(INOUT) :: v
		integer :: err

		if (.not. allocated(v)) then
			allocate(v(st1:en1,st2:en2),STAT=err)
			
			if (err /= 0) then
				call mpiABORT('Allocation of v in allocateArray2D_integer failed ')
			end if
			
		else
			call mpiABORT('Attempt to allocate allocated array ') 
		end if
		
		
	end subroutine
!========================================================================================!

!========================================================================================!	
	subroutine allocateArray3D_integer(v,st1,en1,st2,en2,st3,en3) 
		integer, intent(IN) :: st1, en1, st2, en2, st3, en3
		integer, allocatable, dimension(:,:,:), intent(INOUT) :: v
		integer :: err

		if (.not. allocated(v)) then
			allocate(v(st1:en1,st2:en2,st3:en3),STAT=err)
			
			if (err /= 0) then
				call mpiABORT('Allocation of v in allocateArray3D_integer failed ')
			end if
			
		else
			call mpiABORT('Attempt to allocate allocated array ') 
		end if
		
		
	end subroutine
!========================================================================================!

!========================================================================================!	
	subroutine allocateArray1D_logical(v,st,en) 
		integer, intent(IN) :: st, en
		logical, allocatable, dimension(:), intent(INOUT) :: v
		integer :: err

		if (.not. allocated(v)) then
			allocate(v(st:en),STAT=err)
			
			if (err /= 0) then
				call mpiABORT('Allocation of v in allocateArray1D_logical failed ')
			end if
			
		else
			call mpiABORT('Attempt to allocate allocated array ') 
		end if
		

		
	end subroutine
!========================================================================================!

!========================================================================================!	
	subroutine allocateArray2D_logical(v,st1,en1,st2,en2) 
		integer, intent(IN) :: st1, en1, st2, en2
		logical, allocatable, dimension(:,:), intent(INOUT) :: v
		integer :: err

		if (.not. allocated(v)) then
			allocate(v(st1:en1,st2:en2),STAT=err)
			
			if (err /= 0) then
				call mpiABORT('Allocation of v in allocateArray2D_logical failed ')
			end if
			
		else
			call mpiABORT('Attempt to allocate allocated array ') 
		end if
		
		
	end subroutine
!========================================================================================!

!========================================================================================!	
	subroutine allocateArray3D_logical(v,st1,en1,st2,en2,st3,en3) 
		integer, intent(IN) :: st1, en1, st2, en2, st3, en3
		logical, allocatable, dimension(:,:,:), intent(INOUT) :: v
		integer :: err

		if (.not. allocated(v)) then
			allocate(v(st1:en1,st2:en2,st3:en3),STAT=err)
			
			if (err /= 0) then
				call mpiABORT('Allocation of v in allocateArray3D_logical failed ')
			end if
			
		else
			call mpiABORT('Attempt to allocate allocated array ') 
		end if
		
		
	end subroutine
!========================================================================================!

!========================================================================================!
	subroutine allocateArray1D_char(v,size) 
		integer, intent(IN) :: size
		character(len=:), allocatable, intent(INOUT) :: v
		integer :: err

		if (.not. allocated(v)) then
			allocate(character(len=size) :: v, STAT=err)
			
			if (err /= 0) then
				call mpiABORT('Allocation of v in allocateArray1D_char failed ') 
			end if
			
		else
			call mpiABORT('Attempt to allocate allocated array ') 
		end if
		
		
	end subroutine
!========================================================================================!

!========================================================================================!
	subroutine allocatePointer1D_integer(v,st,en) 
		integer, intent(IN) :: st, en
		integer, pointer, dimension(:), intent(INOUT) :: v
		integer :: err
		
		if (.not. associated(v)) then
			allocate(v(st:en),STAT=err)
			
			if (err /= 0) then
				call mpiABORT('Allocation of v in allocatePointer1D_integer failed ') 
			end if
			
		else
			call mpiABORT('Attempt to allocate allocated pointer ') 
		end if
		
		
	end subroutine
!========================================================================================!

!========================================================================================!
	subroutine allocatePointer1D_double(v,st,en) 
		integer, intent(IN) :: st, en
		real(DP), pointer, dimension(:), intent(INOUT) :: v
		integer :: err
		
		if (.not. associated(v)) then
			allocate(v(st:en),STAT=err)
			
			if (err /= 0) then
				call mpiABORT('Allocation of v in allocatePointer1D_double failed ') 
			end if
			
		else
			call mpiABORT('Attempt to allocate allocated pointer ') 
		end if
		
		
	end subroutine
!========================================================================================!

!========================================================================================!
	subroutine allocatePointer3D_double(v,st1,en1,st2,en2,st3,en3) 
		integer, intent(IN) :: st1,en1,st2,en2,st3,en3
		real(DP), pointer, dimension(:,:,:), intent(INOUT) :: v
		integer :: err
		
		if (.not. associated(v)) then
			allocate(v(st1:en1,st2:en2,st3:en3),STAT=err)
			
			if (err /= 0) then
				call mpiABORT('Allocation of v in allocatePointer3D_double failed ') 
			end if
			
		else
			call mpiABORT('Attempt to allocate allocated pointer ') 
		end if
		
		
	end subroutine
!========================================================================================!

! 										deallocate 
!========================================================================================!
#ifdef FAST_MODE
	subroutine deallocateArray1D_complex(v) 
		complex(DP), allocatable, dimension(:), intent(INOUT) :: v

        if (allocated(v)) then
            deallocate(v)
        end if
		
	end subroutine
	
	subroutine deallocateArray3D_complex(v) 
		complex(DP), allocatable, dimension(:,:,:), intent(INOUT) :: v

        if (allocated(v)) then
            deallocate(v)
        end if
		
	end subroutine
#endif
!========================================================================================!

!========================================================================================!
	subroutine deallocateArray1D_double(v) 
		real(DP), allocatable, dimension(:), intent(INOUT) :: v

        if (allocated(v)) then
            deallocate(v)
        end if
		
	end subroutine
!========================================================================================!

!========================================================================================!	
	subroutine deallocateArray2D_double(v) 
		real(DP), allocatable, dimension(:,:), intent(INOUT) :: v

        if (allocated(v)) then
            deallocate(v)
        end if
		
	end subroutine
!========================================================================================!

!========================================================================================!	
	subroutine deallocateArray3D_double(v) 
		real(DP), allocatable, dimension(:,:,:), intent(INOUT) :: v

        if (allocated(v)) then
            deallocate(v)
        end if
		
	end subroutine
!========================================================================================!

!========================================================================================!	
	subroutine deallocateArray4D_double(v) 
		real(DP), allocatable, dimension(:,:,:,:), intent(INOUT) :: v

        if (allocated(v)) then
            deallocate(v)
        end if
		
	end subroutine
!========================================================================================!

!========================================================================================!
	subroutine deallocateArray1D_integer(v) 
		integer, allocatable, dimension(:), intent(INOUT) :: v

        if (allocated(v)) then
            deallocate(v)
        end if
		
	end subroutine
!========================================================================================!

!========================================================================================!	
	subroutine deallocateArray2D_integer(v) 
		integer, allocatable, dimension(:,:), intent(INOUT) :: v

        if (allocated(v)) then
            deallocate(v)
        end if
		
	end subroutine
!========================================================================================!

!========================================================================================!	
	subroutine deallocateArray3D_integer(v) 
		integer, allocatable, dimension(:,:,:), intent(INOUT) :: v

        if (allocated(v)) then
            deallocate(v)
        end if
		
	end subroutine
!========================================================================================!

!========================================================================================!	
	subroutine deallocateArray1D_logical(v) 
		logical, allocatable, dimension(:), intent(INOUT) :: v

        if (allocated(v)) then
            deallocate(v)
        end if
		
	end subroutine
!========================================================================================!

!========================================================================================!	
	subroutine deallocateArray3D_logical(v) 
		logical, allocatable, dimension(:,:,:), intent(INOUT) :: v

        if (allocated(v)) then
            deallocate(v)
        end if
		
	end subroutine
!========================================================================================!

!========================================================================================!
	subroutine deallocateArray1D_char(v) 
		character(len=:), allocatable, intent(OUT) :: v

        if (allocated(v)) then
            deallocate(v)
        end if
		
	end subroutine
!========================================================================================!

!========================================================================================!
	subroutine deallocatePointer1D_integer(v) 
		integer, pointer, dimension(:), intent(INOUT) :: v

        if (associated(v)) then
            deallocate(v)
            v => NULL()
		end if
		
	end subroutine
!========================================================================================!

!========================================================================================!
	subroutine deallocatePointer1D_double(v) 
		real(DP), pointer, dimension(:), intent(INOUT) :: v

        if (associated(v)) then
            deallocate(v)
            v => NULL()
		end if
		
	end subroutine
!========================================================================================!

!========================================================================================!
	subroutine deallocatePointer3D_double(v) 
		real(DP), pointer, dimension(:,:,:), intent(INOUT) :: v

        if (associated(v)) then
            deallocate(v)
            v => NULL()
		end if
		
	end subroutine
!========================================================================================!

! 										reAllocate 
!========================================================================================!
	subroutine reAllocateArray1D_double(v,st,en) 
		integer, intent(IN) :: st, en
		real(DP), allocatable, dimension(:), intent(INOUT) :: v
		
		if (allocated(v)) then
			call deallocateArray(v)
			call allocateArray(v,st,en)		
		else
			call allocateArray(v,st,en)
		end if
		
	end subroutine
!========================================================================================!

!========================================================================================!
	subroutine reAllocateArray2D_double(v,st1,en1,st2,en2)
		integer, intent(IN) :: st1, en1, st2, en2
		real(DP), allocatable, dimension(:,:), intent(INOUT) :: v
		
		if (allocated(v)) then
			call deallocateArray(v)
			call allocateArray(v,st1,en1,st2,en2)	
		else
			call allocateArray(v,st1,en1,st2,en2)
		end if

		
	end subroutine
!========================================================================================!

!========================================================================================!
	subroutine reAllocateArray3D_double(v,st1,en1,st2,en2,st3,en3)
		integer, intent(IN) :: st1, en1, st2, en2, st3, en3
		real(DP), allocatable, dimension(:,:,:), intent(INOUT) :: v
		
		if (allocated(v)) then
			call deallocateArray(v)
			call allocateArray(v,st1,en1,st2,en2,st3,en3)	
		else
			call allocateArray(v,st1,en1,st2,en2,st3,en3)
		end if

		
	end subroutine
!========================================================================================!

!========================================================================================!
	subroutine reAllocateArray1D_integer(v,st,en) 
		integer, intent(IN) :: st, en
		integer, allocatable, dimension(:), intent(INOUT) :: v
		
		if (allocated(v)) then
			call deallocateArray(v)
			call allocateArray(v,st,en)		
		else
			call allocateArray(v,st,en)
		end if

		
	end subroutine
!========================================================================================!

!========================================================================================!
	subroutine reAllocateArray2D_integer(v,st1,en1,st2,en2)
		integer, intent(IN) :: st1, en1, st2, en2
		integer, allocatable, dimension(:,:), intent(INOUT) :: v
		
		if (allocated(v)) then
			call deallocateArray(v)
			call allocateArray(v,st1,en1,st2,en2)	
		else
			call allocateArray(v,st1,en1,st2,en2)
		end if

		
	end subroutine
!========================================================================================!

!========================================================================================!
	subroutine reAllocateArray3D_integer(v,st1,en1,st2,en2,st3,en3)
		integer, intent(IN) :: st1, en1, st2, en2, st3, en3
		integer, allocatable, dimension(:,:,:), intent(INOUT) :: v
		
		if (allocated(v)) then
			call deallocateArray(v)
			call allocateArray(v,st1,en1,st2,en2,st3,en3)	
		else
			call allocateArray(v,st1,en1,st2,en2,st3,en3)
		end if

		
	end subroutine
!========================================================================================!

!========================================================================================!
	subroutine reAllocateArray1D_logical(v,st1,en1)
		integer, intent(IN) :: st1, en1
		logical, allocatable, dimension(:), intent(INOUT) :: v
		
		if (allocated(v)) then
			call deallocateArray(v)
			call allocateArray(v,st1,en1)	
		else
			call allocateArray(v,st1,en1)
		end if

		
	end subroutine
!========================================================================================!

!========================================================================================!
	subroutine reAllocateArray3D_logical(v,st1,en1,st2,en2,st3,en3)
		integer, intent(IN) :: st1, en1, st2, en2, st3, en3
		logical, allocatable, dimension(:,:,:), intent(INOUT) :: v
		
		if (allocated(v)) then
			call deallocateArray(v)
			call allocateArray(v,st1,en1,st2,en2,st3,en3)	
		else
			call allocateArray(v,st1,en1,st2,en2,st3,en3)
		end if

		
	end subroutine
!========================================================================================!

!========================================================================================!
	subroutine reAllocateArray1D_char(v,size) 
		integer, intent(IN) :: size
		character(len=:), allocatable, intent(INOUT) :: v
		
		if (allocated(v)) then
			call deallocateArray(v)
			call allocateArray(v,size)		
		else
			call allocateArray(v,size)
		end if

		
	end subroutine
!========================================================================================!

!========================================================================================!
	subroutine reAllocatePointer1D_integer(v,st,en) 
		integer, intent(IN) :: st, en
		integer, pointer, dimension(:), intent(INOUT) :: v
		
		if (associated(v)) then
			call deallocateArray(v)
			call allocateArray(v,st,en)		
		else
			call allocateArray(v,st,en)
		end if

		
	end subroutine
!========================================================================================!

!========================================================================================!
	subroutine reAllocatePointer1D_double(v,st,en) 
		integer, intent(IN) :: st, en
		real(DP), pointer, dimension(:), intent(INOUT) :: v
		
		if (associated(v)) then
			call deallocateArray(v)
			call allocateArray(v,st,en)		
		else
			call allocateArray(v,st,en)
		end if

		
	end subroutine
!========================================================================================!

!========================================================================================!
	subroutine reAllocatePointer3D_double(v,st1,en1,st2,en2,st3,en3) 
		integer, intent(IN) :: st1,en1,st2,en2,st3,en3
		real(DP), pointer, dimension(:,:,:), intent(INOUT) :: v
		
		if (associated(v)) then
			call deallocateArray(v)
			call allocateArray(v,st1,en1,st2,en2,st3,en3)		
		else
			call allocateArray(v,st1,en1,st2,en2,st3,en3)
		end if

		
	end subroutine
!========================================================================================!

! 										assign 
!========================================================================================!
	subroutine assignArray3D_double(lhs,rhs) 
		real(DP), allocatable, dimension(:,:,:), intent(INOUT) :: lhs
		real(DP), allocatable, dimension(:,:,:), intent(IN) :: rhs
		integer :: st1,en1,st2,en2,st3,en3
		
		if (allocated(rhs)) then
			
			st1=lbound(rhs,1)
			st2=lbound(rhs,2)
			st3=lbound(rhs,3)
			en1=ubound(rhs,1)
			en2=ubound(rhs,2)
			en3=ubound(rhs,3)
			
			call reAllocateArray(lhs,st1,en1,st2,en2,st3,en3)
			lhs=rhs
			
		else
			call mpiABORT('Attempt to assign to a deallocated array ') 
		end if
		
	end subroutine
!========================================================================================!

!========================================================================================!
	subroutine assignArray1D_integer(lhs,rhs) 
		integer, allocatable, dimension(:), intent(INOUT) :: lhs
		integer, allocatable, dimension(:), intent(IN) :: rhs
		integer :: st1,en1
		
		if (allocated(rhs)) then
			
			st1=lbound(rhs,1)
			en1=ubound(rhs,1)
			
			call reAllocateArray(lhs,st1,en1)
			lhs=rhs
			
		else
			call mpiABORT('Attempt to assign to a deallocated array ') 
		end if
		
	end subroutine
!========================================================================================!

!========================================================================================!
	subroutine assignArray3D_integer(lhs,rhs) 
		integer, allocatable, dimension(:,:,:), intent(INOUT) :: lhs
		integer, allocatable, dimension(:,:,:), intent(IN) :: rhs
		integer :: st1,en1,st2,en2,st3,en3
		
		if (allocated(rhs)) then
			
			st1=lbound(rhs,1)
			st2=lbound(rhs,2)
			st3=lbound(rhs,3)
			en1=ubound(rhs,1)
			en2=ubound(rhs,2)
			en3=ubound(rhs,3)
			
			call reAllocateArray(lhs,st1,en1,st2,en2,st3,en3)
			lhs=rhs
			
		else
			call mpiABORT('Attempt to assign to a deallocated array ') 
		end if
		
	end subroutine
!========================================================================================!

!========================================================================================!
	subroutine assignArray3D_logical(lhs,rhs) 
		logical, allocatable, dimension(:,:,:), intent(INOUT) :: lhs
		logical, allocatable, dimension(:,:,:), intent(IN) :: rhs
		integer :: st1,en1,st2,en2,st3,en3
		
		if (allocated(rhs)) then
			
			st1=lbound(rhs,1)
			st2=lbound(rhs,2)
			st3=lbound(rhs,3)
			en1=ubound(rhs,1)
			en2=ubound(rhs,2)
			en3=ubound(rhs,3)
			
			call reAllocateArray(lhs,st1,en1,st2,en2,st3,en3)
			lhs=rhs
			
		else
			call mpiABORT('Attempt to assign to a deallocated array ') 
		end if
		
	end subroutine
!========================================================================================!
	
end module allocateArraysMod


	



