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

module ompRoutinesMod

	use kindsMod

	implicit none
	
	interface assign_omp
    	module PROCEDURE assign_whole_omp
    	module PROCEDURE assign_section_omp
  	end interface

	interface reduceSum_omp
    	module PROCEDURE reduceSum_whole_omp
    	module PROCEDURE reduceSum_section_omp
  	end interface
  	
	interface reduceMax_omp
    	module PROCEDURE reduceMax_whole_omp
    	module PROCEDURE reduceMax_whole_int_omp
    	module PROCEDURE reduceMax_section_omp
  	end interface
  	
	interface reduceMin_omp
    	module PROCEDURE reduceMin_whole_omp
    	module PROCEDURE reduceMin_section_omp
  	end interface
	
	interface reduceSqrSum_omp
    	module PROCEDURE reduceSqrSum_whole_omp
    	module PROCEDURE reduceSqrSum_section_omp
  	end interface
  	
  	interface unarySum_omp
    	module PROCEDURE unarySum_whole_omp
    	module PROCEDURE unarySum_section_omp
  	end interface
	
contains

!========================================================================================!
	subroutine assign_whole_omp(va,vb)
		real(DP), allocatable, dimension(:,:,:), intent(inout) :: va
		real(DP), allocatable, dimension(:,:,:), intent(in) :: vb
		integer :: is,ie,js,je,ks,ke
		integer :: i,j,k
		
		is = lbound(va,1)
		ie = ubound(va,1)
		js = lbound(va,2)
		je = ubound(va,2)
		ks = lbound(va,3)
		ke = ubound(va,3)
		
		!$OMP PARALLEL DO DEFAULT(none) &
		!$OMP SHARED(va,vb) &
		!$OMP SHARED(is,ie,js,je,ks,ke) &
		!$OMP PRIVATE(i,j,k) 
		do k=ks,ke
			do j=js,je
				do i=is,ie
					va(i,j,k)=vb(i,j,k)
				end do
			end do
		end do
		!$OMP END PARALLEL DO
		
	end subroutine
!========================================================================================!

!========================================================================================!
	subroutine assign_section_omp(va,vb,is,ie,js,je,ks,ke)
		real(DP), allocatable, dimension(:,:,:), intent(inout) :: va
		real(DP), allocatable, dimension(:,:,:), intent(in) :: vb
		integer, intent(in) :: is,ie,js,je,ks,ke
		integer :: i,j,k
		
		!$OMP PARALLEL DO DEFAULT(none) &
		!$OMP SHARED(va,vb) &
		!$OMP SHARED(is,ie,js,je,ks,ke) &
		!$OMP PRIVATE(i,j,k) 
		do k=ks,ke
			do j=js,je
				do i=is,ie
					va(i,j,k)=vb(i,j,k)
				end do
			end do
		end do
		!$OMP END PARALLEL DO
		
	end subroutine
!========================================================================================!

!========================================================================================!
	subroutine set2zero_omp(v)
		real(DP), allocatable, dimension(:,:,:), intent(inout) :: v
		integer :: is,ie,js,je,ks,ke
		integer :: i,j,k
		
		is = lbound(v,1)
		ie = ubound(v,1)
		js = lbound(v,2)
		je = ubound(v,2)
		ks = lbound(v,3)
		ke = ubound(v,3)
		
		!$OMP PARALLEL DO DEFAULT(none) &
		!$OMP SHARED(v) &
		!$OMP SHARED(is,ie,js,je,ks,ke) &
		!$OMP PRIVATE(i,j,k) 
		do k=ks,ke
			do j=js,je
				do i=is,ie
					v(i,j,k)=0.d0
				end do
			end do
		end do
		!$OMP END PARALLEL DO
		
	end subroutine
!========================================================================================!

!========================================================================================!
	subroutine resetLogicalField_omp(v,bool)
		logical, allocatable, dimension(:,:,:), intent(inout) :: v
		logical, intent(in) :: bool
		integer :: is,ie,js,je,ks,ke
		integer :: i,j,k
		
		is = lbound(v,1)
		ie = ubound(v,1)
		js = lbound(v,2)
		je = ubound(v,2)
		ks = lbound(v,3)
		ke = ubound(v,3)
		
		!$OMP PARALLEL DO DEFAULT(none) &
		!$OMP SHARED(v,bool) &
		!$OMP SHARED(is,ie,js,je,ks,ke) &
		!$OMP PRIVATE(i,j,k) 
		do k=ks,ke
			do j=js,je
				do i=is,ie
					v(i,j,k)=bool
				end do
			end do
		end do
		!$OMP END PARALLEL DO
		
	end subroutine
!========================================================================================!

!========================================================================================!
	subroutine reduceSum_whole_omp(v,r)
		real(DP), allocatable, dimension(:,:,:), intent(in) :: v
		real(DP), intent(out) :: r
		integer :: is,ie,js,je,ks,ke
		integer :: i,j,k
		
		is = lbound(v,1)
		ie = ubound(v,1)
		js = lbound(v,2)
		je = ubound(v,2)
		ks = lbound(v,3)
		ke = ubound(v,3)
		
		r = 0.d0
		
		!$OMP PARALLEL DO DEFAULT(none) &
		!$OMP SHARED(v) &
		!$OMP SHARED(is,ie,js,je,ks,ke) &
		!$OMP PRIVATE(i,j,k) &
		!$OMP REDUCTION(+:r)
		do k=ks,ke
			do j=js,je
				do i=is,ie
					r = r + v(i,j,k)
				end do
			end do
		end do
		!$OMP END PARALLEL DO
		
	end subroutine
!========================================================================================!

!========================================================================================!
	subroutine reduceSum_section_omp(v,is,ie,js,je,ks,ke,r)
		real(DP), allocatable, dimension(:,:,:), intent(in) :: v
		real(DP), intent(out) :: r
		integer, intent(in) :: is,ie,js,je,ks,ke
		integer :: i,j,k
		
		r = 0.d0
		
		!$OMP PARALLEL DO DEFAULT(none) &
		!$OMP SHARED(v) &
		!$OMP SHARED(is,ie,js,je,ks,ke) &
		!$OMP PRIVATE(i,j,k) &
		!$OMP REDUCTION(+:r)
		do k=ks,ke
			do j=js,je
				do i=is,ie
					r = r + v(i,j,k)
				end do
			end do
		end do
		!$OMP END PARALLEL DO
		
	end subroutine
!========================================================================================!

!========================================================================================!
	subroutine reduceMax_whole_omp(v,r)
		real(DP), allocatable, dimension(:,:,:), intent(in) :: v
		real(DP), intent(out) :: r
		integer :: is,ie,js,je,ks,ke
		integer :: i,j,k
		
		is = lbound(v,1)
		ie = ubound(v,1)
		js = lbound(v,2)
		je = ubound(v,2)
		ks = lbound(v,3)
		ke = ubound(v,3)
		
		r = v(is,js,ks)
		
		!$OMP PARALLEL DO DEFAULT(none) &
		!$OMP SHARED(v) &
		!$OMP SHARED(is,ie,js,je,ks,ke) &
		!$OMP PRIVATE(i,j,k) &
		!$OMP REDUCTION(max:r)
		do k=ks,ke
			do j=js,je
				do i=is,ie
					r = max(v(i,j,k),r)
				end do
			end do
		end do
		!$OMP END PARALLEL DO
		
	end subroutine
!========================================================================================!

!========================================================================================!
	subroutine reduceMax_whole_int_omp(v,r)
		integer, allocatable, dimension(:,:,:), intent(in) :: v
		integer, intent(out) :: r
		integer :: is,ie,js,je,ks,ke
		integer :: i,j,k
		
		is = lbound(v,1)
		ie = ubound(v,1)
		js = lbound(v,2)
		je = ubound(v,2)
		ks = lbound(v,3)
		ke = ubound(v,3)
		
		r = v(is,js,ks)
		
		!$OMP PARALLEL DO DEFAULT(none) &
		!$OMP SHARED(v) &
		!$OMP SHARED(is,ie,js,je,ks,ke) &
		!$OMP PRIVATE(i,j,k) &
		!$OMP REDUCTION(max:r)
		do k=ks,ke
			do j=js,je
				do i=is,ie
					r = max(v(i,j,k),r)
				end do
			end do
		end do
		!$OMP END PARALLEL DO
		
	end subroutine
!========================================================================================!

!========================================================================================!
	subroutine reduceMax_section_omp(v,is,ie,js,je,ks,ke,r)
		real(DP), allocatable, dimension(:,:,:), intent(in) :: v
		real(DP), intent(out) :: r
		integer, intent(in) :: is,ie,js,je,ks,ke
		integer :: i,j,k
		
		r = v(is,js,ks)
		
		!$OMP PARALLEL DO DEFAULT(none) &
		!$OMP SHARED(v) &
		!$OMP SHARED(is,ie,js,je,ks,ke) &
		!$OMP PRIVATE(i,j,k) &
		!$OMP REDUCTION(max:r)
		do k=ks,ke
			do j=js,je
				do i=is,ie
					r = max(v(i,j,k),r)
				end do
			end do
		end do
		!$OMP END PARALLEL DO
		
	end subroutine
!========================================================================================!

!========================================================================================!
	subroutine reduceMin_whole_omp(v,r)
		real(DP), allocatable, dimension(:,:,:), intent(in) :: v
		real(DP), intent(out) :: r
		integer :: is,ie,js,je,ks,ke
		integer :: i,j,k
		
		is = lbound(v,1)
		ie = ubound(v,1)
		js = lbound(v,2)
		je = ubound(v,2)
		ks = lbound(v,3)
		ke = ubound(v,3)
		
		r = v(is,js,ks)
		
		!$OMP PARALLEL DO DEFAULT(none) &
		!$OMP SHARED(v) &
		!$OMP SHARED(is,ie,js,je,ks,ke) &
		!$OMP PRIVATE(i,j,k) &
		!$OMP REDUCTION(min:r)
		do k=ks,ke
			do j=js,je
				do i=is,ie
					r = min(v(i,j,k),r)
				end do
			end do
		end do
		!$OMP END PARALLEL DO
		
	end subroutine
!========================================================================================!

!========================================================================================!
	subroutine reduceMin_section_omp(v,is,ie,js,je,ks,ke,r)
		real(DP), allocatable, dimension(:,:,:), intent(in) :: v
		real(DP), intent(out) :: r
		integer, intent(in) :: is,ie,js,je,ks,ke
		integer :: i,j,k
		
		r = v(is,js,ks)
		
		!$OMP PARALLEL DO DEFAULT(none) &
		!$OMP SHARED(v) &
		!$OMP SHARED(is,ie,js,je,ks,ke) &
		!$OMP PRIVATE(i,j,k) &
		!$OMP REDUCTION(min:r)
		do k=ks,ke
			do j=js,je
				do i=is,ie
					r = min(v(i,j,k),r)
				end do
			end do
		end do
		!$OMP END PARALLEL DO
		
	end subroutine
!========================================================================================!

!========================================================================================!
	subroutine reduceSqrSum_whole_omp(v,r)
		real(DP), allocatable, dimension(:,:,:), intent(in) :: v
		real(DP), intent(out) :: r
		integer :: is,ie,js,je,ks,ke
		integer :: i,j,k
		
		is = lbound(v,1)
		ie = ubound(v,1)
		js = lbound(v,2)
		je = ubound(v,2)
		ks = lbound(v,3)
		ke = ubound(v,3)
		
		r = 0.d0
		
		!$OMP PARALLEL DO DEFAULT(none) &
		!$OMP SHARED(v) &
		!$OMP SHARED(is,ie,js,je,ks,ke) &
		!$OMP PRIVATE(i,j,k) &
		!$OMP REDUCTION(+:r)
		do k=ks,ke
			do j=js,je
				do i=is,ie
					r = r + v(i,j,k)*v(i,j,k)
				end do
			end do
		end do
		!$OMP END PARALLEL DO
		
	end subroutine
!========================================================================================!

!========================================================================================!
	subroutine reduceSqrSum_section_omp(v,is,ie,js,je,ks,ke,r)
		real(DP), allocatable, dimension(:,:,:), intent(in) :: v
		real(DP), intent(out) :: r
		integer, intent(in) :: is,ie,js,je,ks,ke
		integer :: i,j,k
		
		r = 0.d0
		
		!$OMP PARALLEL DO DEFAULT(none) &
		!$OMP SHARED(v) &
		!$OMP SHARED(is,ie,js,je,ks,ke) &
		!$OMP PRIVATE(i,j,k) &
		!$OMP REDUCTION(+:r)
		do k=ks,ke
			do j=js,je
				do i=is,ie
					r = r + v(i,j,k)*v(i,j,k)
				end do
			end do
		end do
		!$OMP END PARALLEL DO
		
	end subroutine
!========================================================================================!

!========================================================================================!
	subroutine unarySum_whole_omp(va,vb)
		real(DP), allocatable, dimension(:,:,:), intent(inout) :: va
		real(DP), allocatable, dimension(:,:,:), intent(in) :: vb
		integer :: is,ie,js,je,ks,ke
		integer :: i,j,k
		
		is = lbound(va,1)
		ie = ubound(va,1)
		js = lbound(va,2)
		je = ubound(va,2)
		ks = lbound(va,3)
		ke = ubound(va,3)
		
		!$OMP PARALLEL DO DEFAULT(none) &
		!$OMP SHARED(va,vb) &
		!$OMP SHARED(is,ie,js,je,ks,ke) &
		!$OMP PRIVATE(i,j,k)
		do k=ks,ke
			do j=js,je
				do i=is,ie
					va(i,j,k) = va(i,j,k)+vb(i,j,k)
				end do
			end do
		end do
		!$OMP END PARALLEL DO
		
	end subroutine
!========================================================================================!

!========================================================================================!
	subroutine unarySum_section_omp(va,vb,is,ie,js,je,ks,ke)
		real(DP), allocatable, dimension(:,:,:), intent(inout) :: va
		real(DP), allocatable, dimension(:,:,:), intent(in) :: vb
		integer, intent(in) :: is,ie,js,je,ks,ke
		integer :: i,j,k
		
		!$OMP PARALLEL DO DEFAULT(none) &
		!$OMP SHARED(va,vb) &
		!$OMP SHARED(is,ie,js,je,ks,ke) &
		!$OMP PRIVATE(i,j,k)
		do k=ks,ke
			do j=js,je
				do i=is,ie
					va(i,j,k) = va(i,j,k)+vb(i,j,k)
				end do
			end do
		end do
		!$OMP END PARALLEL DO
		
	end subroutine
!========================================================================================!
	
end module ompRoutinesMod


	



