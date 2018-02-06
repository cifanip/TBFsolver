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

module initialConditionsMod

	use vectorFieldMod

	implicit none
	
	!initial velocity profile ch flow coefficients
	!real(DP), parameter :: s_icA = 22.5d0
	real(DP), parameter :: s_icA = 1.d0
	real(DP), parameter :: s_icB = s_icA/10.d0
	real(DP), parameter :: pi = 4.d0*DATAN(1.d0)
	
	public :: initChFlowVelocity
	public :: initShearVelocity
	
	
contains

!========================================================================================!
    subroutine initChFlowVelocity(u,mesh)
    	type(vectorField), intent(inout) :: u
    	type(grid), intent(in) :: mesh
    	real(DP) :: Lx, Ly, Lz
    	integer :: i, j, k
    	real(DP) :: x, y, z
    	
    	Lx = mesh%Lx_
    	Ly = mesh%Ly_
    	Lz = mesh%Lz_
		
		!ux
		do k=u%ux_%ks_,u%ux_%ke_
			do j=u%ux_%js_,u%ux_%je_
				do i=u%ux_%is_,u%ux_%ie_
					
					x = mesh%xf_(i)
					y = mesh%yc_(j)
					z = mesh%zc_(k)
					
					u%ux_%f_(i,j,k) = s_icA*y*(Ly-y) + &
									  s_icB*cos(2.d0*pi*x/Lx)*sin(2.d0*pi*y/Ly)*sin(2.d0*pi*z/Lz) + &
									  s_icB*cos(4.d0*pi*x/Lx)*sin(4.d0*pi*y/Ly)*sin(4.d0*pi*z/Lz)
				
				end do
			end do
		end do
		
		!uy
		do k=u%uy_%ks_,u%uy_%ke_
			do j=u%uy_%js_,u%uy_%je_
				do i=u%uz_%is_,u%uz_%ie_
					
					x = mesh%xc_(i)
					y = mesh%yf_(j)
					z = mesh%zc_(k)
					
					u%uy_%f_(i,j,k) = -(s_icB*Ly)/(2.d0*Lx)* &
					   				  (sin(2.d0*pi*x/Lx)*(-1.d0+cos(2.d0*pi*y/Ly))*sin(2.d0*pi*z/Lz)+ &
					   				   sin(4.d0*pi*x/Lx)*(-1.d0+cos(4.d0*pi*y/Ly))*sin(4.d0*pi*z/Lz))
				
				end do
			end do
		end do
		
		!uz
		do k=u%uz_%ks_,u%uz_%ke_
			do j=u%uz_%js_,u%uz_%je_
				do i=u%uz_%is_,u%uz_%ie_
					
					x = mesh%xc_(i)
					y = mesh%yc_(j)
					z = mesh%zf_(k)
					
					u%uz_%f_(i,j,k) = -(s_icB*Lz)/(2.d0*Lx)* &
					   				  (sin(2.d0*pi*x/Lx)*sin(2.d0*pi*y/Ly)*cos(2.d0*pi*z/Lz)+ &
					   				   sin(4.d0*pi*x/Lx)*sin(4.d0*pi*y/Ly)*cos(4.d0*pi*z/Lz))
				
				end do
			end do
		end do
		
        
    end subroutine
!========================================================================================!

!========================================================================================!
    subroutine initShearVelocity(u,mesh)
    	type(vectorField), intent(inout) :: u
    	type(grid), intent(in) :: mesh
    	integer :: i, j, k
    	real(DP) :: x, y, z
    	
		
		!ux
		do k=u%ux_%ks_,u%ux_%ke_
			do j=u%ux_%js_,u%ux_%je_
				do i=u%ux_%is_,u%ux_%ie_
					
					x = mesh%xf_(i)
					y = mesh%yc_(j)
					z = mesh%zc_(k)
					
					u%ux_%f_(i,j,k) = 2.d0*sin(2.d0*pi*y)*sin(pi*x)*sin(pi*x)*sin(2.d0*pi*z)
				
				end do
			end do
		end do
		
		!uy
		do k=u%uy_%ks_,u%uy_%ke_
			do j=u%uy_%js_,u%uy_%je_
				do i=u%uz_%is_,u%uz_%ie_
					
					x = mesh%xc_(i)
					y = mesh%yf_(j)
					z = mesh%zc_(k)
					
					u%uy_%f_(i,j,k) = -sin(2.d0*pi*x)*sin(pi*y)*sin(pi*y)*sin(2.d0*pi*z)
				
				end do
			end do
		end do
		
		!uz
		do k=u%uz_%ks_,u%uz_%ke_
			do j=u%uz_%js_,u%uz_%je_
				do i=u%uz_%is_,u%uz_%ie_
					
					x = mesh%xc_(i)
					y = mesh%yc_(j)
					z = mesh%zf_(k)
					
					u%uz_%f_(i,j,k) = -sin(2.d0*pi*x)*sin(pi*z)*sin(pi*z)*sin(2.d0*pi*y)
				
				end do
			end do
		end do
		
        
    end subroutine
!========================================================================================!

	
end module initialConditionsMod


	



