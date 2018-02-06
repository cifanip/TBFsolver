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

module interpolationMod

	use kindsMod

	implicit none

!DIR$ IF DEFINED (MG_MODE)	
	public :: restrictionOp2D
	public :: restrictionOp3D
	public :: prolongationOp3D
!DIR$ ENDIF
	public :: cellToVertex
	public :: vertexToCell
	public :: cellToVertexBlock
	public :: vertexToCellBlock

	
contains

! 						2D restriction op: bilinear interp
!DIR$ IF DEFINED (MG_MODE)
!========================================================================================!
	subroutine restrictionOp2D(qf,qc,xf,yf,xc,yc) 
		!c = coarse
		!f = fine
		!q = field; x,y = pos
		! note: the fields and the coords are passed without halo
		real(DP), dimension(:,:), intent(IN) :: qf    
		real(DP), dimension(:,:), intent(INOUT) :: qc 
		real(DP), dimension(:), intent(IN) :: xf, yf   
		real(DP), dimension(:), intent(IN) :: xc, yc 
		real(DP) :: wA, wB, wC, wD, A
		integer :: i, j
		
 		!$OMP PARALLEL DO DEFAULT(none) &
		!$OMP SHARED(qf,qc,xf,yf,xc,yc) &
		!$OMP PRIVATE(i,j) &
		!$OMP PRIVATE(A,wA,wB,wC,wD)		
		do j=1,ubound(qc,2)
			do i=1,ubound(qc,1)
				!compute the weights
				wA = (xf(2*i)-xc(i))*(yc(j)-yf(2*j-1))
				wB = (xc(i)-xf(2*i-1))*(yc(j)-yf(2*j-1))
				wC = (xf(2*i)-xc(i))*(yf(2*j)-yc(j))
				wD = (xc(i)-xf(2*i-1))*(yf(2*j)-yc(j))
				!total area
				A = (xf(2*i)-xf(2*i-1))*(yf(2*j)-yf(2*j-1))
				!interpolate value
				qc(i,j) = ( 									  &  
						    wA*qf(2*i-1,2*j) + wB*qf(2*i,2*j)     &
						  + wC*qf(2*i-1,2*j-1) + wD*qf(2*i,2*j-1) &
						  )/A
			end do
		end do
		!$OMP END PARALLEL DO
		
	end subroutine
!========================================================================================!

! 						3D restriction op: trilinear interp
!========================================================================================!
	subroutine restrictionOp3D(qf,qc,xf,yf,zf,xc,yc,zc) 
		!c = coarse
		!f = fine
		!q = field; x,y,z = pos
		! note: the fields and the coords are passed without halo
		real(DP), dimension(:,:,:), intent(IN) :: qf    
		real(DP), dimension(:,:,:), intent(INOUT) :: qc 
		real(DP), dimension(:), intent(IN) :: xf,yf,zf 
		real(DP), dimension(:), intent(IN) :: xc,yc,zc
		real(DP) :: w1, w2, w3, w4, w5, w6, w7, w8, V
		integer :: i, j, k
		
		
 		!$OMP PARALLEL DO DEFAULT(none) &
		!$OMP SHARED(qf,qc,xf,yf,zf,xc,yc,zc) &
		!$OMP PRIVATE(i,j,k) &
		!$OMP PRIVATE(V,w1,w2,w3,w4,w5,w6,w7,w8)		
		do k=1,ubound(qc,3)
			do j=1,ubound(qc,2)
				do i=1,ubound(qc,1)
				!compute the weights
				w1 = (xf(2*i)-xc(i))*(yf(2*j)-yc(j))*(zf(2*k)-zc(k))
				w2 = (xc(i)-xf(2*i-1))*(yf(2*j)-yc(j))*(zf(2*k)-zc(k))
				w3 = (xf(2*i)-xc(i))*(yc(j)-yf(2*j-1))*(zf(2*k)-zc(k))
				w4 = (xc(i)-xf(2*i-1))*(yc(j)-yf(2*j-1))*(zf(2*k)-zc(k))
				w5 = (xf(2*i)-xc(i))*(yf(2*j)-yc(j))*(zc(k)-zf(2*k-1))
				w6 = (xc(i)-xf(2*i-1))*(yf(2*j)-yc(j))*(zc(k)-zf(2*k-1))
				w7 = (xf(2*i)-xc(i))*(yc(j)-yf(2*j-1))*(zc(k)-zf(2*k-1))
				w8 = (xc(i)-xf(2*i-1))*(yc(j)-yf(2*j-1))*(zc(k)-zf(2*k-1))
				!total volume
				V = (xf(2*i)-xf(2*i-1))*(yf(2*j)-yf(2*j-1))*(zf(2*k)-zf(2*k-1))
				
				!interpolate value
				qc(i,j,k) = ( 									  				&  
						    w1*qf(2*i-1,2*j-1,2*k-1) + w2*qf(2*i,2*j-1,2*k-1)   &
						  + w3*qf(2*i-1,2*j,2*k-1) + w4*qf(2*i,2*j,2*k-1)       &
						  + w5*qf(2*i-1,2*j-1,2*k) + w6*qf(2*i,2*j-1,2*k)       &
						  + w7*qf(2*i-1,2*j,2*k) + w8*qf(2*i,2*j,2*k)			&
						  )/V
				end do
			end do
		end do
		!$OMP END PARALLEL DO
		
	end subroutine
!========================================================================================!

! 						3D prolongation op: trilinear interp
!========================================================================================!
	subroutine prolongationOp3D(qf,qc,xf,yf,zf,xc,yc,zc) 
		!c = coarse
		!f = fine
		!q = field; x,y,z = pos
		! note: the finer field is passed without halo
		!		the corse field is passed entirely 
		real(DP), dimension(:,:,:), intent(INOUT) :: qf    
		real(DP), allocatable, dimension(:,:,:), intent(IN) :: qc 
		real(DP), dimension(:), intent(IN) :: xf,yf,zf 
		real(DP), allocatable, dimension(:), intent(IN) :: xc,yc,zc
		real(DP) :: w1, w2, w3, w4, w5, w6, w7, w8, V
		integer :: i, j, k
		integer :: ih, jh, kh
		
		
 		!$OMP PARALLEL DO DEFAULT(none) &
		!$OMP SHARED(qf,qc,xf,yf,zf,xc,yc,zc) &
		!$OMP PRIVATE(i,j,k,ih,jh,kh) &
		!$OMP PRIVATE(V,w1,w2,w3,w4,w5,w6,w7,w8)		
		do k=1,ubound(qf,3)
			do j=1,ubound(qf,2)
				do i=1,ubound(qf,1)
				
				ih = i/2
				jh = j/2
				kh = k/2
				
				!compute the weights
				w1 = (xc(ih+1)-xf(i))*(yc(jh+1)-yf(j))*(zc(kh+1)-zf(k))
				w2 = (xf(i)-xc(ih))*(yc(jh+1)-yf(j))*(zc(kh+1)-zf(k))
				w3 = (xc(ih+1)-xf(i))*(yf(j)-yc(jh))*(zc(kh+1)-zf(k))
				w4 = (xf(i)-xc(ih))*(yf(j)-yc(jh))*(zc(kh+1)-zf(k))
				w5 = (xc(ih+1)-xf(i))*(yc(jh+1)-yf(j))*(zf(k)-zc(kh))
				w6 = (xf(i)-xc(ih))*(yc(jh+1)-yf(j))*(zf(k)-zc(kh))
				w7 = (xc(ih+1)-xf(i))*(yf(j)-yc(jh))*(zf(k)-zc(kh))
				w8 = (xf(i)-xc(ih))*(yf(j)-yc(jh))*(zf(k)-zc(kh))
				!total volume
				V = (xc(ih+1)-xc(ih))*(yc(jh+1)-yc(jh))*(zc(kh+1)-zc(kh))
				
				!interpolate value
				qf(i,j,k) = ( 									  		  &  
						    w1*qc(ih,jh,kh) + w2*qc(ih+1,jh,kh)   		  &		
						  + w3*qc(ih,jh+1,kh) + w4*qc(ih+1,jh+1,kh)       &
					      + w5*qc(ih,jh,kh+1) + w6*qc(ih+1,jh,kh+1)   	  &		
						  + w7*qc(ih,jh+1,kh+1) + w8*qc(ih+1,jh+1,kh+1)   &	
						  )/V  
					
				end do
			end do
		end do
		!$OMP END PARALLEL DO
		
	end subroutine
!========================================================================================!
!DIR$ ENDIF

! 						3D interpolation to vertces: trilinear interp
!========================================================================================!
	subroutine cellToVertex(q,qv,x,y,z,xv,yv,zv) 
		!v = vertex
		!q = field; x,y,z = pos
		real(DP), allocatable, dimension(:,:,:), intent(IN) :: q    
		real(DP), allocatable, dimension(:,:,:), intent(INOUT) :: qv
		real(DP), allocatable, dimension(:), intent(IN) :: x,y,z
		real(DP), allocatable, dimension(:), intent(IN) :: xv, yv, zv
		real(DP) :: w1, w2, w3, w4, w5, w6, w7, w8, V
		integer :: i, j, k

 		!$OMP PARALLEL DO DEFAULT(none) &
		!$OMP SHARED(q,qv,x,y,z,xv,yv,zv) &
		!$OMP PRIVATE(i,j,k) &
		!$OMP PRIVATE(V,w1,w2,w3,w4,w5,w6,w7,w8)			
		do k=lbound(qv,3),ubound(qv,3)
			do j=lbound(qv,2),ubound(qv,2)
				do i=lbound(qv,1),ubound(qv,1)
				!compute the weights
				w1 = (x(i+1)-xv(i))*(yv(j)-y(j))*(z(k+1)-zv(k))
				w2 = (xv(i)-x(i))*(yv(j)-y(j))*(z(k+1)-zv(k))
				w3 = (x(i+1)-xv(i))*(yv(j)-y(j))*(zv(k)-z(k))
				w4 = (xv(i)-x(i))*(yv(j)-y(j))*(zv(k)-z(k))
				w5 = (x(i+1)-xv(i))*(y(j+1)-yv(j))*(z(k+1)-zv(k))
				w6 = (xv(i)-x(i))*(y(j+1)-yv(j))*(z(k+1)-zv(k))
				w7 = (x(i+1)-xv(i))*(y(j+1)-yv(j))*(zv(k)-z(k))
				w8 = (xv(i)-x(i))*(y(j+1)-yv(j))*(zv(k)-z(k))
				!total volume
				V = (x(i+1)-x(i))*(y(j+1)-y(j))*(z(k+1)-z(k))
				
				!interpolate value
				qv(i,j,k) = ( 										&  
						    w1*q(i,j+1,k) + w2*q(i+1,j+1,k)   		&
						  + w3*q(i,j+1,k+1) + w4*q(i+1,j+1,k+1)     &
						  + w5*q(i,j,k) + w6*q(i+1,j,k) 			&
						  + w7*q(i,j,k+1) + w8*q(i+1,j,k+1)      	&
						  )/V
				end do
			end do
		end do
		!$OMP END PARALLEL DO
		
	end subroutine
!========================================================================================!

! 						vertex average interpolation
!========================================================================================!
	subroutine vertexToCell(q,qv,x,y,z,xv,yv,zv) 
		!v = vertex
		!q = field; x,y,z = pos
		real(DP), allocatable, dimension(:,:,:), intent(INOUT) :: q    
		real(DP), allocatable, dimension(:,:,:), intent(IN) :: qv
		real(DP), allocatable, dimension(:), intent(IN) :: x,y,z
		real(DP), allocatable, dimension(:), intent(IN) :: xv,yv,zv
		real(DP) :: w1, w2, w3, w4, w5, w6, w7, w8, V
		integer :: i, j, k
		
 		!$OMP PARALLEL DO DEFAULT(none) &
		!$OMP SHARED(q,qv,x,y,z,xv,yv,zv) &
		!$OMP PRIVATE(i,j,k) &
		!$OMP PRIVATE(V,w1,w2,w3,w4,w5,w6,w7,w8)		
		do k=1,ubound(qv,3)
			do j=1,ubound(qv,2)
				do i=1,ubound(qv,1)
				!compute the weights
				w1 = (xv(i)-x(i))*(y(j)-yv(j-1))*(zv(k)-z(k))
				w2 = (x(i)-xv(i-1))*(y(j)-yv(j-1))*(zv(k)-z(k))
				w3 = (xv(i)-x(i))*(y(j)-yv(j-1))*(z(k)-zv(k-1))
				w4 = (x(i)-xv(i-1))*(y(j)-yv(j-1))*(z(k)-zv(k-1))
				w5 = (xv(i)-x(i))*(yv(j)-y(j))*(zv(k)-z(k))
				w6 = (x(i)-xv(i-1))*(yv(j)-y(j))*(zv(k)-z(k))
				w7 = (xv(i)-x(i))*(yv(j)-y(j))*(z(k)-zv(k-1))
				w8 = (x(i)-xv(i-1))*(yv(j)-y(j))*(z(k)-zv(k-1))

				!total volume
				V = (xv(i)-xv(i-1))*(yv(j)-yv(j-1))*(zv(k)-zv(k-1))
				
				!interpolate value
				q(i,j,k) = ( 										&  
						    w1*qv(i-1,j,k-1) + w2*qv(i,j,k-1)   	&
						  + w3*qv(i-1,j,k) + w4*qv(i,j,k)   		&
						  + w5*qv(i-1,j-1,k-1) + w6*qv(i,j-1,k-1)   &
						  + w7*qv(i-1,j-1,k) + w8*qv(i,j-1,k)       &
						  )/V
				end do
			end do
		end do
		!$OMP END PARALLEL DO
		
	end subroutine
!========================================================================================!

! 					3D interpolation to vertces (vof block): trilinear interp
!========================================================================================!
	subroutine cellToVertexBlock(q,qv,x,y,z,xv,yv,zv) 
		!v = vertex
		!q = field; x,y,z = pos
		real(DP), allocatable, dimension(:,:,:), intent(IN) :: q    
		real(DP), allocatable, dimension(:,:,:), intent(INOUT) :: qv
		real(DP), pointer, dimension(:), intent(IN) :: x,y,z
		real(DP), pointer, dimension(:), intent(IN) :: xv,yv,zv
!		real(DP), allocatable, dimension(:), intent(IN) :: x,y,z
!		real(DP), allocatable, dimension(:), intent(IN) :: xv,yv,zv
		real(DP) :: w1, w2, w3, w4, w5, w6, w7, w8, V
		integer :: i, j, k
			
		do k=lbound(qv,3),ubound(qv,3)
			do j=lbound(qv,2),ubound(qv,2)
				do i=lbound(qv,1),ubound(qv,1)
				!compute the weights
				w1 = (x(i+1)-xv(i))*(yv(j)-y(j))*(z(k+1)-zv(k))
				w2 = (xv(i)-x(i))*(yv(j)-y(j))*(z(k+1)-zv(k))
				w3 = (x(i+1)-xv(i))*(yv(j)-y(j))*(zv(k)-z(k))
				w4 = (xv(i)-x(i))*(yv(j)-y(j))*(zv(k)-z(k))
				w5 = (x(i+1)-xv(i))*(y(j+1)-yv(j))*(z(k+1)-zv(k))
				w6 = (xv(i)-x(i))*(y(j+1)-yv(j))*(z(k+1)-zv(k))
				w7 = (x(i+1)-xv(i))*(y(j+1)-yv(j))*(zv(k)-z(k))
				w8 = (xv(i)-x(i))*(y(j+1)-yv(j))*(zv(k)-z(k))
				!total volume
				V = (x(i+1)-x(i))*(y(j+1)-y(j))*(z(k+1)-z(k))
				
				!interpolate value
				qv(i,j,k) = ( 										&  
						    w1*q(i,j+1,k) + w2*q(i+1,j+1,k)   		&
						  + w3*q(i,j+1,k+1) + w4*q(i+1,j+1,k+1)     &
						  + w5*q(i,j,k) + w6*q(i+1,j,k) 			&
						  + w7*q(i,j,k+1) + w8*q(i+1,j,k+1)      	&
						  )/V
				end do
			end do
		end do
		
	end subroutine
!========================================================================================!

! 						vertex average interpolation (vof block)
!========================================================================================!
	subroutine vertexToCellBlock(q,qv,x,y,z,xv,yv,zv)
		!v = vertex
		!q = field; x,y,z = pos
		real(DP), allocatable, dimension(:,:,:), intent(INOUT) :: q    
		real(DP), allocatable, dimension(:,:,:), intent(IN) :: qv
		real(DP), pointer, dimension(:), intent(IN) :: x,y,z
		real(DP), pointer, dimension(:), intent(IN) :: xv,yv,zv
!		real(DP), allocatable, dimension(:), intent(IN) :: x,y,z
!		real(DP), allocatable, dimension(:), intent(IN) :: xv,yv,zv
		real(DP) :: w1, w2, w3, w4, w5, w6, w7, w8, V
		integer :: i, j, k
				
		do k=lbound(q,3),ubound(q,3)
			do j=lbound(q,2),ubound(q,2)
				do i=lbound(q,1),ubound(q,1)
				!compute the weights
				w1 = (xv(i)-x(i))*(y(j)-yv(j-1))*(zv(k)-z(k))
				w2 = (x(i)-xv(i-1))*(y(j)-yv(j-1))*(zv(k)-z(k))
				w3 = (xv(i)-x(i))*(y(j)-yv(j-1))*(z(k)-zv(k-1))
				w4 = (x(i)-xv(i-1))*(y(j)-yv(j-1))*(z(k)-zv(k-1))
				w5 = (xv(i)-x(i))*(yv(j)-y(j))*(zv(k)-z(k))
				w6 = (x(i)-xv(i-1))*(yv(j)-y(j))*(zv(k)-z(k))
				w7 = (xv(i)-x(i))*(yv(j)-y(j))*(z(k)-zv(k-1))
				w8 = (x(i)-xv(i-1))*(yv(j)-y(j))*(z(k)-zv(k-1))

				!total volume
				V = (xv(i)-xv(i-1))*(yv(j)-yv(j-1))*(zv(k)-zv(k-1))
				
				!interpolate value
				q(i,j,k) = ( 										&  
						    w1*qv(i-1,j,k-1) + w2*qv(i,j,k-1)   	&
						  + w3*qv(i-1,j,k) + w4*qv(i,j,k)   		&
						  + w5*qv(i-1,j-1,k-1) + w6*qv(i,j-1,k-1)   &
						  + w7*qv(i-1,j-1,k) + w8*qv(i,j-1,k)       &
						  )/V
				end do
			end do
		end do
		
	end subroutine
!========================================================================================!



	
end module interpolationMod


	



