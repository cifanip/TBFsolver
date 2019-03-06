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

	use vfieldMod

	implicit none

	real(DP), parameter :: pi = 4.d0*DATAN(1.d0)
	
	public :: initChFlowVelocity
	public :: initShearVelocity
	public :: init_Bubble_vf
	
	
contains

!========================================================================================!
    subroutine initChFlowVelocity(u,mesh,gmesh)
    	type(vfield), intent(inout) :: u
    	type(grid), intent(in) :: mesh,gmesh
    	real(DP) :: Lx, Ly, Lz
    	integer :: i,j,k,ig,jg,kg,proc,nx,ny,nz
    	real(DP) :: x, y, z, icA, icB
    	type(parFile) :: pfile
    	logical :: found
    	
    	call parFileCTOR(pfile,'initVelocity','specs')
    	call readParameter(pfile,found,'perturbed_parabolic')
    	if (found) then
    		call readParameter(pfile,icA,'icA')
    		icB = icA/10.d0
    	else
    		return
    	end if
    	
    	Lx = gmesh%Lx_
    	Ly = gmesh%Ly_
    	Lz = gmesh%Lz_
    	
    	nx=mesh%nx_
    	ny=mesh%ny_
    	nz=mesh%nz_
    	
    	proc=mesh%ptrMPIC_%rank_
		
		!ux
		do k=u%ux_%ks_,u%ux_%ke_
			do j=u%ux_%js_,u%ux_%je_
				do i=u%ux_%is_,u%ux_%ie_
					
					!global indexes 
					ig = mesh%ptrMPIC_%gCoords_(1,proc)*nx+i
					jg = mesh%ptrMPIC_%gCoords_(2,proc)*ny+j
					kg = mesh%ptrMPIC_%gCoords_(3,proc)*nz+k	
					
					x = gmesh%xf_(ig)
					y = gmesh%yc_(jg)
					z = gmesh%zc_(kg)
					
					u%ux_%f_(i,j,k) = icA*y*(Ly-y) + &
									  icB*cos(2.d0*pi*x/Lx)*sin(2.d0*pi*y/Ly)*sin(2.d0*pi*z/Lz) + &
									  icB*cos(4.d0*pi*x/Lx)*sin(4.d0*pi*y/Ly)*sin(4.d0*pi*z/Lz)
				
				end do
			end do
		end do
		
		!uy
		do k=u%uy_%ks_,u%uy_%ke_
			do j=u%uy_%js_,u%uy_%je_
				do i=u%uz_%is_,u%uz_%ie_
					
					!global indexes 
					ig = mesh%ptrMPIC_%gCoords_(1,proc)*nx+i
					jg = mesh%ptrMPIC_%gCoords_(2,proc)*ny+j
					kg = mesh%ptrMPIC_%gCoords_(3,proc)*nz+k
					
					x = gmesh%xc_(ig)
					y = gmesh%yf_(jg)
					z = gmesh%zc_(kg)
					
					u%uy_%f_(i,j,k) = -(icB*Ly)/(2.d0*Lx)* &
					   				  (sin(2.d0*pi*x/Lx)*(-1.d0+cos(2.d0*pi*y/Ly))*sin(2.d0*pi*z/Lz)+ &
					   				   sin(4.d0*pi*x/Lx)*(-1.d0+cos(4.d0*pi*y/Ly))*sin(4.d0*pi*z/Lz))
				
				end do
			end do
		end do
		
		!uz
		do k=u%uz_%ks_,u%uz_%ke_
			do j=u%uz_%js_,u%uz_%je_
				do i=u%uz_%is_,u%uz_%ie_
				
					!global indexes 
					ig = mesh%ptrMPIC_%gCoords_(1,proc)*nx+i
					jg = mesh%ptrMPIC_%gCoords_(2,proc)*ny+j
					kg = mesh%ptrMPIC_%gCoords_(3,proc)*nz+k
					
					x = gmesh%xc_(ig)
					y = gmesh%yc_(jg)
					z = gmesh%zf_(kg)
					
					u%uz_%f_(i,j,k) = -(icB*Lz)/(2.d0*Lx)* &
					   				  (sin(2.d0*pi*x/Lx)*sin(2.d0*pi*y/Ly)*cos(2.d0*pi*z/Lz)+ &
					   				   sin(4.d0*pi*x/Lx)*sin(4.d0*pi*y/Ly)*cos(4.d0*pi*z/Lz))
				
				end do
			end do
		end do
		
		call updateBoundariesV(u)		
        
    end subroutine
!========================================================================================!

!========================================================================================!
    subroutine initShearVelocity(u,mesh)
    	type(vfield), intent(inout) :: u
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

!========================================================================================!
    subroutine init_Bubble_vf(mesh,cblk,x0,y0,z0,R,nref)
    	type(grid), intent(in) :: mesh
		real(DP), allocatable, dimension(:,:,:), intent(inout) :: cblk
		real(DP), intent(in) :: x0,y0,z0,R
		integer, intent(in) :: nref
		integer :: is,js,ks,ie,je,ke,i,j,k,count,ir,kr,jr
		real(DP), allocatable, dimension(:) :: xcv,ycv,zcv,xfv,yfv,zfv,dxfv,dyfv,dzfv
		real(DP) :: x,y,z,rad,dx,dy,dz,dxref,dyref,dzref,xs,ys,zs,Vex,V
		
		is=lbound(cblk,1)
		ie=ubound(cblk,1)
		js=lbound(cblk,2)
		je=ubound(cblk,2)	
		ks=lbound(cblk,3)
		ke=ubound(cblk,3)
		
		!copy from global mesh
		call allocateArray(xcv,is,ie)
		call allocateArray(ycv,js,je)
		call allocateArray(zcv,ks,ke)
		xcv=mesh%xc_(is:ie)
		ycv=mesh%yc_(js:je)
		zcv=mesh%zc_(ks:ke)
		
		call allocateArray(xfv,is-1,ie)
		call allocateArray(yfv,js-1,je)
		call allocateArray(zfv,ks-1,ke)
		xfv=mesh%xf_(is-1:ie)
		yfv=mesh%yf_(js-1:je)
		zfv=mesh%zf_(ks-1:ke)
		
		call allocateArray(dxfv,is,ie)
		call allocateArray(dyfv,js,je)
		call allocateArray(dzfv,ks,ke)		
		dxfv=mesh%dxf_(is:ie)
		dyfv=mesh%dyf_(js:je)
		dzfv=mesh%dzf_(ks:ke)
		
		!init to zero
		cblk=0.d0
		
		!init gas cells
		do k=ks,ke
			do j=js,je
				do i=is,ie
				
					x = xcv(i)-x0
					y = ycv(j)-y0
					z = zcv(k)-z0
					
					rad= x*x+y*y+z*z
					
					if (rad <= (R*R)) then
						cblk(i,j,k) = 1.d0
					end if
					
				end do
			end do
		end do
		
		!refine vof for interface cells
		do k=ks,ke
			do j=js,je
				do i=is,ie
				
					count = 0
					
					!v1
					x = xfv(i-1)-x0
					y = yfv(j-1)-y0
					z = zfv(k-1)-z0
					
   					rad = x*x+y*y+z*z
   					if (rad <= R*R) then
   						count = count + 1
   					end if 
					
					!v2
					x = xfv(i)-x0
					y = yfv(j-1)-y0
					z = zfv(k-1)-z0
					
   					rad = x*x+y*y+z*z
   					if (rad <= R*R) then
   						count = count + 1
   					end if 
					
					!v3
					x = xfv(i)-x0
					y = yfv(j-1)-y0
					z = zfv(k)-z0
					
   					rad = x*x+y*y+z*z
   					if (rad <= R*R) then
   						count = count + 1
   					end if 
					
					!v4
					x = xfv(i-1)-x0
					y = yfv(j-1)-y0
					z = zfv(k)-z0
					
   					rad = x*x+y*y+z*z
   					if (rad <= R*R) then
   						count = count + 1
   					end if 
					
					!v5
					x = xfv(i-1)-x0
					y = yfv(j)-y0
					z = zfv(k-1)-z0
					
   					rad = x*x+y*y+z*z
   					if (rad <= R*R) then
   						count = count + 1
   					end if 
					
					!v6
					x = xfv(i)-x0
					y = yfv(j)-y0
					z = zfv(k-1)-z0
					
   					rad = x*x+y*y+z*z
   					if (rad <= R*R) then
   						count = count + 1
   					end if 
					
					!v7
					x = xfv(i)-x0
					y = yfv(j)-y0
					z = zfv(k)-z0
					
   					rad = x*x+y*y+z*z
   					if (rad <= R*R) then
   						count = count + 1
   					end if 
					
					!v8
					x = xfv(i-1)-x0
					y = yfv(j)-y0
					z = zfv(k)-z0
					
   					rad = x*x+y*y+z*z
   					if (rad <= R*R) then
   						count = count + 1
   					end if 
   					
   					if ((count >0) .AND. (count < 8)) then 
   					
   						!reset vof interface cell
   						cblk(i,j,k) = 0.d0
   											
   						dx = dxfv(i)
   						dy = dyfv(j)
   						dz = dzfv(k)
   						dxref = dx/nref
   						dyref = dy/nref
   						dzref = dz/nref
   						xs = xfv(i-1)
						ys = yfv(j-1)
						zs = zfv(k-1)
   						
   						do kr=1,nref
   							do jr=1,nref
   								do ir=1,nref

 									x = xs+0.5d0*dxref+(ir-1)*dxref-x0
									y = ys+0.5d0*dyref+(jr-1)*dyref-y0
									z = zs+0.5d0*dzref+(kr-1)*dzref-z0
					
   									rad = x*x+y*y+z*z
   									if (rad <= R*R) then
   										cblk(i,j,k) = cblk(i,j,k)+dxref*dyref*dzref
   									end if  
   									 									
   								end do
   							end do
   						end do
   						
   						cblk(i,j,k) = cblk(i,j,k)/(dx*dy*dz)
   						
   					end if
					
					
				end do
			end do
		end do	


		Vex = (4.d0/3.d0)*pi*R*R*R
		V=0.d0
		do k=ks,ke
			do j=js,je
				do i=is,ie
   					dx = dxfv(i)
   					dy = dyfv(j)
   					dz = dzfv(k)
					V=V+cblk(i,j,k)*dx*dy*dz
				end do
			end do
		end do

		!******** uncomment to print out volume error
		!write(*,'(A,'//s_outputFormat(2:9)//')') &
		!		'VOF bubbles volume error: ', abs(V-Vex)/Vex
        
    end subroutine
!========================================================================================!

end module initialConditionsMod


	



