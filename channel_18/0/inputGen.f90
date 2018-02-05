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


! ************************************************************************************** !
! Example of fields initialisation. The subroutine `writeInternalField' writes the field 
! values on the nodes inside the computational domain. The subroutine `writeBoundaryField' 
! assigns the type of boundary and the boundary value at each of the 6 boundary faces.
! Note: do not run this script for the test-case `channel_18' provided with the code. For 
! this particular case the fields have been already initialised by a precursor simulation.    
! ************************************************************************************** !

program inputGen
	

	implicit none

	integer :: nx, ny, nz, err
	real(8), allocatable, dimension(:,:,:) :: c, p, psi, ux, uy, uz
	real(8), allocatable, dimension(:,:,:) :: phi0x,phi0y,phi0z
	
	
	nx=192
	ny=160
	nz=96
	
	!write c
	allocate(c(nx,ny,nz),STAT=err)
	c = 1.d0
	open(UNIT=50,FILE='c',form='UNFORMATTED',ACCESS='STREAM',STATUS='REPLACE',ACTION='WRITE')
		call writeInternalField(c)
		call writeBoundaryField(1,3,0.d0)
		call writeBoundaryField(2,3,0.d0)
		call writeBoundaryField(3,6,0.d0)
		call writeBoundaryField(4,6,0.d0)
		call writeBoundaryField(5,3,0.d0)
		call writeBoundaryField(6,3,0.d0)
	close(50)
	deallocate(c,STAT=err)

	!write p
	allocate(p(nx,ny,nz),STAT=err)
	p = 0.d0
	open(UNIT=50,FILE='p',form='UNFORMATTED',ACCESS='STREAM',STATUS='REPLACE',ACTION='WRITE')
		call writeInternalField(p)
		call writeBoundaryField(1,3,0.d0)
		call writeBoundaryField(2,3,0.d0)
		call writeBoundaryField(3,5,0.d0)
		call writeBoundaryField(4,5,0.d0)
		call writeBoundaryField(5,3,0.d0)
		call writeBoundaryField(6,3,0.d0)
	close(50)
	deallocate(p,STAT=err)
	
	!write psi
	allocate(psi(nx,ny,nz),STAT=err)
	psi = 0.d0
	open(UNIT=50,FILE='psi',form='UNFORMATTED',ACCESS='STREAM',STATUS='REPLACE',ACTION='WRITE')
		call writeInternalField(psi)
		call writeBoundaryField(1,3,0.d0)
		call writeBoundaryField(2,3,0.d0)
		call writeBoundaryField(3,2,0.d0)
		call writeBoundaryField(4,2,0.d0)
		call writeBoundaryField(5,3,0.d0)
		call writeBoundaryField(6,3,0.d0)
	close(50)
	deallocate(psi,STAT=err)
	
	!write ux
	allocate(ux(0:nx,ny,nz),STAT=err)
	ux = 0.d0
	open(UNIT=50,FILE='ux',form='UNFORMATTED',ACCESS='STREAM',STATUS='REPLACE',ACTION='WRITE')
		call writeInternalField(ux)
		call writeBoundaryField(1,3,0.d0)
		call writeBoundaryField(2,3,0.d0)
		call writeBoundaryField(3,1,0.d0)
		call writeBoundaryField(4,1,0.d0)
		call writeBoundaryField(5,3,0.d0)
		call writeBoundaryField(6,3,0.d0)
	close(50)
	deallocate(ux,STAT=err)
	
	!write uy
	allocate(uy(nx,0:ny,nz),STAT=err)
	uy = 0.d0
	open(UNIT=50,FILE='uy',form='UNFORMATTED',ACCESS='STREAM',STATUS='REPLACE',ACTION='WRITE')
		call writeInternalField(uy)
		call writeBoundaryField(1,3,0.d0)
		call writeBoundaryField(2,3,0.d0)
		call writeBoundaryField(3,1,0.d0)
		call writeBoundaryField(4,1,0.d0)
		call writeBoundaryField(5,3,0.d0)
		call writeBoundaryField(6,3,0.d0)
	close(50)
	deallocate(uy,STAT=err)
	
	!write uz
	allocate(uz(nx,ny,0:nz),STAT=err)
	uz = 0.d0
	open(UNIT=50,FILE='uz',form='UNFORMATTED',ACCESS='STREAM',STATUS='REPLACE',ACTION='WRITE')
		call writeInternalField(uz)
		call writeBoundaryField(1,3,0.d0)
		call writeBoundaryField(2,3,0.d0)
		call writeBoundaryField(3,1,0.d0)
		call writeBoundaryField(4,1,0.d0)
		call writeBoundaryField(5,3,0.d0)
		call writeBoundaryField(6,3,0.d0)
	close(50)
	deallocate(uz,STAT=err)
	
	!write phi0x
	allocate(phi0x(0:nx,ny,nz),STAT=err)
	phi0x = 0.d0
	open(UNIT=50,FILE='phi0x',form='UNFORMATTED',ACCESS='STREAM',STATUS='REPLACE',ACTION='WRITE')
		call writeInternalField(phi0x)
		call writeBoundaryField(1,3,0.d0)
		call writeBoundaryField(2,3,0.d0)
		call writeBoundaryField(3,1,0.d0)
		call writeBoundaryField(4,1,0.d0)
		call writeBoundaryField(5,3,0.d0)
		call writeBoundaryField(6,3,0.d0)
	close(50)
	deallocate(phi0x,STAT=err)
	
	!write phi0y
	allocate(phi0y(nx,0:ny,nz),STAT=err)
	phi0y = 0.d0
	open(UNIT=50,FILE='phi0y',form='UNFORMATTED',ACCESS='STREAM',STATUS='REPLACE',ACTION='WRITE')
		call writeInternalField(phi0y)
		call writeBoundaryField(1,3,0.d0)
		call writeBoundaryField(2,3,0.d0)
		call writeBoundaryField(3,1,0.d0)
		call writeBoundaryField(4,1,0.d0)
		call writeBoundaryField(5,3,0.d0)
		call writeBoundaryField(6,3,0.d0)
	close(50)
	deallocate(phi0y,STAT=err)
	
	!write phi0z
	allocate(phi0z(nx,ny,0:nz),STAT=err)
	phi0z = 0.d0
	open(UNIT=50,FILE='phi0z',form='UNFORMATTED',ACCESS='STREAM',STATUS='REPLACE',ACTION='WRITE')
		call writeInternalField(phi0z)
		call writeBoundaryField(1,3,0.d0)
		call writeBoundaryField(2,3,0.d0)
		call writeBoundaryField(3,1,0.d0)
		call writeBoundaryField(4,1,0.d0)
		call writeBoundaryField(5,3,0.d0)
		call writeBoundaryField(6,3,0.d0)
	close(50)
	deallocate(phi0z,STAT=err)




	
contains

subroutine writeInternalField(q)
	real(8), allocatable, dimension(:,:,:), intent(in) :: q
	integer :: i, j, k, n
	integer :: is, js, ks, ie, je, ke
	
	is = lbound(q,1)
	ie = ubound(q,1)
	js = lbound(q,2)
	je = ubound(q,2)
	ks = lbound(q,3)
	ke = ubound(q,3)
	
	n = (ie-is+1)*(je-js+1)*(ke-ks+1)

	write(50) n
	write(50) q(:,:,:)
	
end subroutine

subroutine writeBoundaryField(bn,bt,bv)
	integer, intent(in) :: bn, bt
	real(8), intent(in) :: bv

	write(50) bn
	write(50) bt
	write(50) bv
		
end subroutine

		
	
end program inputGen