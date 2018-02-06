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

module pencilDecMod
	
	use gridMod
	
	implicit none
	
	!send/recv buffers
	real(DP), protected, allocatable, dimension(:) :: buff_sp,buff_xp_r
	complex(DP), protected, allocatable, dimension(:) :: buff_xp_c,buff_yp
	
	!car sub communicators
	integer, protected :: xy_subcomm,xz_subcomm,yz_subcomm
	
	!static copy of n procs
	integer, protected :: s_npx,s_npy,s_npz

	!node distributions among pencils
	integer, protected, allocatable, dimension(:,:) :: distr_nxg_npx,distr_nyg_npx,distr_nzg_npx
	integer, protected, allocatable, dimension(:,:) :: distr_nxg_npz,distr_nyg_npz,distr_nzg_npz
	
	!MPI counts and displacements for MPI_ALLTOALLV
	integer, protected, allocatable, dimension(:) :: SENDCOUNTS_s2x,RECVCOUNTS_s2x,SDISPLS_s2x,RDISPLS_s2x
	integer, protected, allocatable, dimension(:) :: SENDCOUNTS_x2z,RECVCOUNTS_x2z,SDISPLS_x2z,RDISPLS_x2z
	integer, protected, allocatable, dimension(:) :: SENDCOUNTS_z2y,RECVCOUNTS_z2y,SDISPLS_z2y,RDISPLS_z2y
	integer, protected, allocatable, dimension(:) :: SENDCOUNTS_y2z,RECVCOUNTS_y2z,SDISPLS_y2z,RDISPLS_y2z
	integer, protected, allocatable, dimension(:) :: SENDCOUNTS_z2x,RECVCOUNTS_z2x,SDISPLS_z2x,RDISPLS_z2x
	integer, protected, allocatable, dimension(:) :: SENDCOUNTS_x2s,RECVCOUNTS_x2s,SDISPLS_x2s,RDISPLS_x2s

	type :: pencil
		integer, dimension(3) :: n
		integer, dimension(6) :: idx
	end type
	
	type(pencil), protected :: s_pen, x_pen, y_pen, z_pen
	
	private :: equi_distr
	private :: checkProcPencils
	
	public :: initPencils
	public :: sPen_2_xPen
	public :: xPen_2_zPen
	public :: zPen_2_yPen
	public :: yPen_2_zPen
	public :: zPen_2_xPen
	public :: xPen_2_sPen
	
  	
	
contains


!========================================================================================!
	subroutine initPencils(mesh)
		type(grid), intent(in) :: mesh
		type(mpiControl), pointer :: mpic
		integer :: npx,npy,npz,nxg,nyg,nzg,i
		integer :: ierror
		logical, dimension(3) :: REMAIN_DIMS
		integer, dimension(3) :: procCoord_c1,procCoord_c2,procCoord_c3
		
		mpic => mesh%ptrMPIC_
		
		s_npx = mpic%nProcsAxis_(1)
		s_npy = mpic%nProcsAxis_(2)
		s_npz = mpic%nProcsAxis_(3) 
		!local copy
		npx = s_npx
		npy = s_npy
		npz = s_npz
		
		nxg = mesh%nxg_
		nyg = mesh%nyg_
		nzg = mesh%nzg_
		
		call checkProcPencils(nxg,nyg,nzg,npx,npy,npz)
		
		!pencil config. 1
		procCoord_c1(1)=mpic%procCoord_(2)
		procCoord_c1(2)=mpic%procCoord_(1)
		procCoord_c1(3)=mpic%procCoord_(3)
		!pencil config. 2
		procCoord_c2(1)=procCoord_c1(3)
		procCoord_c2(2)=procCoord_c1(2)
		procCoord_c2(3)=procCoord_c1(1)
		!pencil config. 3
		procCoord_c3(1)=procCoord_c2(1)
		procCoord_c3(2)=procCoord_c2(3)
		procCoord_c3(3)=procCoord_c2(2)
		
		!sub-comm
		!xy
		REMAIN_DIMS(1)=.TRUE.
		REMAIN_DIMS(2)=.FALSE.
		REMAIN_DIMS(3)=.FALSE.
		call MPI_CART_SUB(mpic%cartComm_, REMAIN_DIMS, xy_subcomm, IERROR)
		!yz
		call MPI_CART_SUB(mpic%cartComm_, REMAIN_DIMS, yz_subcomm, IERROR)
		!xz
		REMAIN_DIMS(1)=.FALSE.
		REMAIN_DIMS(2)=.FALSE.
		REMAIN_DIMS(3)=.TRUE.
		call MPI_CART_SUB(mpic%cartComm_, REMAIN_DIMS, xz_subcomm, IERROR)
		
		!set node distributions
		!npx
		call equi_distr(nxg,npx,distr_nxg_npx)
		call equi_distr(nyg,npx,distr_nyg_npx)
		call equi_distr(nzg,npx,distr_nzg_npx)
		!npz
		call equi_distr(nxg,npz,distr_nxg_npz)
		call equi_distr(nyg,npz,distr_nyg_npz)
		call equi_distr(nzg,npz,distr_nzg_npz)
		
		!start pencil
		s_pen%idx(1)=distr_nxg_npx(1,mpic%procCoord_(1))
		s_pen%idx(2)=distr_nxg_npx(2,mpic%procCoord_(1))
		s_pen%idx(3)=1
		s_pen%idx(4)=nyg
		s_pen%idx(5)=distr_nzg_npz(1,mpic%procCoord_(3))
		s_pen%idx(6)=distr_nzg_npz(2,mpic%procCoord_(3))
		s_pen%n(1)=s_pen%idx(2)-s_pen%idx(1)+1
		s_pen%n(2)=s_pen%idx(4)-s_pen%idx(3)+1
		s_pen%n(3)=s_pen%idx(6)-s_pen%idx(5)+1
	
		!index x-pencil
		x_pen%idx(1)=1
		x_pen%idx(2)=nxg
		x_pen%idx(3)=distr_nyg_npx(1,procCoord_c1(2))
		x_pen%idx(4)=distr_nyg_npx(2,procCoord_c1(2))
		x_pen%idx(5)=distr_nzg_npz(1,procCoord_c1(3))
		x_pen%idx(6)=distr_nzg_npz(2,procCoord_c1(3))
		x_pen%n(1)=x_pen%idx(2)-x_pen%idx(1)+1
		x_pen%n(2)=x_pen%idx(4)-x_pen%idx(3)+1
		x_pen%n(3)=x_pen%idx(6)-x_pen%idx(5)+1
	
		!index z-pencil
		z_pen%idx(1)=distr_nxg_npz(1,procCoord_c2(1))
		z_pen%idx(2)=distr_nxg_npz(2,procCoord_c2(1))
		z_pen%idx(3)=distr_nyg_npx(1,procCoord_c2(2))
		z_pen%idx(4)=distr_nyg_npx(2,procCoord_c2(2))
		z_pen%idx(5)=1
		z_pen%idx(6)=nzg
		z_pen%n(1)=z_pen%idx(2)-z_pen%idx(1)+1
		z_pen%n(2)=z_pen%idx(4)-z_pen%idx(3)+1
		z_pen%n(3)=z_pen%idx(6)-z_pen%idx(5)+1
	
		!index y-pencil
		y_pen%idx(1)=distr_nxg_npz(1,procCoord_c3(1))
		y_pen%idx(2)=distr_nxg_npz(2,procCoord_c3(1))
		y_pen%idx(3)=1
		y_pen%idx(4)=nyg
		y_pen%idx(5)=distr_nzg_npx(1,procCoord_c3(3))
		y_pen%idx(6)=distr_nzg_npx(2,procCoord_c3(3))
		y_pen%n(1)=y_pen%idx(2)-y_pen%idx(1)+1
		y_pen%n(2)=y_pen%idx(4)-y_pen%idx(3)+1
		y_pen%n(3)=y_pen%idx(6)-y_pen%idx(5)+1
		
	
		!allocate buffers
		!buffers	
		call allocateArray(buff_sp,1,s_pen%n(1)*s_pen%n(2)*s_pen%n(3))
		call allocateArray(buff_xp_r,1,x_pen%n(1)*x_pen%n(2)*x_pen%n(3))
		call allocateArray(buff_xp_c,1,x_pen%n(1)*x_pen%n(2)*x_pen%n(3))
		call allocateArray(buff_yp,1,y_pen%n(1)*y_pen%n(2)*y_pen%n(3))
		
	
		!allocate counts and displs
		!*******************              start to x
		!send
		call allocateArray(SENDCOUNTS_s2x,0,npx-1)
		do i=0,npx-1
			SENDCOUNTS_s2x(i)=s_pen%n(1)*(distr_nyg_npx(2,i)-distr_nyg_npx(1,i)+1)*s_pen%n(3)
		end do
		call allocateArray(SDISPLS_s2x,0,npx-1)
		SDISPLS_s2x(0)=0
		do i=1,npx-1
			SDISPLS_s2x(i)=SDISPLS_s2x(i-1)+SENDCOUNTS_s2x(i-1)
		end do
		!recv
		call allocateArray(RECVCOUNTS_s2x,0,npx-1)
		do i=0,npx-1
			RECVCOUNTS_s2x(i)=(distr_nxg_npx(2,i)-distr_nxg_npx(1,i)+1)*x_pen%n(2)*x_pen%n(3)
		end do
		call allocateArray(RDISPLS_s2x,0,npx-1)
		RDISPLS_s2x(0)=0
		do i=1,npx-1
			RDISPLS_s2x(i)=RDISPLS_s2x(i-1)+RECVCOUNTS_s2x(i-1)
		end do
	
		!*******************              x to z
		!send
		call allocateArray(SENDCOUNTS_x2z,0,npz-1)
		do i=0,npz-1
			SENDCOUNTS_x2z(i)=(distr_nxg_npz(2,i)-distr_nxg_npz(1,i)+1)*x_pen%n(2)*x_pen%n(3)
		end do
		call allocateArray(SDISPLS_x2z,0,npz-1)
		SDISPLS_x2z(0)=0
		do i=1,npz-1
			SDISPLS_x2z(i)=SDISPLS_x2z(i-1)+SENDCOUNTS_x2z(i-1)
		end do
		!recv
		call allocateArray(RECVCOUNTS_x2z,0,npz-1)
		do i=0,npz-1
			RECVCOUNTS_x2z(i)=z_pen%n(1)*z_pen%n(2)*(distr_nzg_npz(2,i)-distr_nzg_npz(1,i)+1)
		end do
		call allocateArray(RDISPLS_x2z,0,npz-1)
		RDISPLS_x2z(0)=0
		do i=1,npz-1
			RDISPLS_x2z(i)=RDISPLS_x2z(i-1)+RECVCOUNTS_x2z(i-1)
		end do
	
		!*******************              z to y
		!send
		call allocateArray(SENDCOUNTS_z2y,0,npx-1)
		do i=0,npx-1
			SENDCOUNTS_z2y(i)=z_pen%n(1)*z_pen%n(2)*(distr_nzg_npx(2,i)-distr_nzg_npx(1,i)+1)
		end do
		call allocateArray(SDISPLS_z2y,0,npx-1)
		SDISPLS_z2y(0)=0
		do i=1,npx-1
			SDISPLS_z2y(i)=SDISPLS_z2y(i-1)+SENDCOUNTS_z2y(i-1)
		end do	
		!recv
		call allocateArray(RECVCOUNTS_z2y,0,npx-1)
		do i=0,npx-1
			RECVCOUNTS_z2y(i)=y_pen%n(1)*(distr_nyg_npx(2,i)-distr_nyg_npx(1,i)+1)*y_pen%n(3)
		end do
		call allocateArray(RDISPLS_z2y,0,npx-1)
		RDISPLS_z2y(0)=0
		do i=1,npx-1
			RDISPLS_z2y(i)=RDISPLS_z2y(i-1)+RECVCOUNTS_z2y(i-1)
		end do
	
		!*******************              y to z
		!send
		call allocateArray(SENDCOUNTS_y2z,0,npx-1)
		call allocateArray(SDISPLS_y2z,0,npx-1)
		SENDCOUNTS_y2z=RECVCOUNTS_z2y
		SDISPLS_y2z=RDISPLS_z2y
		!recv
		call allocateArray(RECVCOUNTS_y2z,0,npx-1)
		call allocateArray(RDISPLS_y2z,0,npx-1)
		RECVCOUNTS_y2z=SENDCOUNTS_z2y
		RDISPLS_y2z=SDISPLS_z2y
	
		!*******************              z to x
		!send
		call allocateArray(SENDCOUNTS_z2x,0,npz-1)
		call allocateArray(SDISPLS_z2x,0,npz-1)
		SENDCOUNTS_z2x=RECVCOUNTS_x2z
		SDISPLS_z2x=RDISPLS_x2z
		!recv
		call allocateArray(RECVCOUNTS_z2x,0,npz-1)
		call allocateArray(RDISPLS_z2x,0,npz-1)
		RECVCOUNTS_z2x=SENDCOUNTS_x2z
		RDISPLS_z2x=SDISPLS_x2z
	
		!*******************              x to s
		!send
		call allocateArray(SENDCOUNTS_x2s,0,npx-1)
		call allocateArray(SDISPLS_x2s,0,npx-1)
		SENDCOUNTS_x2s=RECVCOUNTS_s2x
		SDISPLS_x2s=RDISPLS_s2x
		!recv
		call allocateArray(RECVCOUNTS_x2s,0,npx-1)
		call allocateArray(RDISPLS_x2s,0,npx-1)
		RECVCOUNTS_x2s=SENDCOUNTS_s2x
		RDISPLS_x2s=SDISPLS_s2x

				  
	end subroutine
!========================================================================================!

!========================================================================================!
	subroutine checkProcPencils(nxg,nyg,nzg,npx,npy,npz)
		integer, intent(in) :: nxg,nyg,nzg,npx,npy,npz
		
		if (nyg<npx) then
			call mpiAbort('Too few grid points along y-direction for given decomposition')
		end if
		
		if (nxg<npz) then
			call mpiAbort('Too few grid points along x-direction for given decomposition')
		end if
		
		if (nzg<npx) then
			call mpiAbort('Too few grid points along z-direction for given decomposition')
		end if
	

	end subroutine
!========================================================================================!


!========================================================================================!
	subroutine equi_distr(n,np,distr)
		integer :: n,np
		integer, allocatable, dimension(:,:) :: distr
		real :: r
		integer :: n_rem,np_rem,n1,n2,i
	
		call allocateArray(distr,1,2,0,np-1)
	
		r=ceiling(real(n)/real(np))
	
		!new node indexes
		n1=1
		n2=r  	
	
		distr(1,0)=n1
		distr(2,0)=n2

		!initial remainder	
		n_rem=n-r
		np_rem=np-1
	
		i=0
		do while (n_rem>0)
		
			r=ceiling(real(n_rem)/real(np_rem))
				
			!increment proc index
			i=i+1
			!increment node indexes
			n1=n2+1
			n2=n2+r
		
			distr(1,i)=n1
			distr(2,i)=n2
		
			!update rem node and procs
			n_rem=n_rem-r
			np_rem=np_rem-1
		
		end do
	
	end subroutine
!========================================================================================!

!========================================================================================!
	subroutine sPen_2_xPen(v_sp,v_xp)
		real(DP), allocatable, dimension(:,:,:), intent(in) :: v_sp
		real(DP), allocatable, dimension(:,:,:), intent(inout) :: v_xp
		integer :: i,j,k,q,ierror,count

		!prepare send buff
		count=1
		do q=0,s_npx-1
	
			do k=s_pen%idx(5),s_pen%idx(6)
				do j=distr_nyg_npx(1,q),distr_nyg_npx(2,q)
					do i=s_pen%idx(1),s_pen%idx(2)
						buff_sp(count)=v_sp(i,j,k)
						count = count + 1
					end do
				end do
			end do
	
		end do
	
	
		call MPI_ALLTOALLV(buff_sp, SENDCOUNTS_s2x, SDISPLS_s2x, MPI_DOUBLE_PRECISION,	 &
    		           	   buff_xp_r, RECVCOUNTS_s2x, RDISPLS_s2x, MPI_DOUBLE_PRECISION, &
    		           	   xy_subcomm, IERROR)

	
		!unpack recv buff
		count=1
		do q=0,s_npx-1
	
			do k=x_pen%idx(5),x_pen%idx(6)
				do j=x_pen%idx(3),x_pen%idx(4)
					do i=distr_nxg_npx(1,q),distr_nxg_npx(2,q)
						v_xp(i,j,k)=buff_xp_r(count)
						count = count + 1
					end do
				end do
			end do
	
		end do	

	end subroutine
!========================================================================================!

!========================================================================================!
	subroutine xPen_2_zPen(v_xp,v_zp)
		complex(DP), allocatable, dimension(:,:,:), intent(in) :: v_xp
		complex(DP), allocatable, dimension(:,:,:), intent(inout) :: v_zp
		integer :: i,j,k,q,ierror,count
	
		!prepare send buff
		count=1
		do q=0,s_npz-1
	
			do k=x_pen%idx(5),x_pen%idx(6)
				do j=x_pen%idx(3),x_pen%idx(4)
					do i=distr_nxg_npz(1,q),distr_nxg_npz(2,q)
						buff_xp_c(count)=v_xp(i,j,k)
						count = count + 1
					end do
				end do
			end do
	
		end do
	
	
		call MPI_ALLTOALLV(buff_xp_c, SENDCOUNTS_x2z, SDISPLS_x2z, MPI_DOUBLE_COMPLEX,	&
    		           	   v_zp, RECVCOUNTS_x2z, RDISPLS_x2z, MPI_DOUBLE_COMPLEX, 	    &
    		           	   xz_subcomm, IERROR)	
	
	end subroutine
!========================================================================================!

!========================================================================================!
	subroutine zPen_2_yPen(v_zp,v_yp)
		complex(DP), allocatable, dimension(:,:,:), intent(in) :: v_zp
		complex(DP), allocatable, dimension(:,:,:), intent(inout) :: v_yp
		integer :: i,j,k,q,ierror,count

	
		call MPI_ALLTOALLV(v_zp, SENDCOUNTS_z2y, SDISPLS_z2y, MPI_DOUBLE_COMPLEX,     &
    		          	   buff_yp, RECVCOUNTS_z2y, RDISPLS_z2y, MPI_DOUBLE_COMPLEX,  &
    		          	   yz_subcomm, IERROR)

	
		!unpack recv buff
		count=1
		do q=0,s_npx-1
	
			do k=y_pen%idx(5),y_pen%idx(6)
				do j=distr_nyg_npx(1,q),distr_nyg_npx(2,q)
					do i=y_pen%idx(1),y_pen%idx(2)
						v_yp(i,j,k)=buff_yp(count)
						count = count + 1
					end do
				end do
			end do
	
		end do

	end subroutine
!========================================================================================!

!========================================================================================!
	subroutine yPen_2_zPen(v_yp,v_zp)
		complex(DP), allocatable, dimension(:,:,:), intent(in) :: v_yp
		complex(DP), allocatable, dimension(:,:,:), intent(inout) :: v_zp
		integer :: i,j,k,q,ierror,count
	
		!prepare send buff
		count=1
		do q=0,s_npx-1
	
			do k=y_pen%idx(5),y_pen%idx(6)
				do j=distr_nyg_npx(1,q),distr_nyg_npx(2,q)
					do i=y_pen%idx(1),y_pen%idx(2)
						buff_yp(count)=v_yp(i,j,k)
						count = count + 1
					end do
				end do
			end do
	
		end do	
	

		call MPI_ALLTOALLV(buff_yp, SENDCOUNTS_y2z, SDISPLS_y2z, MPI_DOUBLE_COMPLEX, &
    		           	   v_zp, RECVCOUNTS_y2z, RDISPLS_y2z, MPI_DOUBLE_COMPLEX, 	 &
    		           	   yz_subcomm, IERROR)

	end subroutine
!========================================================================================!

!========================================================================================!
	subroutine zPen_2_xPen(v_zp,v_xp)
		complex(DP), allocatable, dimension(:,:,:), intent(in) :: v_zp
		complex(DP), allocatable, dimension(:,:,:), intent(inout) :: v_xp
		integer :: i,j,k,q,ierror,count
	
		call MPI_ALLTOALLV(v_zp, SENDCOUNTS_z2x, SDISPLS_z2x, MPI_DOUBLE_COMPLEX,	    &
    		           	   buff_xp_c, RECVCOUNTS_z2x, RDISPLS_z2x, MPI_DOUBLE_COMPLEX,  &
    		           	   xz_subcomm, IERROR)

	
		!unpack recv buff
		count=1
		do q=0,s_npz-1
	
			do k=x_pen%idx(5),x_pen%idx(6)
				do j=x_pen%idx(3),x_pen%idx(4)
					do i=distr_nxg_npz(1,q),distr_nxg_npz(2,q)
						v_xp(i,j,k)=buff_xp_c(count)
						count = count + 1
					end do
				end do
			end do
	
		end do	

	end subroutine
!========================================================================================!	

!========================================================================================!
	subroutine xPen_2_sPen(v_xp,v_sp)
		real(DP), allocatable, dimension(:,:,:), intent(in) :: v_xp
		real(DP), allocatable, dimension(:,:,:), intent(inout) :: v_sp
		integer :: i,j,k,q,ierror,count

		!prepare send buff
		count=1
		do q=0,s_npx-1
	
			do k=x_pen%idx(5),x_pen%idx(6)
				do j=x_pen%idx(3),x_pen%idx(4)
					do i=distr_nxg_npx(1,q),distr_nxg_npx(2,q)
						buff_xp_r(count)=v_xp(i,j,k)
						count = count + 1
					end do
				end do
			end do
	
		end do	


		call MPI_ALLTOALLV(buff_xp_r, SENDCOUNTS_x2s, SDISPLS_x2s, MPI_DOUBLE_PRECISION, &
    		           	   buff_sp, RECVCOUNTS_x2s, RDISPLS_x2s, MPI_DOUBLE_PRECISION,   &
    		           	   xy_subcomm, IERROR)

	
		!unpack recv buff
		count=1
		do q=0,s_npx-1
	
			do k=s_pen%idx(5),s_pen%idx(6)
				do j=distr_nyg_npx(1,q),distr_nyg_npx(2,q)
					do i=s_pen%idx(1),s_pen%idx(2)
						v_sp(i,j,k)=buff_sp(count)
						count = count + 1
					end do
				end do
			end do
	
		end do	

	end subroutine
!========================================================================================!	
	
end module pencilDecMod


