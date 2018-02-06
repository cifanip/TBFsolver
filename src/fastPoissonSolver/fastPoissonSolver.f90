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

module fastPoissonSolverMod
! ************************************************************************************** !
! Part of this module is based on the Poisson solver developed in the code AFiD 
! detailed in:
! van der Poel, Erwin P., et al. "A pencil distributed finite difference code for 
! strongly turbulent wall-bounded flows." Computers & Fluids 116 (2015): 10-16. 
! ************************************************************************************** !
	
	use initialConditionsMod, only: pi
	use scalarFieldMod
	use pencilDecMod
	use, intrinsic :: iso_c_binding
	
	implicit none

	type,public :: fastPoissonSolver
		real(DP), allocatable, private, dimension(:) :: lx,lz
		complex(DP), allocatable, private, dimension(:) :: am,ap,ac
	end type
	
	real(C_DOUBLE), allocatable, dimension(:,:,:) :: ux_ph,uy_ph
	complex(C_DOUBLE_COMPLEX), allocatable, dimension(:,:,:) :: ux_sp,uy_sp,uz_sp
	
    type, bind(C) :: fftw_iodim
        integer(C_INT) n, is, os
    end type fftw_iodim
	
    type(fftw_iodim),dimension(1) :: iodim
    type(fftw_iodim),dimension(2) :: howmany
    
    type(C_PTR) :: guruplan_x_fwd,guruplan_x_bwd
    type(C_PTR) :: guruplan_z_fwd,guruplan_z_bwd
    
    integer, parameter :: FFTW_ESTIMATE=64
    integer, parameter :: FFTW_FORWARD=-1
    integer, parameter :: FFTW_BACKWARD=1
    
    interface
           type(C_PTR) function fftw_plan_guru_dft_r2c(rank,dims, &
             howmany_rank,howmany_dims,in,out,flags) &
             bind(C, name='fftw_plan_guru_dft_r2c')
             import
             integer(C_INT), value :: rank
             type(fftw_iodim), dimension(*), intent(in) :: dims
             integer(C_INT), value :: howmany_rank
             type(fftw_iodim), dimension(*), intent(in) :: howmany_dims
             real(C_DOUBLE), dimension(*), intent(out) :: in
             complex(C_DOUBLE_COMPLEX), dimension(*), intent(out) :: out
             integer(C_INT), value :: flags
           end function fftw_plan_guru_dft_r2c
        
           type(C_PTR) function fftw_plan_guru_dft(rank,dims, &
             howmany_rank,howmany_dims,in,out,sign,flags) &
             bind(C, name='fftw_plan_guru_dft')
             import
             integer(C_INT), value :: rank
             type(fftw_iodim), dimension(*), intent(in) :: dims
             integer(C_INT), value :: howmany_rank
             type(fftw_iodim), dimension(*), intent(in) :: howmany_dims
             complex(C_DOUBLE_COMPLEX), dimension(*), intent(out) :: in
             complex(C_DOUBLE_COMPLEX), dimension(*), intent(out) :: out
             integer(C_INT), value :: sign
             integer(C_INT), value :: flags
           end function fftw_plan_guru_dft

           type(C_PTR) function fftw_plan_guru_dft_c2r(rank,dims, &
             howmany_rank,howmany_dims,in,out,flags) &
             bind(C, name='fftw_plan_guru_dft_c2r')
             import
             integer(C_INT), value :: rank
             type(fftw_iodim), dimension(*), intent(in) :: dims
             integer(C_INT), value :: howmany_rank
             type(fftw_iodim), dimension(*), intent(in) :: howmany_dims
             complex(C_DOUBLE_COMPLEX), dimension(*), intent(out) :: in
             real(C_DOUBLE), dimension(*), intent(out) :: out
             integer(C_INT), value :: flags
           end function fftw_plan_guru_dft_c2r
           
           integer(C_INT) function fftw_init_threads() &
           	 bind(C, name='fftw_init_threads')
           	 import
           end function fftw_init_threads
           
           subroutine fftw_plan_with_nthreads(nthreads) &
           	 bind(C, name='fftw_plan_with_nthreads')
           	 import
           	 integer(C_INT), value :: nthreads
           end subroutine fftw_plan_with_nthreads     
    end interface
	
	private :: buildFFTGuruPlans
	private :: computeEigenvalues
	private :: initOffDiagMatrix
	private :: solveTriDiag
	private :: zeroAveragePressure
	
	public :: fastPoissonSolverCTOR
	public :: solveFPS


contains


!========================================================================================!
	subroutine fastPoissonSolverCTOR(this,mesh,gMesh)
		type(fastPoissonSolver), intent(out) :: this
		type(grid), intent(in) :: mesh,gMesh
		type(mpiControl), pointer :: ptrMPIC
		integer :: nxg,nyg,nzg,init
		
		ptrMPIC => mesh%ptrMPIC_
	
		nxg=mesh%nxg_
		nyg=mesh%nyg_
		nzg=mesh%nzg_
		
  		call initPencils(mesh)
  		
  		call allocateArray(ux_ph,x_pen%idx(1),x_pen%idx(2),x_pen%idx(3),x_pen%idx(4),&
  							     x_pen%idx(5),x_pen%idx(6))
  		call allocateArray(uy_ph,s_pen%idx(1),s_pen%idx(2),s_pen%idx(3),s_pen%idx(4),&
  						         s_pen%idx(5),s_pen%idx(6))
  		call allocateArray(ux_sp,x_pen%idx(1),x_pen%idx(2),x_pen%idx(3),x_pen%idx(4),&
  							     x_pen%idx(5),x_pen%idx(6))
  		call allocateArray(uz_sp,z_pen%idx(1),z_pen%idx(2),z_pen%idx(3),z_pen%idx(4),&
  							     z_pen%idx(5),z_pen%idx(6))
  		call allocateArray(uy_sp,y_pen%idx(1),y_pen%idx(2),y_pen%idx(3),y_pen%idx(4),&
  						         y_pen%idx(5),y_pen%idx(6))
  		
  		!fftw initi
  		init = fftw_init_threads()
  		call buildFFTGuruPlans(nxg,nzg)
  		
  		call computeEigenvalues(this,mesh)
  		call allocateArray(this%ac,1,nyg)
  		call initOffDiagMatrix(this,gMesh)

	end subroutine
!========================================================================================!

!========================================================================================!
	subroutine buildFFTGuruPlans(nxg,nzg)
		integer, intent(in) :: nxg,nzg
		
		!multi-threaded fftw
		call fftw_plan_with_nthreads(N_THREADS)

		!plan c2c
		iodim(1)%n=nzg
   		iodim(1)%is=z_pen%n(1)*z_pen%n(2)
   		iodim(1)%os=z_pen%n(1)*z_pen%n(2)
   		howmany(1)%n=z_pen%n(1)
   		howmany(1)%is=1
   		howmany(1)%os=1
   		howmany(2)%n=z_pen%n(2)
   		howmany(2)%is=z_pen%n(1)
   		howmany(2)%os=z_pen%n(1)
    	guruplan_z_fwd=fftw_plan_guru_dft(1,iodim,2,howmany,uz_sp,uz_sp,&
    								      FFTW_FORWARD,FFTW_ESTIMATE)
    	guruplan_z_bwd=fftw_plan_guru_dft(1,iodim,2,howmany,uz_sp,uz_sp,&
    								      FFTW_BACKWARD,FFTW_ESTIMATE)  
    								        
		if (.not.c_associated(guruplan_z_fwd)) then
			call mpiABORT('guruplan_z_fwd failed ')  
		end if
		if (.not.c_associated(guruplan_z_bwd)) then
			call mpiABORT('guruplan_z_bwd failed ')  
		end if
		
		!plan r2c
		iodim(1)%n=nxg
   		iodim(1)%is=1
   		iodim(1)%os=1
   		howmany(1)%n=x_pen%n(2)
   		howmany(1)%is=x_pen%n(1)
   		howmany(1)%os=x_pen%n(1)
   		howmany(2)%n=x_pen%n(3)
   		howmany(2)%is=x_pen%n(1)*x_pen%n(2)
   		howmany(2)%os=x_pen%n(1)*x_pen%n(2)
    	guruplan_x_fwd=fftw_plan_guru_dft_r2c(1,iodim,2,howmany,ux_ph,ux_sp,FFTW_ESTIMATE)
				
		!plan c2r
    	guruplan_x_bwd=fftw_plan_guru_dft_c2r(1,iodim,2,howmany,ux_sp,ux_ph,FFTW_ESTIMATE)	
    	
		if (.not.c_associated(guruplan_x_bwd)) then
			call mpiABORT('guruplan_x_bwd failed ')  
		end if
		if (.not.c_associated(guruplan_x_fwd)) then
			call mpiABORT('guruplan_x_fwd failed ')  
		end if	
		

	end subroutine
!========================================================================================!

!========================================================================================!
	subroutine computeEigenvalues(this,mesh)
		type(fastPoissonSolver), intent(inout) :: this
		type(grid), intent(in) :: mesh
		real(DP), allocatable, dimension(:) :: vx,vz
		integer :: nxg,nzg,i
		real(DP) :: dx,dz
		
		dx = mesh%xc_(2)-mesh%xc_(1)
		dz = mesh%zc_(2)-mesh%zc_(1)
        
        nxg=mesh%nxg_
        nzg=mesh%nzg_
        
        call allocateArray(this%lx,1,nxg)
        call allocateArray(this%lz,1,nzg)
        
        call allocateArray(vx,1,nxg)
        call allocateArray(vz,1,nzg)
        
		!x
      	do i=1,nxg
        	vx(i)=(i-1)*2.d0*pi
      	enddo
      	do i=1,nxg
        	this%lx(i)=2.d0*(cos(vx(i)/nxg)-1.d0)/((dx*dx))
      	enddo
      	
      	!z
      	do i=1,nzg
        	vz(i)=(i-1)*2.d0*pi
      	enddo
      	do i=1,nzg
        	this%lz(i)=2.d0*(cos(vz(i)/nzg)-1.d0)/((dz*dz))
      	enddo


	end subroutine
!========================================================================================!

!========================================================================================!
	subroutine initOffDiagMatrix(this,mesh)
		type(fastPoissonSolver), intent(inout) :: this
		type(grid), intent(in) :: mesh
		integer :: j,nyg
		
		nyg=mesh%nyg_
		
		call allocateArray(this%am,1,nyg)	
		call allocateArray(this%ap,1,nyg)		
		
		do j=1,nyg
			this%am(j)=1.d0/(mesh%dyc_(j)*mesh%dyf_(j))
		end do
		this%am(1)=0.d0

		do j=1,nyg
			this%ap(j)=1.d0/(mesh%dyc_(j+1)*mesh%dyf_(j))
		end do	
		this%ap(nyg)=0.d0	
		

	end subroutine
!========================================================================================!

!========================================================================================!
	subroutine solveTriDiag(this,f)
		type(fastPoissonSolver), intent(inout) :: this
		complex(DP), allocatable, dimension(:,:,:), intent(inout) :: f
		integer :: i,j,k,n,nyg,info
		complex(DP), allocatable, dimension(:) :: du2,fy,am_scal,ap_scal,ac_scal
		integer, allocatable, dimension(:) :: ipiv
		complex(DP) :: r
		character(9) :: err_msg
	
		n=y_pen%n(2)
		call allocateArray(du2,1,n-2)
		call allocateArray(ipiv,1,n)
		call allocateArray(fy,1,n)
		call allocateArray(am_scal,1,n)
		call allocateArray(ap_scal,1,n)
		call allocateArray(ac_scal,1,n)

		!$OMP PARALLEL DO COLLAPSE(2)	&
		!$OMP DEFAULT(none)	&
		!$OMP SHARED(this,y_pen,f,n)	&
		!$OMP PRIVATE(ac_scal,am_scal,ap_scal)	&
		!$OMP PRIVATE(i,j,k) &
		!$OMP PRIVATE(r,fy,du2,ipiv,info,err_msg)
      	do k=y_pen%idx(5),y_pen%idx(6)
        	do i=y_pen%idx(1),y_pen%idx(2)
        	
        		fy=f(i,:,k)
        		
        		do j=y_pen%idx(3),y_pen%idx(4)
        			r=1.d0/(this%lx(i)+this%lz(k)-this%am(j)-this%ap(j))
        			am_scal(j)=this%am(j)*r
        			ap_scal(j)=this%ap(j)*r
        			ac_scal(j)=1.d0
        			fy(j)=fy(j)*r
        		end do
        		
			call zgttrf(n,am_scal(2),ac_scal,ap_scal(1),du2,ipiv,info)
			
			if (info.ne.0) then
				ac_scal(n)=epsilon(0.d0)
				!write(err_msg,'(I9)') info
				!call mpiAbort('zgttrf solver returned: '//err_msg)
			end if
		
			call zgttrs('N',n,1,am_scal(2),ac_scal,ap_scal(1),du2,ipiv, &
						fy,n,info)
			
			f(i,:,k)=fy		

        end do
      end do
      !$OMP END PARALLEL DO
		
	end subroutine
!========================================================================================!

!========================================================================================!
	subroutine solveFPS(this,p,s)
		type(fastPoissonSolver), intent(inout) :: this
		type(scalarField), intent(inout) :: p
		type(scalarField), intent(in) :: s
		integer :: nx,nz,nxg,nyg,nzg

		nx=p%nx_
		nz=p%nz_
		nxg=p%ptrMesh_%nxg_
		nyg=p%ptrMesh_%nyg_
		nzg=p%ptrMesh_%nzg_
		
		uy_ph=s%f_(1:nx,1:nyg,1:nz)	
		call sPen_2_xPen(uy_ph,ux_ph)

		!fft x-pencil (periodic)
		call dfftw_execute_dft_r2c(guruplan_x_fwd,ux_ph,ux_sp)

		!fft z-pensil (periodic)
		call xPen_2_zPen(ux_sp,uz_sp)
		call dfftw_execute_dft(guruplan_z_fwd,uz_sp,uz_sp)
		
		!transpose wall-normal direction
		call zPen_2_yPen(uz_sp,uy_sp)
		
		!solve tridiagonal system
		call solveTriDiag(this,uy_sp)
		
		!inv-fft z-pencil
		call yPen_2_zPen(uy_sp,uz_sp)
		call dfftw_execute_dft(guruplan_z_bwd,uz_sp,uz_sp)
		
		!inv-fft x-pencil
		call zPen_2_xPen(uz_sp,ux_sp)
		call dfftw_execute_dft_c2r(guruplan_x_bwd,ux_sp,ux_ph)
		
		!scaling factor
		ux_ph=ux_ph/(nxg*nzg)
		
		call xPen_2_sPen(ux_ph,uy_ph)
		
		p%f_(1:nx,1:nyg,1:nz)=uy_ph
		
		call zeroAveragePressure(p)
		call updateBoundaries(p)
		

	end subroutine
!========================================================================================!

!========================================================================================!
    subroutine zeroAveragePressure(p)
    	type(scalarField), intent(inout) :: p
    	type(grid), pointer :: mesh
        integer :: nx,ny,nz
        integer :: i,j,k,ierror
        real(DP) :: av,avg
        
        mesh => p%ptrMesh_
        
        nx=mesh%nx_
        ny=mesh%ny_
        nz=mesh%nz_
        
        av=0.d0
        
 		!$OMP PARALLEL DO DEFAULT(none) &
		!$OMP SHARED(p,mesh) &
		!$OMP SHARED(nx,ny,nz) &
		!$OMP PRIVATE(i,j,k) &
		!$OMP REDUCTION(+:av)
        do k = 1,nz
        	do j = 1,ny
        		do i = 1,nx
        			av = av + p%f_(i,j,k)*mesh%V_(i,j,k)
        		end do
        	end do
        end do
        !$OMP END PARALLEL DO 
        
        call Mpi_Allreduce(av, avg, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
					       mesh%ptrMPIC_%cartComm_, ierror)
			
	    avg=avg/s_Vg
		p%f_=p%f_-avg
        
    end subroutine
!========================================================================================!


end module fastPoissonSolverMod

