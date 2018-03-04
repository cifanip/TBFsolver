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

module gridMod
	
	use mpiControlMod
	
	implicit none

	
	! gamma factor hyperbolic profile wall-normal
	real(DP), private, parameter :: s_gamma = 1.4d0
	
	!total volume
	real(DP), protected :: s_Vg

	type, public :: grid

#ifdef MG_MODE
		!class pointer for multi-grids
		type(grid), pointer :: ptrGrid_ => NULL()
#endif	
	
		! domain
		real(DP) :: Lx_, Ly_, Lz_ 
		! domain global
		real(DP) :: Lxg_, Lyg_, Lzg_
		
		! number of cells (excluding the halo)
		integer :: nx_, ny_, nz_ 
		! number of global cells (excluding the halo)
		integer :: nxg_, nyg_, nzg_ 
		
		!global indexes
		integer :: i0g_, i1g_, j0g_, j1g_, k0g_, k1g_ 
		
		!halo dim of the grid (maximum hd of all fields)
		integer :: hd_

		!logical for grid ref 		
		logical, private, dimension(3) :: isDirUnif_
	
		!mesh parFile
		type(parFile), private :: pfile_
		
#ifdef MG_MODE	
		!multi-grid level		
		integer :: level_
#endif	
		
		!keep a pointer to mpiControl
		type(mpiControl), pointer :: ptrMPIC_ => NULL()

		! coordinates of the centroids (global coord)
		real(DP), allocatable, dimension(:) :: xc_
		real(DP), allocatable, dimension(:) :: yc_
		real(DP), allocatable, dimension(:) :: zc_
		! group coord c_
		real(DP), allocatable, dimension(:,:) :: posc_
		
		! coordinates of the faces	(global coord)
		real(DP), allocatable, dimension(:) :: xf_
		real(DP), allocatable, dimension(:) :: yf_
		real(DP), allocatable, dimension(:) :: zf_
		! group coord f_
		real(DP), allocatable, dimension(:,:) :: posf_
		
		! distance between centroids
		real(DP), allocatable, dimension(:) :: dxc_
		real(DP), allocatable, dimension(:) :: dyc_
		real(DP), allocatable, dimension(:) :: dzc_
		! group dc_
		real(DP), allocatable, dimension(:,:) :: dc_
		
		! distance between faces
		real(DP), allocatable, dimension(:) :: dxf_
		real(DP), allocatable, dimension(:) :: dyf_
		real(DP), allocatable, dimension(:) :: dzf_
		! group df_
		real(DP), allocatable, dimension(:,:) :: df_
		
		!cell volume
		real(DP), allocatable, dimension(:,:,:) :: V_
		real(DP), allocatable, dimension(:,:,:) :: Vsx_
		
		contains
		
		final :: delete_grid

	end type
	
#ifdef FAST_MODE
	private :: checkGridSpacing
#endif
	private :: coord_unif_c
	private :: coord_unif_f
	private :: coord_hyper_c
	private :: coord_hyper_f
	private :: metrics_centroids
	private :: metrics_faces
#ifdef MG_MODE
	private :: coarsenGrid
	private :: allocatePtrGrid
	private :: deallocatePtrGrid
#endif
	private :: computeGlobalIndexes
	private :: computeCoordinates
	private :: computeMetrics
	private :: groupCoordAndMetrics
	private :: computeVolume
	private :: totalVolume
	
	public :: gridCTOR
	public :: delete_grid
	public :: decomposeGrid
	public :: broadCastGlobalGrid
#ifdef MG_MODE
	public :: coarsenGrids
	public :: setMgLevels
#endif
	public :: globalIndexesFromAll
	


  	
contains

!========================================================================================!
    subroutine delete_grid(this)
        type(grid), intent(inout) :: this

#ifdef MG_MODE        
        !deallocate ptrGrid
        call deallocatePtrGrid(this)
#endif	
        
    end subroutine
!========================================================================================!

!========================================================================================!
	subroutine gridCTOR(this,mpiCTRL)
		type(grid), intent(out) :: this
		type(mpiControl), intent(in), target :: mpiCTRL

		call parFileCTOR(this%pfile_ ,'mesh','specs')
		
		this%ptrMPIC_ => mpiCTRL
			
		call readParameter(this%pfile_,this%Lx_,'Lx',bcast=.FALSE.)
		call readParameter(this%pfile_,this%Ly_,'Ly',bcast=.FALSE.)
		call readParameter(this%pfile_,this%Lz_,'Lz',bcast=.FALSE.)
		
		call readParameter(this%pfile_,this%nx_,'nx',bcast=.FALSE.)
		call readParameter(this%pfile_,this%ny_,'ny',bcast=.FALSE.)
		call readParameter(this%pfile_,this%nz_,'nz',bcast=.FALSE.)
		
		!grid halo dim
		this%hd_=3
		
		!global = local (calling constructor only on MASTER)
		this%Lxg_ = this%Lx_
		this%Lyg_ = this%Ly_
		this%Lzg_ = this%Lz_
		
		this%nxg_ = this%nx_
		this%nyg_ = this%ny_
		this%nzg_ = this%nz_
		
		!global indexes excluding halo
		this%i0g_ = 1
		this%i1g_ = this%nx_
		this%j0g_ = 1
		this%j1g_ = this%ny_
		this%k0g_ = 1
		this%k1g_ = this%nz_

	
		call readParameter(this%pfile_,this%isDirUnif_(1),'isXunif',bcast=.FALSE.)
		call readParameter(this%pfile_,this%isDirUnif_(2),'isYunif',bcast=.FALSE.)
		call readParameter(this%pfile_,this%isDirUnif_(3),'isZunif',bcast=.FALSE.)
#ifdef FAST_MODE
		call checkGridSpacing(this)
#endif
			
		!init coordinates
		call computeCoordinates(this)
	    
	    !compute metrics
	    call computeMetrics(this)
	    
	    !assemble multi-dim arrays
		call groupCoordAndMetrics(this)
		
		!compute volume
		call computeVolume(this)
		!comute total volume (only internal, no halo region)
		call totalVolume(this)
	    
	    
	end subroutine
!========================================================================================!

!========================================================================================!
    subroutine checkGridSpacing(this)
        type(grid), intent(in) :: this
        
		if (.NOT.(this%isDirUnif_(1))) then
			call mpiAbort('Grid non-uniform along x direction')
		end if
		
		if (.NOT.(this%isDirUnif_(3))) then
			call mpiAbort('Grid non-uniform along z direction')
		end if
		
    end subroutine
!========================================================================================!

!								uniform points distribution
!========================================================================================!
    subroutine coord_unif_c(q,Lg,ncg,i0,gHD)
    	real(DP), allocatable, dimension(:), intent(INOUT) :: q
		real(DP), intent(IN) :: Lg    !global mesh length
		integer, intent(IN) :: ncg 	  !global mesh cells along this dir
		integer, intent(IN) :: i0     !start global index
		integer, intent(in) :: gHD
		real(DP) :: delta
		integer :: i, lb, ub, n
		
		lb = lbound(q,1)
		ub = ubound(q,1)
		
		!start gIndex including halo
		n = i0 - gHD
		
		delta = Lg/ncg
		
		do i=lb,ub
			q(i) = 0.5d0*delta + (n-1)*delta
			n = n + 1
		end do
		
    end subroutine
!========================================================================================!

!========================================================================================!
    subroutine coord_unif_f(q,Lg,ncg,i0,gHD)
    	real(DP), allocatable, dimension(:), intent(INOUT) :: q
		real(DP), intent(IN) :: Lg    !global mesh length
		integer, intent(IN) :: ncg 	  !global mesh cells along this dir
		integer, intent(IN) :: i0     !start global index
		integer, intent(in) :: gHD
		real(DP) :: delta
		integer :: i, lb, ub, n

		lb = lbound(q,1)
		ub = ubound(q,1)
		
		!start gIndex including halo
		n = i0-1-gHD
		
		delta = Lg/ncg
		
		do i=lb,ub
			q(i) = n*delta
			n = n+1
		end do
		
    end subroutine
!========================================================================================!	

!						hyperbolic tangent points distribution
!========================================================================================!
    subroutine coord_hyper_c(q,Lg,ncg,i0,gHD)
    	real(DP), allocatable, dimension(:), intent(INOUT) :: q
		real(DP), intent(IN) :: Lg    !global mesh length
		integer, intent(IN) :: ncg 	  !global mesh cells along this dir
		integer, intent(IN) :: i0     !start global index
		integer, intent(in) :: gHD
		integer :: i, lb, ub, n
		
		lb = lbound(q,1)
		ub = ubound(q,1)
		
		!start gIndex including halo
		n = i0 - gHD
		
		do i=lb,ub
			q(i) = 0.5d0*Lg*(1.d0 + tanh(-s_gamma+2.d0*s_gamma*(n-0.5d0)/ncg)/tanh(s_gamma))
			n = n+1
		end do
		
    end subroutine
!========================================================================================!

!========================================================================================!
    subroutine coord_hyper_f(q,Lg,ncg,i0,gHD)
    	real(DP), allocatable, dimension(:), intent(INOUT) :: q
		real(DP), intent(IN) :: Lg    !global mesh length
		integer, intent(IN) :: ncg 	  !global mesh cells along this dir
		integer, intent(IN) :: i0     !start global index
		integer, intent(in) :: gHD
		integer :: i, lb, ub, n
		
		lb = lbound(q,1)
		ub = ubound(q,1)
		
		!start gIndex including halo
		n = i0-1-gHD
		
		do i=lb,ub
			q(i) = 0.5d0*Lg*(1d0 + tanh(-s_gamma+2*s_gamma*n/ncg)/tanh(s_gamma))
			n = n+1
		end do
		
    end subroutine
!========================================================================================!

!========================================================================================!
    subroutine metrics_centroids(this)
        type(grid), intent(inout) :: this
        integer :: i, xlb, xub, ylb, yub, zlb, zub
        
        !note: example for x: distance between node i and i+1 is dxc_(i+1)

		!x metrics
		xlb = lbound(this%xc_,1)
		xub = ubound(this%xc_,1)
		call allocateArray(this%dxc_,xlb+1,xub)
		do i=xlb+1,xub
			this%dxc_(i) = this%xc_(i)-this%xc_(i-1)
		end do
		
		!y metrics
		ylb = lbound(this%yc_,1)
		yub = ubound(this%yc_,1)
		call allocateArray(this%dyc_,ylb+1,yub)
		do i=ylb+1,yub
			this%dyc_(i) = this%yc_(i)-this%yc_(i-1)
		end do
		
		!z metrics
		zlb = lbound(this%zc_,1)
		zub = ubound(this%zc_,1)
		call allocateArray(this%dzc_,zlb+1,zub)
		do i=zlb+1,zub
			this%dzc_(i) = this%zc_(i)-this%zc_(i-1)
		end do
		
    end subroutine
!========================================================================================!

!========================================================================================!
    subroutine metrics_faces(this)
        type(grid), intent(inout) :: this
        integer :: i, xlb, xub, ylb, yub, zlb, zub
        
        !note: example for x: distance between face i and i+1 is dxf_(i+1)
        
		!x metrics
		xlb = lbound(this%xf_,1)
		xub = ubound(this%xf_,1)
		call allocateArray(this%dxf_,xlb+1,xub)
		do i=xlb+1,xub
			this%dxf_(i) = this%xf_(i)-this%xf_(i-1)
		end do
	
		!y metrics
		ylb = lbound(this%yf_,1)
		yub = ubound(this%yf_,1)		
		call allocateArray(this%dyf_,ylb+1,yub)	
		do i=ylb+1,yub
			this%dyf_(i) = this%yf_(i)-this%yf_(i-1)
		end do
		
		!z metrics
		zlb = lbound(this%zf_,1)
		zub = ubound(this%zf_,1)
		call allocateArray(this%dzf_,zlb+1,zub)
		do i=zlb+1,zub
			this%dzf_(i) = this%zf_(i)-this%zf_(i-1)
		end do
		
		
    end subroutine
!========================================================================================!

!========================================================================================!
    subroutine decomposeGrid(this,gLoc,mpiCTRL)
        type(grid), intent(inout) :: this    		   		!global mesh
        type(grid), intent(inout) :: gLoc					!local mesh
		type(mpiControl), intent(in), target :: mpiCTRL
        real(DP), dimension(3) :: nProcsAxis
        integer :: ierror
        
        
        !keep a pointer to mpiCTRL
		gLoc%ptrMPIC_ => mpiCTRL

		!mesh parFile
		call parFileCTOR(gLoc%pfile_,'mesh','specs')
		
		!mesh geom
		call MPI_BCAST(this%nx_, 1, MPI_INTEGER, 0, mpiCTRL%cartComm_, ierror)
		call MPI_BCAST(this%ny_, 1, MPI_INTEGER, 0, mpiCTRL%cartComm_, ierror)
		call MPI_BCAST(this%nz_, 1, MPI_INTEGER, 0, mpiCTRL%cartComm_, ierror)
		call MPI_BCAST(this%Lx_, 1, MPI_DOUBLE_PRECISION, 0, mpiCTRL%cartComm_, ierror)
		call MPI_BCAST(this%Ly_, 1, MPI_DOUBLE_PRECISION, 0, mpiCTRL%cartComm_, ierror)
		call MPI_BCAST(this%Lz_, 1, MPI_DOUBLE_PRECISION, 0, mpiCTRL%cartComm_, ierror)
		call MPI_BCAST(this%isDirUnif_, 3, MPI_LOGICAL, 0, mpiCTRL%cartComm_, ierror)
		call MPI_BCAST(this%hd_, 1, MPI_INTEGER, 0, mpiCTRL%cartComm_, ierror)
		
		
		nProcsAxis = mpiCTRL%nProcsAxis_
		!local domain length
		gLoc%Lx_ = this%Lx_/nProcsAxis(1)
		gLoc%Ly_ = this%Ly_/nProcsAxis(2)
		gLoc%Lz_ = this%Lz_/nProcsAxis(3)
		
		!keep a copy to global domain length
		gLoc%Lxg_ = this%Lx_
		gLoc%Lyg_ = this%Ly_
		gLoc%Lzg_ = this%Lz_
		
		!set logical array for coord stretching
		gLoc%isDirUnif_ = this%isDirUnif_
		
		!local domain size
		gLoc%nx_ = this%nx_/nProcsAxis(1)
		gLoc%ny_ = this%ny_/nProcsAxis(2)
		gLoc%nz_ = this%nz_/nProcsAxis(3)
		
		!set halo dim
		gLoc%hd_ = this%hd_
		
		if ((gLoc%nx_ < gLoc%hd_) .OR. &
			(gLoc%ny_ < gLoc%hd_) .OR. &
			(gLoc%nz_ < gLoc%hd_)) then
			call mpiABORT('Too coarse grid found in decomposeGrid for the given halo')
		end if
		
		!keep a copy to global domain size
		gLoc%nxg_ = this%nx_
		gLoc%nyg_ = this%ny_
		gLoc%nzg_ = this%nz_
		
		!set global indexes
		call computeGlobalIndexes(gLoc)
		
		!init coordinates
		call computeCoordinates(gLoc)
		
		!compute metrics
		call computeMetrics(gLoc)
		
		!assemble multi-dim arrays
		call groupCoordAndMetrics(gLoc)
		
		!compute volume
		call computeVolume(gLoc)
		call computeVolume_sx(gLoc)
		
		!broadcast total volume
		call MPI_BCAST(s_Vg, 1, MPI_DOUBLE_PRECISION, 0, mpiCTRL%cartComm_, ierror)
		
		!clean up memory
		call deallocateArray(this%V_)
				
        
    end subroutine
!========================================================================================!

#ifdef MG_MODE

!========================================================================================!
	subroutine coarsenGrid(this)
		type(grid), intent(inout) :: this 
        integer, dimension(3) :: procCoord

		procCoord = this%ptrMPIC_%procCoord_

		!allocate coarse grid
		call allocatePtrGrid(this)
		
		!copy member vars
		this%ptrGrid_%Lx_ = this%Lx_
		this%ptrGrid_%Ly_ = this%Ly_
		this%ptrGrid_%Lz_ = this%Lz_
		this%ptrGrid_%Lxg_ = this%Lxg_
		this%ptrGrid_%Lyg_ = this%Lyg_
		this%ptrGrid_%Lzg_ = this%Lzg_	
		this%ptrGrid_%isDirUnif_ = this%isDirUnif_
		this%ptrGrid_%pfile_ = this%pfile_
		this%ptrGrid_%ptrMPIC_ => this%ptrMPIC_
		
		!new grid size
		!local
		this%ptrGrid_%nx_ = this%nx_/2
		this%ptrGrid_%ny_ = this%ny_/2
		this%ptrGrid_%nz_ = this%nz_/2
		!global
		this%ptrGrid_%nxg_ = this%nxg_/2
		this%ptrGrid_%nyg_ = this%nyg_/2
		this%ptrGrid_%nzg_ = this%nzg_/2
		
		!set halo dim
		this%ptrGrid_%hd_ = 2
		
		if ((this%ptrGrid_%nx_ < this%ptrGrid_%hd_) .OR. &
			(this%ptrGrid_%ny_ < this%ptrGrid_%hd_) .OR. &
			(this%ptrGrid_%nz_ < this%ptrGrid_%hd_)) then
			call mpiABORT('Too coarse grid found in coarsenGrid for the given halo')
		end if
		
		!set global indexes
		call computeGlobalIndexes(this%ptrGrid_)
		
		!init coordinates
		call computeCoordinates(this%ptrGrid_)
		
		!compute metrics
		call computeMetrics(this%ptrGrid_)
		
		!assemble multi-dim arrays
		call groupCoordAndMetrics(this%ptrGrid_)
		
		!compute volume
		call computeVolume(this%ptrGrid_)
	    
	end subroutine
!========================================================================================!

!========================================================================================!
    recursive subroutine setMgLevels(this,k)
        type(grid), intent(inout) :: this
        integer, intent(in) :: k
			
		this%level_ = k
		
		if ( associated(this%ptrGrid_) ) then 
			call setMgLevels(this%ptrGrid_,k-1)
		else
			return
		end if
		
        
    end subroutine
!========================================================================================!

!========================================================================================!
    recursive subroutine coarsenGrids(this,n)
        type(grid), intent(inout) :: this
        integer, intent(in) :: n
		
		if (n > 1) then 
			call coarsenGrid(this)
			call coarsenGrids(this%ptrGrid_,n-1)
		else
			return
		end if
		
        
    end subroutine
!========================================================================================!

!========================================================================================!
	subroutine allocatePtrGrid(this) 
		type(grid), intent(inout) :: this 
		integer :: err
		
		
		if (.not. associated(this%ptrGrid_)) then
			
			allocate(this%ptrGrid_,STAT=err)

			if (err /= 0) then
				call mpiABORT('Allocation of ptrGrid failed ') 
			end if
		else
			call mpiABORT('Attempt to allocate an already associated ptrGrid ')
		end if

		
	end subroutine
!========================================================================================!

!========================================================================================!
	subroutine deallocatePtrGrid(this) 
		type(grid), intent(inout) :: this 

        if (associated(this%ptrGrid_)) then
            deallocate(this%ptrGrid_)
        end if
	
	end subroutine
!========================================================================================!

#endif

!========================================================================================!
	subroutine computeGlobalIndexes(this) 
		type(grid), intent(inout) :: this 
		integer, dimension(3) :: procCoord

        procCoord = this%ptrMPIC_%procCoord_
        
		this%i0g_ = procCoord(1)*this%nx_+1
		this%i1g_ = this%i0g_ + this%nx_ - 1		
		this%j0g_ = procCoord(2)*this%ny_+1
		this%j1g_ = this%j0g_ + this%ny_ - 1	
		this%k0g_ = procCoord(3)*this%nz_+1
		this%k1g_ = this%k0g_ + this%nz_ - 1
        
	
	end subroutine
!========================================================================================!

!========================================================================================!
	subroutine computeCoordinates(this) 
		type(grid), intent(inout) :: this 
		integer :: gHD
		
		gHD = this%hd_

		call allocateArray(this%xc_,1-gHD,this%nx_+gHD)
		call allocateArray(this%yc_,1-gHD,this%ny_+gHD)
		call allocateArray(this%zc_,1-gHD,this%nz_+gHD)
		
		call allocateArray(this%xf_,0-gHD,this%nx_+gHD)
		call allocateArray(this%yf_,0-gHD,this%ny_+gHD)
		call allocateArray(this%zf_,0-gHD,this%nz_+gHD)
		
		
		!init coordinates
		!centroids
		if (this%isDirUnif_(1)) then
			call coord_unif_c(this%xc_,this%Lxg_,this%nxg_, this%i0g_, gHD)
		else
			call coord_hyper_c(this%xc_,this%Lxg_,this%nxg_, this%i0g_, gHD)
		end if
		
	    if (this%isDirUnif_(2)) then
	    	call coord_unif_c(this%yc_,this%Lyg_,this%nyg_, this%j0g_, gHD)
	    else
	    	call coord_hyper_c(this%yc_,this%Lyg_,this%nyg_, this%j0g_, gHD)
	    end if
	    
	    if (this%isDirUnif_(3)) then
	    	call coord_unif_c(this%zc_,this%Lzg_,this%nzg_, this%k0g_, gHD)
	    else
	    	call coord_hyper_c(this%zc_,this%Lzg_,this%nzg_, this%k0g_, gHD)
	    end if
	    
	    !faces
	    if (this%isDirUnif_(1)) then
	    	call coord_unif_f(this%xf_,this%Lxg_,this%nxg_, this%i0g_, gHD)
	    else
	    	call coord_hyper_f(this%xf_,this%Lxg_,this%nxg_, this%i0g_, gHD)
	    end if
	    
	    if (this%isDirUnif_(2)) then
	    	call coord_unif_f(this%yf_,this%Lyg_,this%nyg_, this%j0g_, gHD)
	    else
	    	call coord_hyper_f(this%yf_,this%Lyg_,this%nyg_, this%j0g_, gHD)
	    end if
	    
	    if (this%isDirUnif_(3)) then
	    	call coord_unif_f(this%zf_,this%Lzg_,this%nzg_, this%k0g_, gHD)
	    else
	    	call coord_hyper_f(this%zf_,this%Lzg_,this%nzg_, this%k0g_, gHD)
	    end if
        
	
	end subroutine
!========================================================================================!

!========================================================================================!
	subroutine computeMetrics(this) 
		type(grid), intent(inout) :: this 

	    call metrics_centroids(this)
	    call metrics_faces(this)
        
	end subroutine
!========================================================================================!

!========================================================================================!
	subroutine groupCoordAndMetrics(this) 
		type(grid), intent(inout) :: this 
		integer :: nmax
		integer :: nx, ny, nz
		integer :: gHD
		
		gHD = this%hd_
		
		nx = this%nx_
		ny = this%ny_
		nz = this%nz_

	    nmax = max(nx,ny,nz)
	    !collect coords c
	    call allocateArray(this%posc_,1,3,1-gHD,nmax+gHD)
	    this%posc_ = 0.d0
	    this%posc_(1,1-gHD:nx+gHD) = this%xc_
	    this%posc_(2,1-gHD:ny+gHD) = this%yc_
	    this%posc_(3,1-gHD:nz+gHD) = this%zc_
	    
	    !collect coords f
	    call allocateArray(this%posf_,1,3,0-gHD,nmax+gHD)
	    this%posf_ = 0.d0
	    this%posf_(1,0-gHD:nx+gHD) = this%xf_
	    this%posf_(2,0-gHD:ny+gHD) = this%yf_
	    this%posf_(3,0-gHD:nz+gHD) = this%zf_
	    
	    !collect deltas c
	    call allocateArray(this%dc_,1,3,2-gHD,nmax+gHD)
	    this%dc_ = 0.d0
	    this%dc_(1,2-gHD:nx+gHD) = this%dxc_
	    this%dc_(2,2-gHD:ny+gHD) = this%dyc_
	    this%dc_(3,2-gHD:nz+gHD) = this%dzc_
	    
	    !collect deltas f
	    call allocateArray(this%df_,1,3,1-gHD,nmax+gHD)
	    this%df_ = 0.d0
	    this%df_(1,1-gHD:nx+gHD) = this%dxf_
	    this%df_(2,1-gHD:ny+gHD) = this%dyf_
	    this%df_(3,1-gHD:nz+gHD) = this%dzf_
        
	end subroutine
!========================================================================================!

!========================================================================================!
	subroutine globalIndexesFromAll(this,i0g,i1g,j0g,j1g,k0g,k1g,n,tp) 
		type(grid), intent(in) :: this 
		integer, intent(in) :: n
		character(len=2), intent(in) :: tp
		integer, intent(out) :: i0g,i1g,j0g,j1g,k0g,k1g
		integer :: nx, ny, nz
		
		nx = this%nx_
		ny = this%ny_
		nz = this%nz_
		
		!default case = 'cl'
		i0g = this%ptrMPIC_%gCoords_(1,n)*nx+1
		i1g = i0g + nx - 1
				
		j0g = this%ptrMPIC_%gCoords_(2,n)*ny+1
		j1g = j0g + ny - 1
				
		k0g = this%ptrMPIC_%gCoords_(3,n)*nz+1
		k1g = k0g + nz - 1
		

		SELECT CASE (tp)
			CASE('sx')
				i0g = i0g - 1
			CASE('sy')
				j0g = j0g - 1
			CASE('sz')
				k0g = k0g - 1
   			CASE DEFAULT
			END SELECT
	    
        
	end subroutine
!========================================================================================!

!========================================================================================!
	subroutine computeVolume(this) 
		type(grid), intent(inout) :: this 
		integer :: nx, ny, nz
		integer :: i, j, k
		integer :: gHD
		
		gHD = this%hd_
		
		nx = this%nx_
		ny = this%ny_
		nz = this%nz_
		
		call allocateArray(this%V_,1-gHD,nx+gHD,1-gHD,ny+gHD,1-gHD,nz+gHD)

		do k=1-gHD,nz+gHD
			do j=1-gHD,ny+gHD
				do i=1-gHD,nx+gHD
					this%V_(i,j,k) = this%dxf_(i)*this%dyf_(j)*this%dzf_(k)
				end do
			end do
		end do
		
        
	end subroutine
!========================================================================================!

!========================================================================================!
	subroutine computeVolume_sx(this) 
		type(grid), intent(inout) :: this 
		real(DP) :: dx,dy,dz,dV
		integer :: nx, ny, nz
		integer :: i, j, k
		
		!stagg-x volume (no halo)
		
		nx = this%nx_
		ny = this%ny_
		nz = this%nz_
		
		call allocateArray(this%Vsx_,0,nx,1,ny,1,nz)

		do k=1,nz
			do j=1,ny
				do i=0,nx
				
        			dx=this%dxc_(i+1)
        			dy=this%dyf_(j)
        			dz=this%dzf_(k)
        			dV=dx*dy*dz	 
        			
        			this%Vsx_(i,j,k)=dV	
        			
				end do
			end do
		end do
		
        
	end subroutine
!========================================================================================!

!========================================================================================!
	subroutine broadCastGlobalGrid(gmesh,lmesh) 
		type(grid), intent(inout) :: gmesh
		type(grid), intent(in) :: lmesh
		integer :: nx,ny,nz,nmax,gHD
		integer :: ierror
		
		gmesh%ptrMPIC_ => lmesh%ptrMPIC_	
		
		gHD = lmesh%hd_
		
		nx = lmesh%nxg_
		ny = lmesh%nyg_
		nz = lmesh%nzg_
		nmax = max(nx,ny,nz)

		if (.NOT.IS_MASTER) then
			call allocateArray(gmesh%xc_,1-gHD,nx+gHD)
			call allocateArray(gmesh%yc_,1-gHD,ny+gHD)
			call allocateArray(gmesh%zc_,1-gHD,nz+gHD)
			call allocateArray(gmesh%posc_,1,3,1-gHD,nmax+gHD)
		
			call allocateArray(gmesh%xf_,0-gHD,nx+gHD)
			call allocateArray(gmesh%yf_,0-gHD,ny+gHD)
			call allocateArray(gmesh%zf_,0-gHD,nz+gHD)
			call allocateArray(gmesh%posf_,1,3,0-gHD,nmax+gHD)
			
			call allocateArray(gmesh%dxc_,lbound(gmesh%xc_,1)+1,ubound(gmesh%xc_,1))
			call allocateArray(gmesh%dyc_,lbound(gmesh%yc_,1)+1,ubound(gmesh%yc_,1))
			call allocateArray(gmesh%dzc_,lbound(gmesh%zc_,1)+1,ubound(gmesh%zc_,1))
			call allocateArray(gmesh%dc_,1,3,2-gHD,nmax+gHD)

			call allocateArray(gmesh%dxf_,lbound(gmesh%xf_,1)+1,ubound(gmesh%xf_,1))		
			call allocateArray(gmesh%dyf_,lbound(gmesh%yf_,1)+1,ubound(gmesh%yf_,1))			
			call allocateArray(gmesh%dzf_,lbound(gmesh%zf_,1)+1,ubound(gmesh%zf_,1))
			call allocateArray(gmesh%df_,1,3,1-gHD,nmax+gHD)			
			
		end if
		
		! number of global grid points
		gmesh%nxg_=nx
		gmesh%nyg_=ny
		gmesh%nzg_=nz
		
		!coord c
		call MPI_BCAST(gmesh%xc_, size(gmesh%xc_), MPI_DOUBLE_PRECISION, 0, &
					   gmesh%ptrMPIC_%cartComm_, ierror)
		call MPI_BCAST(gmesh%yc_, size(gmesh%yc_), MPI_DOUBLE_PRECISION, 0, &
					   gmesh%ptrMPIC_%cartComm_, ierror)
		call MPI_BCAST(gmesh%zc_, size(gmesh%zc_), MPI_DOUBLE_PRECISION, 0, &
					   gmesh%ptrMPIC_%cartComm_, ierror)
		call MPI_BCAST(gmesh%posc_, size(gmesh%posc_), MPI_DOUBLE_PRECISION, 0, &
					   gmesh%ptrMPIC_%cartComm_, ierror)
					   
		!coord f			   
		call MPI_BCAST(gmesh%xf_, size(gmesh%xf_), MPI_DOUBLE_PRECISION, 0, &
					   gmesh%ptrMPIC_%cartComm_, ierror)
		call MPI_BCAST(gmesh%yf_, size(gmesh%yf_), MPI_DOUBLE_PRECISION, 0, &
					   gmesh%ptrMPIC_%cartComm_, ierror)
		call MPI_BCAST(gmesh%zf_, size(gmesh%zf_), MPI_DOUBLE_PRECISION, 0, &
					   gmesh%ptrMPIC_%cartComm_, ierror)
		call MPI_BCAST(gmesh%posf_, size(gmesh%posf_), MPI_DOUBLE_PRECISION, 0, &
					   gmesh%ptrMPIC_%cartComm_, ierror)
		
		!metrics c
		call MPI_BCAST(gmesh%dxc_, size(gmesh%dxc_), MPI_DOUBLE_PRECISION, 0, &
					   gmesh%ptrMPIC_%cartComm_, ierror)
		call MPI_BCAST(gmesh%dyc_, size(gmesh%dyc_), MPI_DOUBLE_PRECISION, 0, &
					   gmesh%ptrMPIC_%cartComm_, ierror)
		call MPI_BCAST(gmesh%dzc_, size(gmesh%dzc_), MPI_DOUBLE_PRECISION, 0, &
					   gmesh%ptrMPIC_%cartComm_, ierror)
		call MPI_BCAST(gmesh%dc_, size(gmesh%dc_), MPI_DOUBLE_PRECISION, 0, &
					   gmesh%ptrMPIC_%cartComm_, ierror)
		
		!metrics f
		call MPI_BCAST(gmesh%dxf_, size(gmesh%dxf_), MPI_DOUBLE_PRECISION, 0, &
					   gmesh%ptrMPIC_%cartComm_, ierror)
		call MPI_BCAST(gmesh%dyf_, size(gmesh%dyf_), MPI_DOUBLE_PRECISION, 0, &
					   gmesh%ptrMPIC_%cartComm_, ierror)
		call MPI_BCAST(gmesh%dzf_, size(gmesh%dzf_), MPI_DOUBLE_PRECISION, 0, &
					   gmesh%ptrMPIC_%cartComm_, ierror)
		call MPI_BCAST(gmesh%df_, size(gmesh%df_), MPI_DOUBLE_PRECISION, 0, &
					   gmesh%ptrMPIC_%cartComm_, ierror)		
					   
        
	end subroutine
!========================================================================================!

!========================================================================================!
    subroutine totalVolume(this)
        type(grid), intent(in) :: this
        integer :: i,j,k,nx,ny,nz
        
        nx=this%nx_
        ny=this%ny_
        nz=this%nz_
        
        s_Vg=0.d0
        
        do k = 1,nz
        	do j = 1,ny
        		do i = 1,nx
        			s_Vg = s_Vg + this%V_(i,j,k)
        		end do
        	end do
        end do
        
    end subroutine
!========================================================================================!


	
end module gridMod


