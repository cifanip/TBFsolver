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

module vofBlocksMod

	use initialConditionsMod
	use timeMod
	
	implicit none

	integer, protected :: s_nb
	integer, protected :: s_nblk
	
	logical, allocatable, dimension(:), protected :: s_exchange_g
	logical, allocatable, dimension(:), protected :: s_exchange_b
    
    !off-set block
    integer, parameter :: offset_c = 3
    integer, parameter :: offset_u = 2
    integer, parameter :: offset_st = 1
	
	integer, allocatable, dimension(:,:,:), protected :: s_gbList
	
	!array for initialisation of bubbles
	integer, allocatable, dimension(:,:), protected :: s_idx_init
	real(DP), allocatable, dimension(:,:), protected :: s_pos_init
	
	!lists for bookkeeping
	integer, allocatable, dimension(:,:), protected :: s_blk_data
	integer, allocatable, dimension(:,:), protected :: s_blk_proc

	!comm param
	integer, parameter :: UNPACK_MAX = 1,UNPACK_SUM=2
	integer, parameter :: PACK_BOX_VF = 1
	integer, parameter :: PACK_BOX_K = 2
	integer, parameter :: PACK_BOX_ST = 3
	integer, parameter :: PACK_BOX_U = 4
	
	!init param
	integer, parameter :: TIME_LEVEL_0_BLK=1
	integer, parameter :: UPDATE_BLK=2
	integer, parameter :: REDISTRIBUTION_BLK=3
	integer, parameter :: REINIT_BLK=4
	
	type, public :: vofBlock
		integer :: master, bn
		integer, dimension(6) :: idx
		real(DP), pointer, dimension(:) :: xc => NULL()
		real(DP), pointer, dimension(:) :: yc => NULL()
		real(DP), pointer, dimension(:) :: zc => NULL()
		real(DP), pointer, dimension(:) :: xf => NULL()
		real(DP), pointer, dimension(:) :: yf => NULL()
		real(DP), pointer, dimension(:) :: zf => NULL()
		real(DP), pointer, dimension(:) :: dxc => NULL()
		real(DP), pointer, dimension(:) :: dyc => NULL()
		real(DP), pointer, dimension(:) :: dzc => NULL()
		real(DP), pointer, dimension(:) :: dxf => NULL()
		real(DP), pointer, dimension(:) :: dyf => NULL()
		real(DP), pointer, dimension(:) :: dzf => NULL()
		integer, allocatable, dimension(:) :: idx_mx,idx_my,idx_mz
		real(DP), allocatable, dimension(:,:,:) :: c, cs, c0
		real(DP), allocatable, dimension(:,:,:) :: k
		real(DP), allocatable, dimension(:,:,:) :: stx, sty, stz
		real(DP), allocatable, dimension(:,:,:) :: ux, uy, uz	
		real(DP), allocatable, dimension(:,:,:) :: cv
		real(DP), allocatable, dimension(:,:,:) :: nx, ny, nz
		logical, allocatable, dimension(:,:,:) :: isMixed, isFull
		real(DP), allocatable, dimension(:,:,:) :: q
		real(DP), allocatable, dimension(:,:,:) :: geoFlux, corrFlux, corrTerm
	end type
	
    interface assignment(=)
       module procedure assignBlock
    end interface
	
	type(vofBlock), allocatable, dimension(:) :: vofBlocks
	
	
	public :: vofBlocksCTOR
	public :: grid_2_boxes_u
	public :: boxes_2_grid_f
	public :: boxes_2_grid_vf
	public :: printVOFblocks
	public :: gatherLogicalExchange
	public :: excLists
	public :: updateBlock
	public :: reInitBlockDistribution
	public :: deallocateBlocks
	
	private :: build_VOF_blocks
	private :: initVOFblocks
	private :: measureBlocksDistr
	private :: readVOFblocks
	private :: initBlock
	private :: allocate_blk_arrays
	private :: addNewBlock
	private :: copyBlocks
	private :: assignBlock
	private :: checkBlockSize
	private :: copyBlockMesh
	private :: computeModuloIdx
	private :: bubbles_distribution_space
	private :: init_indexes_box
	private :: readSingleBubble
	private :: readTwoBubbles
	private :: readBubblesArray
	private :: allocateBlockSendRecvBuffers
	private :: packSendToSlaveBuff
	private :: packSendToMasterBuff
	private :: unPackRecvFromMasterBuff
	private :: unPackRecvFromSlaveBuff
	private :: fillUpMasterList
	private :: isBubbleCenterHere
	private :: isBubbleInThisCore
	private :: bubblesLagrQ
	private :: max_box_size

	
contains

!========================================================================================!
    subroutine vofBlocksCTOR(mesh,gmesh,rt)
    	type(grid), intent(in) :: mesh,gmesh
    	type(time), intent(in) :: rt
    	integer :: nprocs
    
		!init blocks
		call build_VOF_blocks(mesh,gmesh,rt)

		nprocs=mesh%ptrMPIC_%nProcs_
    	call allocateArray(s_gbList,1,8,0,nprocs-1,1,s_nb)
    	
    	call reAllocateArray(s_blk_data,1,8,1,s_nblk)
		call reAllocateArray(s_blk_proc,1,s_nblk,0,nprocs-1)
    	
		call excLists(mesh)
		
		call allocateArray(s_exchange_g,0,nprocs-1)
		if (s_nblk>0) then
			call allocateArray(s_exchange_b,1,s_nblk)
		else
			call allocateArray(s_exchange_b,1,1)		
		end if
		s_exchange_g(mesh%ptrMPIC_%rank_)=.FALSE.
		s_exchange_b = .FALSE.	
		
    end subroutine   	
!========================================================================================!

!========================================================================================!
    subroutine build_VOF_blocks(mesh,gmesh,rt)
    	type(grid), intent(in) :: mesh,gmesh
    	type(time), intent(in) :: rt
    	
    	if ((rt%t_>0.d0).OR.(rt%restart_boxes_)) then
    		call readVOFblocks(mesh,gmesh,rt)
    	else
    		call initVOFblocks(mesh,gmesh)
    	end if
    	
    end subroutine   	
!========================================================================================!

!========================================================================================!
    subroutine measureBlocksDistr(mpic)
    	type(mpiControl), intent(in) :: mpic
    	integer :: n_max,n_min,nproc,ierror
    	real(DP) :: av_sq,av_sq_sc,q

		call Mpi_Reduce(s_nblk, n_max, 1, MPI_INTEGER, MPI_MAX, 0, &
					    mpic%cartComm_, ierror)
		call Mpi_Reduce(s_nblk, n_min, 1, MPI_INTEGER, MPI_MIN, 0, &
					    mpic%cartComm_, ierror)
	
		nproc=mpic%nProcs_
		q=s_nblk*s_nblk

		call Mpi_Reduce(q, av_sq, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, &
					    mpic%cartComm_, ierror)		
		
		if (IS_MASTER) then
			av_sq=sqrt(av_sq/nproc)
			av_sq_sc=av_sq/(s_nb/nproc)
		end if  
		
		if (IS_MASTER) then	
			write(*,'(A,I5,A,I5,A,ES11.4E2)') '	Block distr.: MAX = ',n_max, &
										      ' MIN = ',n_min, ' AVG = ',av_sq_sc
		end if	 
		
    end subroutine   	
!========================================================================================!

!========================================================================================!
    subroutine readVOFblocks(mesh,gmesh,rt)
    	type(grid), intent(in) :: mesh,gmesh
    	type(time), intent(in) :: rt
    	type(mpiControl), pointer :: mpic
    	real(DP), allocatable, dimension(:,:,:) :: c,c0
    	integer, dimension(6) :: idx
    	character(len=10) :: dirName,bname
    	real(DP) :: dummy
    	integer :: ierror,nb_rmp,i,master,bn,nb_tmp
    	
    	mpic => mesh%ptrMPIC_
    	
		if (IS_MASTER) then
			write(*,'(A)') 'READ VOF BLOCKS'
		end if
    	
    	write(dirName,s_intFormat) rt%inputFold_
    	
    	!read s_nb
    	if (IS_MASTER) then
        	open(UNIT=s_IOunitNumber,FILE=adjustl(trim(dirName)//'/info_restart'),&
        		 STATUS='old',ACTION='read')
				read(s_IOunitNumber,s_doubleFormat) dummy
				read(s_IOunitNumber,s_intFormat) s_nb
			close(s_IOunitNumber) 
		end if
		call MPI_BCAST(s_nb, 1, MPI_INTEGER, 0, mpic%cartComm_, ierror)
    	
    	do i=1,s_nb
    		
			write(bname,s_intFormat) i
        	open(UNIT=s_IOunitNumber,FILE=trim(adjustl(dirName))//'/b'//trim(adjustl(bname)), &
        	 	      form='UNFORMATTED',ACCESS='STREAM',STATUS='OLD',ACTION='READ') 
        	
        	read(s_IOunitNumber) master
        	
        	if (master==mpic%rank_) then
        		read(s_IOunitNumber) bn
        		read(s_IOunitNumber) idx
        		call reAllocateArray(c,idx(1)-offset_c,idx(2)+offset_c,&
        							   idx(3)-offset_c,idx(4)+offset_c,&
        							   idx(5)-offset_c,idx(6)+offset_c)
        		call reAllocateArray(c0,idx(1)-offset_c,idx(2)+offset_c,&
        							    idx(3)-offset_c,idx(4)+offset_c,&
        							    idx(5)-offset_c,idx(6)+offset_c)
        		read(s_IOunitNumber) c
        		c0=c
        		
        		call addNewBlock(vofBlocks,nb_tmp)
        		call initBlock(mesh,gmesh,vofBlocks(nb_tmp),bn,master,idx(1),idx(2),&
        				       idx(3),idx(4),idx(5),idx(6),REINIT_BLK,c,c0) 
        		
        	end if	
    		
    	end do
    	
		!set number of blocks
    	if (allocated(vofBlocks)) then
			s_nblk = size(vofBlocks)
		else
			s_nblk=0
		end if

		if (IS_MASTER) then
			write(*,'(A)') 'END READ VOF BLOCKS'
		end if
		
		call measureBlocksDistr(mpic)
    	
    end subroutine   	
!========================================================================================!

!========================================================================================!
    subroutine initVOFblocks(mesh,gmesh)
    	type(grid), intent(in) :: mesh,gmesh
    	type(mpiControl), pointer :: mpiCTRL
    	integer, parameter :: single_bubbles = 1, array_bubbles = 2, two_bubbles=3
    	logical, allocatable, dimension(:) :: b_proc_bool
    	type(parFile) :: pfile
    	integer :: method,nref,ierror,bi,is,js,ks,ie,je,ke,nb_tmp
    	real(DP) :: x0,y0,z0,R
    	logical :: present
    	!dummies for initBlock routine
		real(DP), allocatable, dimension(:,:,:) :: dummy_rdp
    	
    	mpiCTRL => mesh%ptrMPIC_
    	
		if (IS_MASTER) then
			write(*,'(A)') 'INIT VOF BLOCKS'
		end if
			

		if (IS_MASTER) then		
		
			call parFileCTOR(pfile,'initBubbles','specs')
			call readParameter(pfile,method,'method',bcast=.FALSE.)
		
			select case(method)
				case(single_bubbles)
					call readSingleBubble(gmesh,pfile)
				case(array_bubbles)
					call readBubblesArray(gmesh,pfile)
				case(two_bubbles)
					call readTwoBubbles(gmesh,pfile)
				case default
			end select
		
		end if
    	
    	!broadcast total number of bubbles
    	call MPI_BCAST(s_nb, 1, MPI_INTEGER, 0, mpiCTRL%cartComm_, ierror)
    	!refinements parameter
		call readParameter(pfile,nref,'nref')

    	!broadcast initial bubble indexes and pos
    	if (.NOT.(IS_MASTER)) then
    		call allocateArray(s_idx_init,1,6,1,s_nb)
    		call allocateArray(s_pos_init,1,4,1,s_nb)
    	end if
    	call MPI_BCAST(s_idx_init, size(s_idx_init), MPI_INTEGER, 0, mpiCTRL%cartComm_, ierror)
    	call MPI_BCAST(s_pos_init, size(s_pos_init), MPI_DOUBLE_PRECISION, 0, &
    			       mpiCTRL%cartComm_, ierror)

		!assign bubble distribution
		call bubbles_distribution_space(mesh,b_proc_bool,mpiCTRL%rank_,s_nb)
    	
    	!init box
    	do bi=1,s_nb
    	
    		present=b_proc_bool(bi)
    		
    		if (present) then
    			
				is=s_idx_init(1,bi)
				ie=s_idx_init(2,bi)
				js=s_idx_init(3,bi)
				je=s_idx_init(4,bi)
				ks=s_idx_init(5,bi)
				ke=s_idx_init(6,bi)
    			
				call addNewBlock(vofBlocks,nb_tmp)
				
    			x0=s_pos_init(1,bi)
    			y0=s_pos_init(2,bi)
    			z0=s_pos_init(3,bi)
    			R=s_pos_init(4,bi)
				
				call initBlock(mesh,gmesh,vofBlocks(nb_tmp),bi,mpiCTRL%rank_,is,ie,js,je,ks,ke,&
    					       TIME_LEVEL_0_BLK,dummy_rdp,dummy_rdp,x0,y0,z0,R,nref)
    					       
    		end if
    	end do
    	    	
		!set number of blocks
    	if (allocated(vofBlocks)) then
			s_nblk = size(vofBlocks)
		else
			s_nblk=0
		end if
    	
    	!clean up	
		call deallocateArray(s_idx_init)
		call deallocateArray(s_pos_init)
		
		!initialise c0 field
		do bi=1,s_nblk
			vofBlocks(bi)%c0=vofBlocks(bi)%c
		end do  

		if (IS_MASTER) then
			write(*,'(A)') 'END INIT VOF BLOCKS'
		end if
		
		call measureBlocksDistr(mpiCTRL)
	
    	
     end subroutine   	
!========================================================================================!

!========================================================================================!
    subroutine initBlock(mesh,gmesh,vofb,bubble_id,box_master,is,ie,js,je,ks,ke,&
    					 INIT_TYPE,tmp_c,tmp_c0,&
    					 x0_opt,y0_opt,z0_opt,R_opt,nref_opt)
   		type(grid), intent(in) :: mesh,gmesh
    	type(vofBlock), intent(inout) :: vofb
    	integer, intent(in) :: bubble_id,box_master
    	integer, intent(inout) :: is,ie,js,je,ks,ke
    	real(DP), intent(in), optional :: x0_opt,y0_opt,z0_opt,R_opt
    	integer, intent(in), optional :: nref_opt
    	integer, intent(in) :: INIT_TYPE
    	real(DP), allocatable, dimension(:,:,:), intent(inout) :: tmp_c,tmp_c0
    	real(DP), allocatable, dimension(:,:,:) :: c_blk
    	real(DP) :: x0,y0,z0,R
    	integer :: nxg,nyg,nzg,nref
		
		nxg = mesh%nxg_
		nyg = mesh%nyg_
		nzg = mesh%nzg_
		
		if (INIT_TYPE==UPDATE_BLK) then		
			!tmp copy vof field
			call reAllocateArray(tmp_c,is,ie,js,je,ks,ke)
			call reAllocateArray(tmp_c0,is,ie,js,je,ks,ke)
			tmp_c = vofb%c(is:ie,js:je,ks:ke)
			tmp_c0 = vofb%c0(is:ie,js:je,ks:ke)
			
			!reset periodic 
			if ((is>nxg).OR.(ie<1)) then
				is = modulo(is-1,nxg)+1
				ie = modulo(ie-1,nxg)+1	
			end if
			if ((js>nyg).OR.(je<1)) then
				js = modulo(js-1,nyg)+1
				je = modulo(je-1,nyg)+1	
			end if
			if ((ks>nzg).OR.(ke<1)) then
				ks = modulo(ks-1,nzg)+1
				ke = modulo(ke-1,nzg)+1	
			end if	
		end if	
							
		vofb%bn=bubble_id
		vofb%master=box_master
		vofb%idx(1)=is
		vofb%idx(2)=ie
		vofb%idx(3)=js
		vofb%idx(4)=je
		vofb%idx(5)=ks
		vofb%idx(6)=ke	
				
		call allocate_blk_arrays(vofb)	
		
		call computeModuloIdx(mesh,vofb,is-offset_c,ie+offset_c,&
								        js-offset_c,je+offset_c,&
								        ks-offset_c,ke+offset_c)

		call copyBlockMesh(mesh,gmesh,vofb)
	
		call checkBlockSize(mesh,vofb)
		
		!reset vf fields
		vofb%c=0.d0
		vofb%c0=0.d0
		
    	select case(INIT_TYPE)
    	
    		case(TIME_LEVEL_0_BLK)		
				x0=x0_opt
				y0=y0_opt
				z0=z0_opt
				R=R_opt
				nref=nref_opt				
			    call reAllocateArray(c_blk,is,ie,js,je,ks,ke)
    			call init_Bubble_vf(gmesh,c_blk,x0,y0,z0,R,nref)
				vofb%c(is:ie,js:je,ks:ke)=c_blk

			case(UPDATE_BLK,REDISTRIBUTION_BLK)
				vofb%c(is:ie,js:je,ks:ke) = tmp_c
				vofb%c0(is:ie,js:je,ks:ke) = tmp_c0
				
			case(REINIT_BLK)
				vofb%c = tmp_c
				vofb%c0 = tmp_c0

			case default
		end select


    end subroutine   	
!========================================================================================!

!========================================================================================!
    subroutine updateBlock(mesh,gmesh,vofb,b)
    	type(grid), intent(in) :: mesh,gmesh
    	type(mpiControl), pointer :: mpic
    	type(vofBlock), intent(inout) :: vofb
    	integer, intent(in) :: b
    	integer :: i,j,k,q
    	integer :: is,ie,js,je,ks,ke
    	integer :: imin,imax,jmin,jmax,kmin,kmax
    	integer :: imin_blk,imax_blk,jmin_blk,jmax_blk,kmin_blk,kmax_blk
		real(DP), allocatable, dimension(:,:,:) :: tmp_c,tmp_c0
    	
    	mpic => mesh%ptrMPIC_
    	
 		is = lbound(vofb%c,1)
		ie = ubound(vofb%c,1)
		js = lbound(vofb%c,2)
		je = ubound(vofb%c,2)
		ks = lbound(vofb%c,3)
		ke = ubound(vofb%c,3) 
			
    	imax=min(is,js,ks)-1
    	jmax=min(is,js,ks)-1
    	kmax=min(is,js,ks)-1
    	imin=ie+je+ke
    	jmin=ie+je+ke
    	kmin=ie+je+ke
    		 		
		do k=ks,ke
			do j=js,je
				do i=is,ie
					if (vofb%isMixed(i,j,k).OR.&
						vofb%isFull(i,j,k)) then
						imax=max(imax,i)
						jmax=max(jmax,j)
						kmax=max(kmax,k)
						imin=min(imin,i)
						jmin=min(jmin,j)
						kmin=min(kmin,k)
					end if
				end do
			end do
		end do

		imax_blk = imax+offset_c
		imin_blk = imin-offset_c
		jmax_blk = jmax+offset_c
		jmin_blk = jmin-offset_c
		kmax_blk = kmax+offset_c
		kmin_blk = kmin-offset_c

			
		if ( (imax_blk /= ie) .OR. &
			 (imin_blk /= is) .OR. &
			 (jmax_blk /= je) .OR. &
			 (jmin_blk /= js) .OR. &
			 (kmax_blk /= ke) .OR. &
			 (kmin_blk /= ks) ) then
				
			 s_exchange_b(b) = .TRUE.
			 
			 call initBlock(mesh,gmesh,vofb,vofb%bn,vofb%master,imin,imax,jmin,jmax,&
			 				kmin,kmax,UPDATE_BLK,tmp_c,tmp_c0) 			    
			 
		else
			s_exchange_b(b) = .FALSE.
		end if

    end subroutine   	
!========================================================================================!

!========================================================================================!
    subroutine allocate_blk_arrays(vofb)
    	type(vofBlock), intent(inout) :: vofb
    	integer :: is,ie,js,je,ks,ke
    	
		is=vofb%idx(1)
		ie=vofb%idx(2)
		js=vofb%idx(3)
		je=vofb%idx(4)
		ks=vofb%idx(5)
		ke=vofb%idx(6)
    	
		call reAllocateArray(vofb%idx_mx,is-offset_c,ie+offset_c)
		call reAllocateArray(vofb%idx_my,js-offset_c,je+offset_c)
		call reAllocateArray(vofb%idx_mz,ks-offset_c,ke+offset_c)
		!mesh
		call reAllocateArray(vofb%xc,is-offset_c,ie+offset_c)
		call reAllocateArray(vofb%yc,js-offset_c,je+offset_c)
		call reAllocateArray(vofb%zc,ks-offset_c,ke+offset_c)
		call reAllocateArray(vofb%xf,is-offset_c-1,ie+offset_c)
		call reAllocateArray(vofb%yf,js-offset_c-1,je+offset_c)
		call reAllocateArray(vofb%zf,ks-offset_c-1,ke+offset_c)
		call reAllocateArray(vofb%dxc,is-offset_c+1,ie+offset_c)
		call reAllocateArray(vofb%dyc,js-offset_c+1,je+offset_c)
		call reAllocateArray(vofb%dzc,ks-offset_c+1,ke+offset_c)
		call reAllocateArray(vofb%dxf,is-offset_c,ie+offset_c)
		call reAllocateArray(vofb%dyf,js-offset_c,je+offset_c)
		call reAllocateArray(vofb%dzf,ks-offset_c,ke+offset_c)
		!fields
		call reAllocateArray(vofb%c,is-offset_c,ie+offset_c,		&
									js-offset_c,je+offset_c,		&
									ks-offset_c,ke+offset_c)
		call reAllocateArray(vofb%c0,is-offset_c,ie+offset_c,		&
									 js-offset_c,je+offset_c,		&
									 ks-offset_c,ke+offset_c)
		call reAllocateArray(vofb%cs,is-offset_u,ie+offset_u,		&
									 js-offset_u,je+offset_u,		&
									 ks-offset_u,ke+offset_u)
		call reAllocateArray(vofb%ux,is-offset_u-1,ie+offset_u,		&
									 js-offset_u,je+offset_u,		&
									 ks-offset_u,ke+offset_u)
		call reAllocateArray(vofb%uy,is-offset_u,ie+offset_u,		&
									 js-offset_u-1,je+offset_u,		&
									 ks-offset_u,ke+offset_u)
		call reAllocateArray(vofb%uz,is-offset_u,ie+offset_u,		&	
									 js-offset_u,je+offset_u,		&
									 ks-offset_u-1,ke+offset_u)
		call reAllocateArray(vofb%stx,is-offset_st-1,ie+offset_st,  &
									 js-offset_st,je+offset_st,		&
									 ks-offset_st,ke+offset_st)
		call reAllocateArray(vofb%sty,is-offset_st,ie+offset_st,	&
									 js-offset_st-1,je+offset_st, 	&
									 ks-offset_st,ke+offset_st)
		call reAllocateArray(vofb%stz,is-offset_st,ie+offset_st,	&
									 js-offset_st,je+offset_st,		&
									 ks-offset_st-1,ke+offset_st)
		call reAllocateArray(vofb%k,is-1,ie+1,js-1,je+1,ks-1,ke+1)		
		call reAllocateArray(vofb%nx,is-1,ie+1,js-1,je+1,ks-1,ke+1)
		call reAllocateArray(vofb%ny,is-1,ie+1,js-1,je+1,ks-1,ke+1)
		call reAllocateArray(vofb%nz,is-1,ie+1,js-1,je+1,ks-1,ke+1)
		call reAllocateArray(vofb%isMixed,is-offset_c,ie+offset_c,	&
										  js-offset_c,je+offset_c,	&
										  ks-offset_c,ke+offset_c)
		call reAllocateArray(vofb%isFull,is-offset_c,ie+offset_c,	&
										 js-offset_c,je+offset_c,	&
										 ks-offset_c,ke+offset_c)
		call reAllocateArray(vofb%cv,is-offset_u-1,ie+offset_u,		&
									 js-offset_u-1,je+offset_u,		&
									 ks-offset_u-1,ke+offset_u)
		call reAllocateArray(vofb%q,is-1,ie+1,js-1,je+1,ks-1,ke+1)
		call reAllocateArray(vofb%geoFlux,is-offset_u,ie+offset_u,	&
										  js-offset_u,je+offset_u,	&
										  ks-offset_u,ke+offset_u)
		call reAllocateArray(vofb%corrFlux,is-offset_u,ie+offset_u,	&
										   js-offset_u,je+offset_u,	&
										   ks-offset_u,ke+offset_u)
		call reAllocateArray(vofb%corrTerm,is-offset_u,ie+offset_u,	&
										   js-offset_u,je+offset_u,	&
										   ks-offset_u,ke+offset_u)   	

     end subroutine   	
!========================================================================================!

!========================================================================================!
    subroutine addNewBlock(blocks,nb)
    	type(vofBlock), allocatable, dimension(:), intent(inout) :: blocks
    	integer, intent(out) :: nb 
    	type(vofBlock), allocatable, dimension(:) :: tmp
    	integer :: i
    		
		if (allocated(blocks)) then
		    nb = size(blocks)
			call copyBlocks(tmp,blocks)
			call deallocateBlocks(blocks)
			allocate(blocks(nb+1))
			do i=1,nb
    				blocks(i)=tmp(i)
    		end do 
			nb=nb+1
		else
			allocate(blocks(1))
			nb=1
		end if
		
		if (allocated(tmp)) then
			call deallocateBlocks(tmp)
		end if
									
    end subroutine   	
!========================================================================================!

!========================================================================================!
    subroutine copyBlocks(blks1,blks2)
    	!copy from 2 to 1
    	type(vofBlock), allocatable, dimension(:), intent(inout) :: blks1
    	type(vofBlock), allocatable, dimension(:), intent(in) :: blks2
    	integer :: i,nb1,nb2
    	
    	if (allocated(blks1)) then
    	
    		nb1 = size(blks1)
    		nb2 = size(blks2)
    		
    		if (nb1==nb2) then
    			do i=1,nb2
    				blks1(i)=blks2(i)
    			end do    			
    		else
    			call deallocateBlocks(blks1)
    			allocate(blks1(nb2))
      			do i=1,nb2
    				blks1(i)=blks2(i)
    			end do   			
    		end if
    		
    	else
    	
    		nb2 = size(blks2)
    		allocate(blks1(nb2))	
    		do i=1,nb2
    			blks1(i)=blks2(i)
    		end do
    		
    	end if
									
    end subroutine   	
!========================================================================================!

!========================================================================================!
    subroutine assignBlock(lhs,rhs)
    	type(vofBlock), intent(inout) :: lhs
    	type(vofBlock), intent(in) :: rhs
    	integer :: lb,ub

		lhs%master=rhs%master
		lhs%bn=rhs%bn
		lhs%idx=rhs%idx
		
		call assignArray(lhs%idx_mx,rhs%idx_mx)
		call assignArray(lhs%idx_my,rhs%idx_my)
		call assignArray(lhs%idx_mz,rhs%idx_mz)
		call assignArray(lhs%c,rhs%c)
		call assignArray(lhs%c0,rhs%c0)
		call assignArray(lhs%cs,rhs%cs)
		call assignArray(lhs%ux,rhs%ux)
		call assignArray(lhs%uy,rhs%uy)
		call assignArray(lhs%uz,rhs%uz)
		call assignArray(lhs%stx,rhs%stx)
		call assignArray(lhs%sty,rhs%sty)
		call assignArray(lhs%stz,rhs%stz)
		call assignArray(lhs%nx,rhs%nx)
		call assignArray(lhs%ny,rhs%ny)
		call assignArray(lhs%nz,rhs%nz)	
		call assignArray(lhs%isMixed,rhs%isMixed)
		call assignArray(lhs%isFull,rhs%isFull)
		call assignArray(lhs%cv,rhs%cv)
		call assignArray(lhs%q,rhs%q)
		call assignArray(lhs%k,rhs%k)
		call assignArray(lhs%geoFlux,rhs%geoFlux)
		call assignArray(lhs%corrFlux,rhs%corrFlux)
		call assignArray(lhs%corrTerm,rhs%corrTerm)

		lb=lbound(rhs%xc,1)
		ub=ubound(rhs%xc,1)
		call reAllocateArray(lhs%xc,lb,ub)
		lhs%xc=rhs%xc
		
		lb=lbound(rhs%yc,1)
		ub=ubound(rhs%yc,1)
		call reAllocateArray(lhs%yc,lb,ub)
		lhs%yc=rhs%yc
		
		lb=lbound(rhs%zc,1)
		ub=ubound(rhs%zc,1)
		call reAllocateArray(lhs%zc,lb,ub)
		lhs%zc=rhs%zc
		
		lb=lbound(rhs%xf,1)
		ub=ubound(rhs%xf,1)
		call reAllocateArray(lhs%xf,lb,ub)
		lhs%xf=rhs%xf
		
		lb=lbound(rhs%yf,1)
		ub=ubound(rhs%yf,1)
		call reAllocateArray(lhs%yf,lb,ub)
		lhs%yf=rhs%yf
		
		lb=lbound(rhs%zf,1)
		ub=ubound(rhs%zf,1)
		call reAllocateArray(lhs%zf,lb,ub)
		lhs%zf=rhs%zf
		
		lb=lbound(rhs%dxc,1)
		ub=ubound(rhs%dxc,1)
		call reAllocateArray(lhs%dxc,lb,ub)
		lhs%dxc=rhs%dxc
		
		lb=lbound(rhs%dyc,1)
		ub=ubound(rhs%dyc,1)
		call reAllocateArray(lhs%dyc,lb,ub)
		lhs%dyc=rhs%dyc
		
		lb=lbound(rhs%dzc,1)
		ub=ubound(rhs%dzc,1)
		call reAllocateArray(lhs%dzc,lb,ub)
		lhs%dzc=rhs%dzc
		
		lb=lbound(rhs%dxf,1)
		ub=ubound(rhs%dxf,1)
		call reAllocateArray(lhs%dxf,lb,ub)
		lhs%dxf=rhs%dxf
		
		lb=lbound(rhs%dyf,1)
		ub=ubound(rhs%dyf,1)
		call reAllocateArray(lhs%dyf,lb,ub)
		lhs%dyf=rhs%dyf
		
		lb=lbound(rhs%dzf,1)
		ub=ubound(rhs%dzf,1)
		call reAllocateArray(lhs%dzf,lb,ub)
		lhs%dzf=rhs%dzf
    						
    end subroutine   	
!========================================================================================!

!========================================================================================!
    subroutine deallocateBlocks(blks)
    	type(vofBlock), allocatable, dimension(:), intent(inout) :: blks
    	integer :: i,nb

		if (allocated(blks)) then

			nb=size(blks)
			
			do i=1,nb
				call deallocateArray(blks(i)%xc)
				call deallocateArray(blks(i)%yc)
				call deallocateArray(blks(i)%zc)
				call deallocateArray(blks(i)%xf)
				call deallocateArray(blks(i)%yf)
				call deallocateArray(blks(i)%zf)
				call deallocateArray(blks(i)%dxc)
				call deallocateArray(blks(i)%dyc)
				call deallocateArray(blks(i)%dzc)
				call deallocateArray(blks(i)%dxf)
				call deallocateArray(blks(i)%dyf)
				call deallocateArray(blks(i)%dzf)	
			end do	
		
			deallocate(blks)
		
		end if
									
    end subroutine   	
!========================================================================================!

!========================================================================================!
    subroutine checkBlockSize(mesh,vofb)
    	type(grid), intent(in) :: mesh
    	type(vofBlock), intent(in) :: vofb
    	integer :: is,ie,js,je,ks,ke
		integer :: nxg,nyg,nzg,sx,sy,sz
		logical :: wrap_x,wrap_y,wrap_z

		
		is=vofb%idx(1)
		ie=vofb%idx(2)
		js=vofb%idx(3)
		je=vofb%idx(4)
		ks=vofb%idx(5)
		ke=vofb%idx(6)		
    	
		nxg = mesh%nxg_
		nyg = mesh%nyg_
		nzg = mesh%nzg_
		
		wrap_x = mesh%ptrMPIC_%wrapAround_(1)
		wrap_y = mesh%ptrMPIC_%wrapAround_(2)
		wrap_z = mesh%ptrMPIC_%wrapAround_(3)
	
		!check box size
		sx=ie-is+1
		sy=je-js+1
		sz=ke-ks+1

		if ((wrap_x).AND.(sx>=nxg)) then
			call mpiABORT('vofBlock size along x is too big ') 
		end if
		if ((wrap_y).AND.(sy>=nyg)) then
			call mpiABORT('vofBlock size along y is too big ') 
		end if
		if ((wrap_z).AND.(sz>=nzg)) then
			call mpiABORT('vofBlock size along z is too big ') 
		end if
		
    end subroutine   	
!========================================================================================!

!========================================================================================!
    subroutine copyBlockMesh(mesh,gmesh,vofb)
    	type(grid), intent(in) :: mesh,gmesh
    	type(vofBlock), intent(inout) :: vofb
    	integer :: i,n,lbi,ubi,lbj,ubj,lbk,ubk,sw
    	integer :: nxg,nyg,nzg
    	real(DP) :: Lx,Ly,Lz,r
    	logical :: wrap_x,wrap_y,wrap_z
    	
    	
    	nxg=mesh%nxg_
    	nyg=mesh%nyg_
    	nzg=mesh%nzg_
    	
    	Lx=mesh%Lxg_
    	Ly=mesh%Lyg_
    	Lz=mesh%Lzg_
    	
    	wrap_x = mesh%ptrMPIC_%wrapAround_(1)
    	wrap_y = mesh%ptrMPIC_%wrapAround_(2)
    	wrap_z = mesh%ptrMPIC_%wrapAround_(3)
    
    	lbi=lbound(vofb%c,1)
    	ubi=ubound(vofb%c,1)
    	lbj=lbound(vofb%c,2)
    	ubj=ubound(vofb%c,2)
    	lbk=lbound(vofb%c,3)
    	ubk=ubound(vofb%c,3)
		
		!pos x
		if (wrap_x) then
			do i=lbi,ubi
				r=i-nxg
				n=vofb%idx_mx(i)
				sw=min(abs(i-n),1)
				vofb%xc(i)=gmesh%xc_(n)+sw*sign(Lx,r)
				vofb%xf(i-1)=gmesh%xf_(n-1)+sw*sign(Lx,r)
			end do
			r=(i-1)-nxg
			vofb%xf(i-1)=gmesh%xf_(n)+sw*sign(Lx,r)
		else
			vofb%xc=gmesh%xc_(lbi:ubi)
			vofb%xf=gmesh%xf_(lbi-1:ubi)
		end if
		
		!pos y
		if (wrap_y) then
			do i=lbj,ubj
				r=i-nyg
				n=vofb%idx_my(i)
				sw=min(abs(i-n),1)
				vofb%yc(i)=gmesh%yc_(n)+sw*sign(Ly,r)
				vofb%yf(i-1)=gmesh%yf_(n-1)+sw*sign(Ly,r)
			end do
			r=(i-1)-nyg
			vofb%yf(i-1)=gmesh%yf_(n)+sw*sign(Ly,r)
		else
			vofb%yc=gmesh%yc_(lbj:ubj)
			vofb%yf=gmesh%yf_(lbj-1:ubj)
		end if
		
		!pos z
		if (wrap_z) then
			do i=lbk,ubk
				r=i-nzg
				n=vofb%idx_mz(i)
				sw=min(abs(i-n),1)
				vofb%zc(i)=gmesh%zc_(n)+sw*sign(Lz,r)
				vofb%zf(i-1)=gmesh%zf_(n-1)+sw*sign(Lz,r)
			end do
			r=(i-1)-nzg
			vofb%zf(i-1)=gmesh%zf_(n)+sw*sign(Lz,r)
		else
			vofb%zc=gmesh%zc_(lbk:ubk)
			vofb%zf=gmesh%zf_(lbk-1:ubk)
		end if
			
		
		!delta
		do i=lbi+1,ubi
			vofb%dxc(i) = vofb%xc(i)-vofb%xc(i-1)
		end do
		do i=lbj+1,ubj
			vofb%dyc(i) = vofb%yc(i)-vofb%yc(i-1)
		end do
		do i=lbk+1,ubk
			vofb%dzc(i) = vofb%zc(i)-vofb%zc(i-1)
		end do
		
		do i=lbi,ubi
			vofb%dxf(i) = vofb%xf(i)-vofb%xf(i-1)
		end do
		do i=lbj,ubj
			vofb%dyf(i) = vofb%yf(i)-vofb%yf(i-1)
		end do
		do i=lbk,ubk
			vofb%dzf(i) = vofb%zf(i)-vofb%zf(i-1)
		end do

    	
     end subroutine   	
!========================================================================================!

!========================================================================================!
    subroutine computeModuloIdx(mesh,vofb,is,ie,js,je,ks,ke)
    	type(grid), intent(in) :: mesh
    	type(vofBlock), intent(inout) :: vofb
    	integer, intent(in) :: is,ie,js,je,ks,ke
    	integer :: nxg,nyg,nzg
    	type(mpiControl), pointer :: mpic
    	logical :: wrap_x,wrap_y,wrap_z
		integer :: i
		
		mpic => mesh%ptrMPIC_
		wrap_x = mpic%wrapAround_(1)
		wrap_y = mpic%wrapAround_(2)
		wrap_z = mpic%wrapAround_(3)
		
		nxg = mesh%nxg_
		nyg = mesh%nyg_
		nzg = mesh%nzg_
		
		if (wrap_x) then
			do i=is,ie
				vofb%idx_mx(i) = modulo(i-1,nxg)+1
			end do
		else
			do i=is,ie
				vofb%idx_mx(i) = i
			end do			
		end if
		
		if (wrap_y) then
			do i=js,je
				vofb%idx_my(i) = modulo(i-1,nyg)+1
			end do
		else
			do i=js,je
				vofb%idx_my(i) = i
			end do
		end if
		
		if (wrap_z) then		
			do i=ks,ke
				vofb%idx_mz(i) = modulo(i-1,nzg)+1
			end do
		else
			do i=ks,ke
				vofb%idx_mz(i) = i
			end do			
		end if
    	
     end subroutine   	
!========================================================================================!
	
!========================================================================================!
    subroutine bubbles_distribution_space(mesh,b_proc_bool,rank,nb)
    	type(grid), intent(in) :: mesh
    	logical, allocatable, dimension(:), intent(out) :: b_proc_bool
    	integer, intent(in) :: rank,nb
    	integer :: n,b,i0g,i1g,j0g,j1g,k0g,k1g,nx,ny,nz,ic,jc,kc
    	logical :: present
    	
    	
    	if (nb==0) then
    		call mpiAbort('Zero bubbles initialised')
    	end if
    	
    	nx=mesh%nxg_
    	ny=mesh%nyg_
    	nz=mesh%nzg_
    	
    	call allocateArray(b_proc_bool,1,nb)
    	
		call globalIndexesFromAll(mesh,i0g,i1g,j0g,j1g,k0g,k1g,rank,'cl')
			
		do b=1,nb
			ic=int(0.5d0*(s_idx_init(1,b)+s_idx_init(2,b)))
			jc=int(0.5d0*(s_idx_init(3,b)+s_idx_init(4,b)))
			kc=int(0.5d0*(s_idx_init(5,b)+s_idx_init(6,b)))
			call isBubbleInThisCore(ic,jc,kc,i0g,i1g,j0g,j1g,k0g,k1g,present)
			b_proc_bool(b)=present
		end do
			
			
    end subroutine
!========================================================================================!

!========================================================================================!
    subroutine init_indexes_box(mesh,bi,x0,y0,z0,R)
    	type(grid), intent(in) :: mesh
    	integer, intent(in) :: bi
    	real(DP), intent(in) :: x0,y0,z0,R
    	integer :: nx,ny,nz,i,j,k,nmin_x,nmax_x,nmin_y,nmax_y,nmin_z,nmax_z
    	real(DP) :: p0,p1
    
		nx = mesh%nx_
		ny = mesh%ny_
		nz = mesh%nz_
		
		!range x
		p0=x0-R
		p1=x0+R
		do i=1,nx
			if ((p0>=mesh%xc_(i)).AND.(p0<=mesh%xc_(i+1))) then
				nmin_x=i
			end if
			if ((p1>=mesh%xc_(i)).AND.(p1<=mesh%xc_(i+1))) then
				nmax_x=i+1
			end if
		end do
		
		!range y
		p0=y0-R
		p1=y0+R
		do j=1,ny
			if ((p0>=mesh%yc_(j)).AND.(p0<=mesh%yc_(j+1))) then
				nmin_y=j
			end if
			if ((p1>=mesh%yc_(j)).AND.(p1<=mesh%yc_(j+1))) then
				nmax_y=j+1
			end if
		end do
		
		!range z
		p0=z0-R
		p1=z0+R
		do k=1,nz
			if ((p0>=mesh%zc_(k)).AND.(p0<=mesh%zc_(k+1))) then
				nmin_z=k
			end if
			if ((p1>=mesh%zc_(k)).AND.(p1<=mesh%zc_(k+1))) then
				nmax_z=k+1
			end if
		end do
		
		s_idx_init(1,bi)=nmin_x
		s_idx_init(2,bi)=nmax_x
		s_idx_init(3,bi)=nmin_y
		s_idx_init(4,bi)=nmax_y
		s_idx_init(5,bi)=nmin_z
		s_idx_init(6,bi)=nmax_z
			
    end subroutine
!========================================================================================!

!========================================================================================!
    subroutine readSingleBubble(mesh,pfile)
    	type(grid), intent(in) :: mesh
		type(parFile), intent(in) :: pfile
		real(DP) :: x0,y0,z0,R
		
		call readParameter(pfile,x0,'x0',bcast=.FALSE.)
		call readParameter(pfile,y0,'y0',bcast=.FALSE.)
		call readParameter(pfile,z0,'z0',bcast=.FALSE.)
		call readParameter(pfile,R,'R',bcast=.FALSE.)
					
		!set static bubble number
		s_nb = 1
					
		!allocate bubble indexes array
		call allocateArray(s_idx_init,1,6,1,1)
		call init_indexes_box(mesh,1,x0,y0,z0,R)
		
		!allocate bubble position
		call allocateArray(s_pos_init,1,4,1,1)
		s_pos_init(1,1)=x0
		s_pos_init(2,1)=y0
		s_pos_init(3,1)=z0
		s_pos_init(4,1)=R

			
    end subroutine
!========================================================================================!

!========================================================================================!
    subroutine readTwoBubbles(mesh,pfile)
    	type(grid), intent(in) :: mesh
		type(parFile), intent(in) :: pfile
		real(DP) :: x0_b1,y0_b1,z0_b1,R_b1
		real(DP) :: x0_b2,y0_b2,z0_b2,R_b2
		
		!b1
		call readParameter(pfile,x0_b1,'x0_b1',bcast=.FALSE.)
		call readParameter(pfile,y0_b1,'y0_b1',bcast=.FALSE.)
		call readParameter(pfile,z0_b1,'z0_b1',bcast=.FALSE.)
		call readParameter(pfile,R_b1,'R_b1',bcast=.FALSE.)
		!b2
		call readParameter(pfile,x0_b2,'x0_b2',bcast=.FALSE.)
		call readParameter(pfile,y0_b2,'y0_b2',bcast=.FALSE.)
		call readParameter(pfile,z0_b2,'z0_b2',bcast=.FALSE.)
		call readParameter(pfile,R_b2,'R_b2',bcast=.FALSE.)
					
		!set static bubble number
		s_nb = 2
					
		!allocate bubble indexes array
		call allocateArray(s_idx_init,1,6,1,s_nb)
		!b1
		call init_indexes_box(mesh,1,x0_b1,y0_b1,z0_b1,R_b1)
		!b2
		call init_indexes_box(mesh,2,x0_b2,y0_b2,z0_b2,R_b2)
		
		!allocate bubble position
		call allocateArray(s_pos_init,1,4,1,s_nb)
		!b1
		s_pos_init(1,1)=x0_b1
		s_pos_init(2,1)=y0_b1
		s_pos_init(3,1)=z0_b1
		s_pos_init(4,1)=R_b1
		!b2
		s_pos_init(1,2)=x0_b2
		s_pos_init(2,2)=y0_b2
		s_pos_init(3,2)=z0_b2
		s_pos_init(4,2)=R_b2

			
    end subroutine
!========================================================================================!

!========================================================================================!
    subroutine readBubblesArray(mesh,pfile)
    	type(grid), intent(in) :: mesh
		type(parFile), intent(in) :: pfile
		integer :: nbx,nby,nbz
		real(DP) :: Lx,Ly,Lz,R,dx,dy,dz,dxb,dyb,dzb
		real(DP) :: sx,sy,sz,x0,y0,z0,sx_eps,sy_eps,sz_eps
		real(DP) :: rnd(3)
		integer :: i,j,k,bi
		integer :: nx,ny,nz
		logical :: random_distr
		
		nx = mesh%nx_
		ny = mesh%ny_
		nz = mesh%nz_
		
		Lx = mesh%Lx_
		Ly = mesh%Ly_
		Lz = mesh%Lz_
		
		dx = Lx/nx
		dy = Ly/ny
		dz = Lz/nz
		
		call readParameter(pfile,nbx,'nbx',bcast=.FALSE.)
		call readParameter(pfile,nby,'nby',bcast=.FALSE.)
		call readParameter(pfile,nbz,'nbz',bcast=.FALSE.)
		call readParameter(pfile,R,'R',bcast=.FALSE.)
		call readParameter(pfile,random_distr,'random_distr',bcast=.FALSE.)

		!bubble displacement 
		dxb = Lx/nbx
		dyb = Ly/nby
		dzb = Lz/nbz

		!min spacing
		sx = 0.5d0*(dxb)-R
		sy = 0.5d0*(dyb)-R
		sz = 0.5d0*(dzb)-R
		
		!check if bubbles fit
		if ( sx < 1.1d0*dx ) then
			call mpiAbort('Too many bubbles in x direction ')
		else if ( sy < 1.1d0*dy ) then
			call mpiAbort('Too many bubbles in y direction ')
		else if ( sz < 1.1d0*dz ) then
			call mpiAbort('Too many bubbles in z direction ')
		end if
		
		!set static bubble number
		s_nb = nbx*nby*nbz
		
		!allocate bubble indexes array
		call allocateArray(s_idx_init,1,6,1,s_nb)
		
		!allocate bubble potion
		call allocateArray(s_pos_init,1,4,1,s_nb)
		
		bi = 1
		do k=1,nbz
			do j=1,nby
				do i=1,nbx
				
					x0 = 0.5d0*dxb + (i-1)*dxb
					y0 = 0.5d0*dyb + (j-1)*dyb
					z0 = 0.5d0*dzb + (k-1)*dzb
					
					if (random_distr) then
						call RANDOM_NUMBER(rnd)
						x0 = x0 + (rnd(1)-0.5d0)*0.3d0*sx
						y0 = y0 + (rnd(2)-0.5d0)*0.3d0*sy
						z0 = z0 + (rnd(3)-0.5d0)*0.3d0*sz
					end if
								
					call init_indexes_box(mesh,bi,x0,y0,z0,R)
					
					s_pos_init(1,bi)=x0
					s_pos_init(2,bi)=y0
					s_pos_init(3,bi)=z0
					s_pos_init(4,bi)=R
					
					bi=bi+1
					
				end do
			end do
		end do

			
    end subroutine
!========================================================================================!

!========================================================================================!
    subroutine grid_2_boxes_u(mesh,u)
    	type(grid), intent(in) :: mesh
    	type(vfield), intent(in) :: u
    	type(mpiControl), pointer :: mpic
		integer :: nprocs,b,master,slave,n,bl
		integer :: i0g,j0g,k0g,i1g,j1g,k1g,nxg,nyg,nzg
    	integer :: is,ie,js,je,ks,ke
		integer :: tag,ierror
    	integer, dimension(3) :: requests
    	integer, dimension(MPI_STATUS_SIZE,3) :: status
    	integer, allocatable, dimension(:,:) :: ux_idx,uy_idx,uz_idx
    	real(DP), allocatable, dimension(:) :: ux_blk,uy_blk,uz_blk
    	logical :: wrap_x,wrap_y,wrap_z
    	character(len=:), allocatable :: buff_ux,buff_uy,buff_uz
    	integer :: sizeBuff_ux,sizeBuff_uy,sizeBuff_uz,position
    	

		mpic => mesh%ptrMPIC_
		nprocs = mpic%nProcs_
		
		nxg=mesh%nxg_
		nyg=mesh%nyg_
		nzg=mesh%nzg_
		
		wrap_x=mpic%wrapAround_(1)
		wrap_y=mpic%wrapAround_(2)
		wrap_z=mpic%wrapAround_(3)
    	
    	!send velocity field
    	do b=1,s_nb
    	  	
    		do n=0,nprocs-1
		
    			!manage own bubbls
    			if ((mpic%rank_==s_gbList(1,n,b)).AND.(s_gbList(1,n,b)==s_gbList(2,n,b))) then
    			
    				do bl=1,s_nblk
    					
    					if (vofBlocks(bl)%bn==b) then
    				
    						i0g = mesh%i0g_
    						j0g = mesh%j0g_
    						k0g = mesh%k0g_
    						i1g = mesh%i1g_
    						j1g = mesh%j1g_
    						k1g = mesh%k1g_
    						
    						is = s_gbList(3,n,b)-offset_u
    						ie = s_gbList(4,n,b)+offset_u
    						js = s_gbList(5,n,b)-offset_u
    						je = s_gbList(6,n,b)+offset_u
    						ks = s_gbList(7,n,b)-offset_u
    						ke = s_gbList(8,n,b)+offset_u 
    						
    						
							call allocateBlockSendRecvBuffers(ux_idx,ux_blk,1,is,ie,js,je,ks,ke)
    						call packSendToMasterBuff(ux_idx,ux_blk,u%ux_%f_,1,is,ie,js,je,ks,ke,	&
    									     	      offset_u,i0g,i1g,j0g,j1g,k0g,k1g,nxg,nyg,nzg,	&
    									      	      wrap_x,wrap_y,wrap_z)
							
							call allocateBlockSendRecvBuffers(uy_idx,uy_blk,2,is,ie,js,je,ks,ke)			      
    						call packSendToMasterBuff(uy_idx,uy_blk,u%uy_%f_,2,is,ie,js,je,ks,ke,	&
    									      		  offset_u,i0g,i1g,j0g,j1g,k0g,k1g,nxg,nyg,nzg,	&
    									      		  wrap_x,wrap_y,wrap_z)

							call allocateBlockSendRecvBuffers(uz_idx,uz_blk,3,is,ie,js,je,ks,ke)    									
    						call packSendToMasterBuff(uz_idx,uz_blk,u%uz_%f_,3,is,ie,js,je,ks,ke,	&
    									      		  offset_u,i0g,i1g,j0g,j1g,k0g,k1g,nxg,nyg,nzg,	&
    									      	      wrap_x,wrap_y,wrap_z)
    									      		      
    						call unPackRecvFromSlaveBuff(vofBlocks(bl)%ux,ux_idx,ux_blk)
    						call unPackRecvFromSlaveBuff(vofBlocks(bl)%uy,uy_idx,uy_blk)
    						call unPackRecvFromSlaveBuff(vofBlocks(bl)%uz,uz_idx,uz_blk)
  
    						
    					end if 
    				
    				end do	
    				
    			
    			!send to master
    			else if (mpic%rank_==s_gbList(2,n,b)) then
    			
					master = s_gbList(1,n,b)
					
    				i0g = mesh%i0g_
    				j0g = mesh%j0g_
    				k0g = mesh%k0g_
    				i1g = mesh%i1g_
    				j1g = mesh%j1g_
    				k1g = mesh%k1g_
    			
    				is = s_gbList(3,n,b)-offset_u
    				ie = s_gbList(4,n,b)+offset_u
    				js = s_gbList(5,n,b)-offset_u
    				je = s_gbList(6,n,b)+offset_u
    				ks = s_gbList(7,n,b)-offset_u
    				ke = s_gbList(8,n,b)+offset_u 
    				
    				!pack ux
					call allocateBlockSendRecvBuffers(ux_idx,ux_blk,1,is,ie,js,je,ks,ke)		
    				call packSendToMasterBuff(ux_idx,ux_blk,u%ux_%f_,1,is,ie,js,je,ks,ke,	&
    							       	 	  offset_u,i0g,i1g,j0g,j1g,k0g,k1g,nxg,nyg,nzg,	&
    							      		  wrap_x,wrap_y,wrap_z)
    				sizeBuff_ux = integer_size*size(ux_idx) + realDP_size*size(ux_blk)	
    				call reAllocateArray(buff_ux,sizeBuff_ux)      
					position = 0
					call MPI_PACK(ux_idx, size(ux_idx), MPI_INTEGER, buff_ux, sizeBuff_ux, position, &
							      mpic%cartComm_, ierror)
					call MPI_PACK(ux_blk, size(ux_blk), MPI_DOUBLE_PRECISION, buff_ux, sizeBuff_ux, &
								  position, mpic%cartComm_, ierror)
    				
    				!pack uy
					call allocateBlockSendRecvBuffers(uy_idx,uy_blk,2,is,ie,js,je,ks,ke)						      
    				call packSendToMasterBuff(uy_idx,uy_blk,u%uy_%f_,2,is,ie,js,je,ks,ke,	&
    							      	 	  offset_u,i0g,i1g,j0g,j1g,k0g,k1g,nxg,nyg,nzg,	&
    							      	 	  wrap_x,wrap_y,wrap_z)
    				sizeBuff_uy = integer_size*size(uy_idx) + realDP_size*size(uy_blk)	
    				call reAllocateArray(buff_uy,sizeBuff_uy)      
					position = 0
					call MPI_PACK(uy_idx, size(uy_idx), MPI_INTEGER, buff_uy, sizeBuff_uy, position, &
							      mpic%cartComm_, ierror)
					call MPI_PACK(uy_blk, size(uy_blk), MPI_DOUBLE_PRECISION, buff_uy, sizeBuff_uy, &
								  position, mpic%cartComm_, ierror)
    				
    				!pack uz	
					call allocateBlockSendRecvBuffers(uz_idx,uz_blk,3,is,ie,js,je,ks,ke)				
    				call packSendToMasterBuff(uz_idx,uz_blk,u%uz_%f_,3,is,ie,js,je,ks,ke,	&
    							      	 	  offset_u,i0g,i1g,j0g,j1g,k0g,k1g,nxg,nyg,nzg,	&
    							      	 	  wrap_x,wrap_y,wrap_z)
    				sizeBuff_uz = integer_size*size(uz_idx) + realDP_size*size(uz_blk)	
    				call reAllocateArray(buff_uz,sizeBuff_uz)      
					position = 0
					call MPI_PACK(uz_idx, size(uz_idx), MPI_INTEGER, buff_uz, sizeBuff_uz, position, &
							      mpic%cartComm_, ierror)
					call MPI_PACK(uz_blk, size(uz_blk), MPI_DOUBLE_PRECISION, buff_uz, sizeBuff_uz, &
								  position, mpic%cartComm_, ierror)
    							      

					tag = 0	
					call MPI_ISEND(buff_ux, sizeBuff_ux, MPI_CHARACTER, master, &
						       	   tag, mpic%cartComm_, requests(1), ierror) 
					tag = 1
					call MPI_ISEND(buff_uy, sizeBuff_uy, MPI_CHARACTER, master, &
						       	   tag, mpic%cartComm_, requests(2), ierror) 
					tag = 2
					call MPI_ISEND(buff_uz, sizeBuff_uz, MPI_CHARACTER, master, &
						       	   tag, mpic%cartComm_, requests(3), ierror) 
	           
					call MPI_WAITALL(3, requests, status, ierror)  	

    			!recv from slaves
    			else if (mpic%rank_==s_gbList(1,n,b)) then
    			
    				slave = s_gbList(2,n,b)
    				
    				call globalIndexesFromAll(mesh,i0g,i1g,j0g,j1g,k0g,k1g,slave,'cl')
    				
    				is = s_gbList(3,n,b)-offset_u
    				ie = s_gbList(4,n,b)+offset_u
    				js = s_gbList(5,n,b)-offset_u
    				je = s_gbList(6,n,b)+offset_u
    				ks = s_gbList(7,n,b)-offset_u
    				ke = s_gbList(8,n,b)+offset_u 
	
  
					call allocateBlockSendRecvBuffers(ux_idx,ux_blk,1,is,ie,js,je,ks,ke)
					call allocateBlockSendRecvBuffers(uy_idx,uy_blk,2,is,ie,js,je,ks,ke)
					call allocateBlockSendRecvBuffers(uz_idx,uz_blk,3,is,ie,js,je,ks,ke)
    				
 					tag = 0
					sizeBuff_ux = integer_size*size(ux_idx) + realDP_size*size(ux_blk)
					call reAllocateArray(buff_ux,sizeBuff_ux) 
					call MPI_IRECV(buff_ux,sizeBuff_ux,MPI_CHARACTER, slave, &
					               tag, mpic%cartComm_, requests(1), ierror)	
					               
					tag = 1
					sizeBuff_uy = integer_size*size(uy_idx) + realDP_size*size(uy_blk)
					call reAllocateArray(buff_uy,sizeBuff_uy)
					call MPI_IRECV(buff_uy,sizeBuff_uy,MPI_CHARACTER, slave, &
					               tag, mpic%cartComm_, requests(2), ierror)	
					               
					tag = 2
					sizeBuff_uz = integer_size*size(uz_idx) + realDP_size*size(uz_blk)
					call reAllocateArray(buff_uz,sizeBuff_uz)
					call MPI_IRECV(buff_uz,sizeBuff_uz,MPI_CHARACTER, slave, &
					               tag, mpic%cartComm_, requests(3), ierror)
  

					call MPI_WAITALL(3, requests, status, ierror) 
					
					!unpack ux
					position = 0
					call MPI_UNPACK(buff_ux, sizeBuff_ux, position, ux_idx, size(ux_idx), &
								    MPI_INTEGER, mpic%cartComm_, ierror)
					call MPI_UNPACK(buff_ux, sizeBuff_ux, position, ux_blk, size(ux_blk), &
								    MPI_DOUBLE_PRECISION, mpic%cartComm_, ierror)	
								    
					!unpack uy
					position = 0
					call MPI_UNPACK(buff_uy, sizeBuff_uy, position, uy_idx, size(uy_idx), &
								    MPI_INTEGER, mpic%cartComm_, ierror)
					call MPI_UNPACK(buff_uy, sizeBuff_uy, position, uy_blk, size(uy_blk), &
								    MPI_DOUBLE_PRECISION, mpic%cartComm_, ierror)
								    
					!unpack uz
					position = 0
					call MPI_UNPACK(buff_uz, sizeBuff_uz, position, uz_idx, size(uz_idx), &
								    MPI_INTEGER, mpic%cartComm_, ierror)
					call MPI_UNPACK(buff_uz, sizeBuff_uz, position, uz_blk, size(uz_blk), &
								    MPI_DOUBLE_PRECISION, mpic%cartComm_, ierror)				 
  
    				do bl=1,s_nblk
    					if (vofBlocks(bl)%bn==b) then
    					
    						call unPackRecvFromSlaveBuff(vofBlocks(bl)%ux,ux_idx,ux_blk)
    						call unPackRecvFromSlaveBuff(vofBlocks(bl)%uy,uy_idx,uy_blk)
    						call unPackRecvFromSlaveBuff(vofBlocks(bl)%uz,uz_idx,uz_blk)   						
    					
    					end if
    				end do
			
    				
    			end if
    		
    		end do
    	end do


    end subroutine   	
!========================================================================================!

!========================================================================================!
    subroutine boxes_2_grid_f(mesh,c,field_2_pack,unpack_op)
    	type(grid), intent(in) :: mesh
    	type(field), intent(inout) :: c
    	integer, intent(in) :: field_2_pack,unpack_op
    	type(mpiControl), pointer :: mpic
		integer :: nprocs,b,master,slave,n,bl
		integer :: i0g,j0g,k0g,i1g,j1g,k1g
		integer :: is,js,ks,ie,je,ke,nxg,nyg,nzg
		integer :: tag,ierror,position,sizeBuff,offset
    	integer, dimension(MPI_STATUS_SIZE) :: status
    	integer, allocatable, dimension(:,:) :: c_idx
    	real(DP), allocatable, dimension(:) :: c_blk
    	character(len=:), allocatable :: buff
    	logical :: wrap_x,wrap_y,wrap_z
    	
    	
		select case(field_2_pack)
			case(PACK_BOX_VF)
				offset=0
			case(PACK_BOX_K)
				offset=0
			case default
		end select

		mpic => mesh%ptrMPIC_
		nprocs = mpic%nProcs_
		tag=0
		
		nxg=mesh%nxg_
		nyg=mesh%nyg_
		nzg=mesh%nzg_
		
		wrap_x=mpic%wrapAround_(1)
		wrap_y=mpic%wrapAround_(2)
		wrap_z=mpic%wrapAround_(3)
		
		!reset field
		call set2zero_omp(c%f_)
    	
    	do b=1,s_nb
    	  	
    		do n=0,nprocs-1
		
    			!manage own bubbls
    			if ((mpic%rank_==s_gbList(1,n,b)).AND.(s_gbList(1,n,b)==s_gbList(2,n,b))) then
    			
    				i0g = mesh%i0g_
    				i1g = mesh%i1g_
    				j0g = mesh%j0g_
    				j1g = mesh%j1g_
    				k0g = mesh%k0g_
    				k1g = mesh%k1g_
    			
    				do bl=1,s_nblk
    					
    					if (vofBlocks(bl)%bn==b) then

    						is = vofBlocks(bl)%idx(1)-offset
    						ie = vofBlocks(bl)%idx(2)+offset
    						js = vofBlocks(bl)%idx(3)-offset
    						je = vofBlocks(bl)%idx(4)+offset
    						ks = vofBlocks(bl)%idx(5)-offset
    						ke = vofBlocks(bl)%idx(6)+offset
    					
    						call allocateBlockSendRecvBuffers(c_idx,c_blk,0,is,ie,js,je,ks,ke)
							call packSendToSlaveBuff(c_idx,c_blk,vofBlocks(bl),&
													 field_2_pack,0,i0g,i1g,j0g,j1g,k0g,k1g,&
													 is,ie,js,je,ks,ke,nxg,nyg,nzg,&
    								  				 wrap_x,wrap_y,wrap_z)
													 
							call unPackRecvFromMasterBuff(c_idx,c_blk,c%f_,unpack_op)
						
						end if
						
					end do
						
    			
    			!recv from master
    			else if (mpic%rank_==s_gbList(2,n,b)) then
    			
					master = s_gbList(1,n,b)

 					is = s_gbList(3,n,b)-offset
					ie = s_gbList(4,n,b)+offset
					js = s_gbList(5,n,b)-offset
					je = s_gbList(6,n,b)+offset
					ks = s_gbList(7,n,b)-offset
					ke = s_gbList(8,n,b)+offset 						
    				
    				call allocateBlockSendRecvBuffers(c_idx,c_blk,0,is,ie,js,je,ks,ke)
    				
					sizeBuff = integer_size*size(c_idx) + realDP_size*size(c_blk)
					call reAllocateArray(buff,sizeBuff)
					call MPI_RECV(buff,sizeBuff,MPI_CHARACTER, master, &
					               tag, mpic%cartComm_, status, ierror)
					               
					position = 0
					call MPI_UNPACK(buff, sizeBuff, position, c_idx, size(c_idx), &
								    MPI_INTEGER, mpic%cartComm_, ierror)
					call MPI_UNPACK(buff, sizeBuff, position, c_blk, size(c_blk), &
								    MPI_DOUBLE_PRECISION, mpic%cartComm_, ierror)					
					               
					call unPackRecvFromMasterBuff(c_idx,c_blk,c%f_,unpack_op)
				 			
    	
    			!send to slaves
    			else if (mpic%rank_==s_gbList(1,n,b)) then
    			
    				slave = s_gbList(2,n,b)

    				call globalIndexesFromAll(mesh,i0g,i1g,j0g,j1g,k0g,k1g,slave,'cl')
    			
    				is = s_gbList(3,n,b)-offset
    				ie = s_gbList(4,n,b)+offset
    				js = s_gbList(5,n,b)-offset
    				je = s_gbList(6,n,b)+offset
    				ks = s_gbList(7,n,b)-offset
    				ke = s_gbList(8,n,b)+offset
    				
    				do bl=1,s_nblk
    					if (vofBlocks(bl)%bn==b) then
							
							call allocateBlockSendRecvBuffers(c_idx,c_blk,0,is,ie,js,je,ks,ke)
							call packSendToSlaveBuff(c_idx,c_blk,vofBlocks(bl),&
													 field_2_pack,0,i0g,i1g,j0g,j1g,k0g,k1g,&
													 is,ie,js,je,ks,ke,nxg,nyg,nzg,&
    								  				 wrap_x,wrap_y,wrap_z)
							
    					end if
    				end do
    				
    				
    				sizeBuff = integer_size*size(c_idx) + realDP_size*size(c_blk)	
    				call reAllocateArray(buff,sizeBuff) 
    				
					position = 0
					call MPI_PACK(c_idx, size(c_idx), MPI_INTEGER, buff, sizeBuff, position, &
							      mpic%cartComm_, ierror)
					call MPI_PACK(c_blk, size(c_blk), MPI_DOUBLE_PRECISION, buff, sizeBuff, &
								  position, mpic%cartComm_, ierror)
   				
					call MPI_SEND(buff, sizeBuff, MPI_CHARACTER, slave, tag, mpic%cartComm_, ierror)     				
    											
    				
    			end if
    		
    		end do
    	end do
    	
    	call updateBoundaries(c)


    end subroutine   	
!========================================================================================!

!========================================================================================!
    subroutine boxes_2_grid_vf(mesh,vf,field_2_pack,unpack_op)
    	type(grid), intent(in) :: mesh
    	type(vfield), intent(inout) :: vf
    	integer, intent(in) :: field_2_pack,unpack_op
    	type(mpiControl), pointer :: mpic
		integer :: nprocs,b,master,slave,n,bl,offset
		integer :: i0g,j0g,k0g,i1g,j1g,k1g
		integer :: is,js,ks,ie,je,ke,nxg,nyg,nzg
		integer :: tag,ierror,position
		integer :: sizeBuff_vfx,sizeBuff_vfy,sizeBuff_vfz
    	integer, dimension(3) :: requests
    	integer, dimension(MPI_STATUS_SIZE,3) :: status
    	integer, allocatable, dimension(:,:) :: vfx_idx,vfy_idx,vfz_idx
    	real(DP), allocatable, dimension(:) :: vfx_blk,vfy_blk,vfz_blk
    	character(len=:), allocatable :: buff_vfx,buff_vfy,buff_vfz
    	logical :: wrap_x,wrap_y,wrap_z
    	

		select case(field_2_pack)
			case(PACK_BOX_U)
				offset=offset_u
			case(PACK_BOX_ST)
				offset=offset_st
			case default
		end select
    	

		mpic => mesh%ptrMPIC_
		nprocs = mpic%nProcs_
		
		!reset field
		call set2zero_omp(vf%ux_%f_)
		call set2zero_omp(vf%uy_%f_)
		call set2zero_omp(vf%uz_%f_)
		
		nxg=mesh%nxg_
		nyg=mesh%nyg_
		nzg=mesh%nzg_
		
		wrap_x=mpic%wrapAround_(1)
		wrap_y=mpic%wrapAround_(2)
		wrap_z=mpic%wrapAround_(3)
    	
    	do b=1,s_nb
    	  	
    		do n=0,nprocs-1
		
    			!manage own bubbls
    			if ((mpic%rank_==s_gbList(1,n,b)).AND.(s_gbList(1,n,b)==s_gbList(2,n,b))) then
    			
    				i0g = mesh%i0g_
    				i1g = mesh%i1g_
    				j0g = mesh%j0g_
    				j1g = mesh%j1g_
    				k0g = mesh%k0g_
    				k1g = mesh%k1g_
    			
    				do bl=1,s_nblk
    					
    					if (vofBlocks(bl)%bn==b) then
    					
    						is = vofBlocks(bl)%idx(1)-offset
    						ie = vofBlocks(bl)%idx(2)+offset
    						js = vofBlocks(bl)%idx(3)-offset
    						je = vofBlocks(bl)%idx(4)+offset
    						ks = vofBlocks(bl)%idx(5)-offset
    						ke = vofBlocks(bl)%idx(6)+offset
    					
    						
    						call allocateBlockSendRecvBuffers(vfx_idx,vfx_blk,1,is,ie,js,je,ks,ke)
    						call packSendToSlaveBuff(vfx_idx,vfx_blk,vofBlocks(bl),&
    										         field_2_pack,1,i0g,i1g,j0g,j1g,k0g,k1g,&
    								  		         is,ie,js,je,ks,ke,nxg,nyg,nzg,&
    								  				 wrap_x,wrap_y,wrap_z)
											            
    						call allocateBlockSendRecvBuffers(vfy_idx,vfy_blk,2,is,ie,js,je,ks,ke)
    						call packSendToSlaveBuff(vfy_idx,vfy_blk,vofBlocks(bl),&
    										         field_2_pack,2,i0g,i1g,j0g,j1g,k0g,k1g,&
    								  				 is,ie,js,je,ks,ke,nxg,nyg,nzg,&
    								  				 wrap_x,wrap_y,wrap_z)
											            
    						call allocateBlockSendRecvBuffers(vfz_idx,vfz_blk,3,is,ie,js,je,ks,ke)
    						call packSendToSlaveBuff(vfz_idx,vfz_blk,vofBlocks(bl),&
    										         field_2_pack,3,i0g,i1g,j0g,j1g,k0g,k1g,&
    								  				 is,ie,js,je,ks,ke,nxg,nyg,nzg,&
    								  				 wrap_x,wrap_y,wrap_z)
							
							call unPackRecvFromMasterBuff(vfx_idx,vfx_blk,vf%ux_%f_,unpack_op)
							call unPackRecvFromMasterBuff(vfy_idx,vfy_blk,vf%uy_%f_,unpack_op)
							call unPackRecvFromMasterBuff(vfz_idx,vfz_blk,vf%uz_%f_,unpack_op)
							
						end if
						
					end do
						
    			
    			!recv from master
    			else if (mpic%rank_==s_gbList(2,n,b)) then
    			
					master = s_gbList(1,n,b)

 					is = s_gbList(3,n,b)-offset
					ie = s_gbList(4,n,b)+offset
					js = s_gbList(5,n,b)-offset
					je = s_gbList(6,n,b)+offset
					ks = s_gbList(7,n,b)-offset
					ke = s_gbList(8,n,b)+offset 						
    				
    				call allocateBlockSendRecvBuffers(vfx_idx,vfx_blk,1,is,ie,js,je,ks,ke)
    				call allocateBlockSendRecvBuffers(vfy_idx,vfy_blk,2,is,ie,js,je,ks,ke)
    				call allocateBlockSendRecvBuffers(vfz_idx,vfz_blk,3,is,ie,js,je,ks,ke)
    				
    				tag=0
					sizeBuff_vfx = integer_size*size(vfx_idx) + realDP_size*size(vfx_blk)
					call reAllocateArray(buff_vfx,sizeBuff_vfx)
					call MPI_IRECV(buff_vfx,sizeBuff_vfx,MPI_CHARACTER, master, &
					               tag, mpic%cartComm_, requests(1), ierror)
					               
    				tag=1
					sizeBuff_vfy = integer_size*size(vfy_idx) + realDP_size*size(vfy_blk)
					call reAllocateArray(buff_vfy,sizeBuff_vfy)
					call MPI_IRECV(buff_vfy,sizeBuff_vfy,MPI_CHARACTER, master, &
					               tag, mpic%cartComm_, requests(2), ierror)
					               
    				tag=2
					sizeBuff_vfz = integer_size*size(vfz_idx) + realDP_size*size(vfz_blk)
					call reAllocateArray(buff_vfz,sizeBuff_vfz)
					call MPI_IRECV(buff_vfz,sizeBuff_vfz,MPI_CHARACTER, master, &
					               tag, mpic%cartComm_, requests(3), ierror)
					               
					call MPI_WAITALL(3, requests, status, ierror) 
					
					!unpack vfx               
					position = 0
					call MPI_UNPACK(buff_vfx, sizeBuff_vfx, position, vfx_idx, size(vfx_idx), &
								    MPI_INTEGER, mpic%cartComm_, ierror)
					call MPI_UNPACK(buff_vfx, sizeBuff_vfx, position, vfx_blk, size(vfx_blk), &
								    MPI_DOUBLE_PRECISION, mpic%cartComm_, ierror)					

					!unpack vfy               
					position = 0
					call MPI_UNPACK(buff_vfy, sizeBuff_vfy, position, vfy_idx, size(vfy_idx), &
								    MPI_INTEGER, mpic%cartComm_, ierror)
					call MPI_UNPACK(buff_vfy, sizeBuff_vfy, position, vfy_blk, size(vfy_blk), &
								    MPI_DOUBLE_PRECISION, mpic%cartComm_, ierror)	
								    
					!unpack vfz               
					position = 0
					call MPI_UNPACK(buff_vfz, sizeBuff_vfz, position, vfz_idx, size(vfz_idx), &
								    MPI_INTEGER, mpic%cartComm_, ierror)
					call MPI_UNPACK(buff_vfz, sizeBuff_vfz, position, vfz_blk, size(vfz_blk), &
								    MPI_DOUBLE_PRECISION, mpic%cartComm_, ierror)				               
					
					call unPackRecvFromMasterBuff(vfx_idx,vfx_blk,vf%ux_%f_,unpack_op)
					call unPackRecvFromMasterBuff(vfy_idx,vfy_blk,vf%uy_%f_,unpack_op)
					call unPackRecvFromMasterBuff(vfz_idx,vfz_blk,vf%uz_%f_,unpack_op)
					

    			!send to slaves
    			else if (mpic%rank_==s_gbList(1,n,b)) then
    			
    				slave = s_gbList(2,n,b)

    				call globalIndexesFromAll(mesh,i0g,i1g,j0g,j1g,k0g,k1g,slave,'cl')
    			
    				is = s_gbList(3,n,b)-offset
    				ie = s_gbList(4,n,b)+offset
    				js = s_gbList(5,n,b)-offset
    				je = s_gbList(6,n,b)+offset
    				ks = s_gbList(7,n,b)-offset
    				ke = s_gbList(8,n,b)+offset
    				
    				do bl=1,s_nblk
    					if (vofBlocks(bl)%bn==b) then
							
							call allocateBlockSendRecvBuffers(vfx_idx,vfx_blk,1,is,ie,js,je,ks,ke)
							call allocateBlockSendRecvBuffers(vfy_idx,vfy_blk,2,is,ie,js,je,ks,ke)
							call allocateBlockSendRecvBuffers(vfz_idx,vfz_blk,3,is,ie,js,je,ks,ke)
							
    						call packSendToSlaveBuff(vfx_idx,vfx_blk,vofBlocks(bl),&
    												 field_2_pack,1,i0g,i1g,j0g,j1g,k0g,k1g,&
    								  				 is,ie,js,je,ks,ke,nxg,nyg,nzg,&
    								  				 wrap_x,wrap_y,wrap_z)
    						call packSendToSlaveBuff(vfy_idx,vfy_blk,vofBlocks(bl),&
    											     field_2_pack,2,i0g,i1g,j0g,j1g,k0g,k1g,&
    								  				 is,ie,js,je,ks,ke,nxg,nyg,nzg,&
    								  				 wrap_x,wrap_y,wrap_z)
    						call packSendToSlaveBuff(vfz_idx,vfz_blk,vofBlocks(bl),&
    											     field_2_pack,3,i0g,i1g,j0g,j1g,k0g,k1g,&
    								  				 is,ie,js,je,ks,ke,nxg,nyg,nzg,&
    								  				 wrap_x,wrap_y,wrap_z)
							
    					end if
    				end do
    				
    				!pack vfx
    				sizeBuff_vfx = integer_size*size(vfx_idx) + realDP_size*size(vfx_blk)	
    				call reAllocateArray(buff_vfx,sizeBuff_vfx) 
					position = 0
					call MPI_PACK(vfx_idx, size(vfx_idx), MPI_INTEGER, buff_vfx, sizeBuff_vfx, position, &
							      mpic%cartComm_, ierror)
					call MPI_PACK(vfx_blk, size(vfx_blk), MPI_DOUBLE_PRECISION, buff_vfx, sizeBuff_vfx, &
								  position, mpic%cartComm_, ierror)
					
					!pack vfy		  
    				sizeBuff_vfy = integer_size*size(vfy_idx) + realDP_size*size(vfy_blk)	
    				call reAllocateArray(buff_vfy,sizeBuff_vfy) 
					position = 0
					call MPI_PACK(vfy_idx, size(vfy_idx), MPI_INTEGER, buff_vfy, sizeBuff_vfy, position, &
							      mpic%cartComm_, ierror)
					call MPI_PACK(vfy_blk, size(vfy_blk), MPI_DOUBLE_PRECISION, buff_vfy, sizeBuff_vfy, &
								  position, mpic%cartComm_, ierror)
					
					!pack vfz			  
    				sizeBuff_vfz = integer_size*size(vfz_idx) + realDP_size*size(vfz_blk)	
    				call reAllocateArray(buff_vfz,sizeBuff_vfz) 
					position = 0
					call MPI_PACK(vfz_idx, size(vfz_idx), MPI_INTEGER, buff_vfz, sizeBuff_vfz, position, &
							      mpic%cartComm_, ierror)
					call MPI_PACK(vfz_blk, size(vfz_blk), MPI_DOUBLE_PRECISION, buff_vfz, sizeBuff_vfz, &
								  position, mpic%cartComm_, ierror)
		  
					tag = 0	
					call MPI_ISEND(buff_vfx, sizeBuff_vfx, MPI_CHARACTER, slave, &
						       	   tag, mpic%cartComm_, requests(1), ierror) 
					tag = 1
					call MPI_ISEND(buff_vfy, sizeBuff_vfy, MPI_CHARACTER, slave, &
						       	   tag, mpic%cartComm_, requests(2), ierror) 
					tag = 2
					call MPI_ISEND(buff_vfz, sizeBuff_vfz, MPI_CHARACTER, slave, &
						       	   tag, mpic%cartComm_, requests(3), ierror) 
	           
					call MPI_WAITALL(3, requests, status, ierror)		    				   											
    				
    			end if
    		
    		end do
    	end do
    	
    	call updateBoundariesV(vf)


    end subroutine   	
!========================================================================================!

!========================================================================================!
    subroutine allocateBlockSendRecvBuffers(idx_buff,val_buff,ftype,is,ie,js,je,ks,ke)
    	integer, allocatable, dimension(:,:), intent(inout) :: idx_buff
    	real(DP), allocatable, dimension(:), intent(inout) :: val_buff
    	integer, intent(in) :: is,ie,js,je,ks,ke,ftype
    	integer :: sizeBuff
    	
		select case(ftype)
			case(1)
				sizeBuff=(ie-is+2)*(je-js+1)*(ke-ks+1)
			case(2)
				sizeBuff=(ie-is+1)*(je-js+2)*(ke-ks+1)
			case(3)
				sizeBuff=(ie-is+1)*(je-js+1)*(ke-ks+2)
			case default
				sizeBuff=(ie-is+1)*(je-js+1)*(ke-ks+1)
		end select
		
		call reAllocateArray(idx_buff,1,4,1,sizeBuff)
		call reAllocateArray(val_buff,1,sizeBuff) 
		val_buff=0.d0	
    	    		
	end subroutine
!========================================================================================!

!========================================================================================!
    subroutine packSendToSlaveBuff(idx_buff,val_buff,vofb,field_2_pack,ftype,i0g,i1g,&
    						       j0g,j1g,k0g,k1g,is_org,ie_org,js_org,je_org,ks_org,&
    						       ke_org,nxg,nyg,nzg,wrap_x,wrap_y,wrap_z)
    	integer, allocatable, dimension(:,:), intent(inout) :: idx_buff
    	real(DP), allocatable, dimension(:), intent(inout) :: val_buff
    	type(vofBlock), intent(in), target :: vofb
    	integer, intent(in) :: field_2_pack,ftype,i0g,i1g,j0g,j1g,k0g,k1g,nxg,nyg,nzg
    	integer, intent(in) :: is_org,ie_org,js_org,je_org,ks_org,ke_org
    	integer :: i,j,k,i0gm,j0gm,k0gm,i1gm,j1gm,k1gm,il,jl,kl,n,is,js,ks,ie,je,ke
    	integer, allocatable, dimension(:) :: im,jm,km
    	logical, intent(in) :: wrap_x,wrap_y,wrap_z
    	real(DP), pointer, dimension(:,:,:) :: ptrf => NULL()

    	select case(field_2_pack)
    	
    		case(PACK_BOX_ST)
    			select case(ftype)
    				case(1)
    					ptrf => vofb%stx
    				case(2)
    					ptrf => vofb%sty
    				case(3)
    					ptrf => vofb%stz
    				case default
    			end select

    		case(PACK_BOX_U)
    			select case(ftype)
    				case(1)
    					ptrf => vofb%ux
    				case(2)
    					ptrf => vofb%uy
    				case(3)
    					ptrf => vofb%uz
    				case default
    			end select
    			
    		case(PACK_BOX_VF)
    			ptrf => vofb%c
    			
    		case(PACK_BOX_K)
    			ptrf => vofb%k
    			
			case default
			
		end select
    
    
    	!copy ranges
    	is=is_org
    	js=js_org
    	ks=ks_org
    	ie=ie_org
    	je=je_org
    	ke=ke_org

    	select case(ftype)
    		case(0)
    		    call allocateArray(im,is,ie)
				call allocateArray(jm,js,je)
				call allocateArray(km,ks,ke)
    		case(1)
    		    call allocateArray(im,is-1,ie)
				call allocateArray(jm,js,je)
				call allocateArray(km,ks,ke)
    		case(2)
    		    call allocateArray(im,is,ie)
				call allocateArray(jm,js-1,je)
				call allocateArray(km,ks,ke)
    		case(3)
    		    call allocateArray(im,is,ie)
				call allocateArray(jm,js,je)
				call allocateArray(km,ks-1,ke)
    		case default
    	end select
    

    	!exception periodic bc
    	if (wrap_x) then
    		do i=is,ie
			im(i)=modulo(i-1,nxg)+1
		end do
    	else
    		do i=is,ie
			im(i)=i
		end do	
    	end if
    	if (wrap_y) then
    		do j=js,je
			jm(j)=modulo(j-1,nyg)+1
		end do
    	else
    		do j=js,je
			jm(j)=j
		end do	
    	end if
    	if (wrap_z) then
    		do k=ks,ke
			km(k)=modulo(k-1,nzg)+1
		end do
    	else
    		do k=ks,ke
			km(k)=k
		end do	
    	end if

    	i0gm=i0g
    	i1gm=i1g
    	j0gm=j0g
    	j1gm=j1g
    	k0gm=k0g
    	k1gm=k1g
    		
    	
    	idx_buff(1,:)=-1
    	n=1
    	
    	select case(ftype)
    		case(1)		
	   			is=is-1	
				i0gm=i0gm-1
				im(is)=im(is+1)-1
			case(2)
	   			js=js-1	
				j0gm=j0gm-1
				jm(js)=jm(js+1)-1
			case(3)
	   			ks=ks-1	
				k0gm=k0gm-1
				km(ks)=km(ks+1)-1
			case default
		end select	


		do k=ks,ke
			do j=js,je
				do i=is,ie
					
					il=im(i)-i0g+1
					jl=jm(j)-j0g+1
					kl=km(k)-k0g+1
					
					if ((im(i)<=i1gm).AND.(im(i)>=i0gm)) then
						if ((jm(j)<=j1gm).AND.(jm(j)>=j0gm)) then
							if ((km(k)<=k1gm).AND.(km(k)>=k0gm)) then
								idx_buff(1,n)=0
								idx_buff(2,n)=il
								idx_buff(3,n)=jl
								idx_buff(4,n)=kl
								val_buff(n)=ptrf(i,j,k)
								n=n+1							
							end if
						end if
					end if
					
				end do
			end do
		end do		
			
    						
	end subroutine
!========================================================================================!

!========================================================================================!
    subroutine packSendToMasterBuff(idx_buff,val_buff,f,ftype,is_org,ie_org,js_org,&
    						        je_org,ks_org,ke_org,offset,i0g,i1g,j0g,j1g,k0g,k1g,&
    						        nxg,nyg,nzg,wrap_x,wrap_y,wrap_z)
    	integer, allocatable, dimension(:,:), intent(inout) :: idx_buff
    	real(DP), allocatable, dimension(:), intent(inout) :: val_buff
    	real(DP), allocatable, dimension(:,:,:), intent(in) :: f
    	integer, intent(in) :: ftype,offset
    	logical, intent(in) :: wrap_x,wrap_y,wrap_z
    	integer, intent(in) :: is_org,ie_org,js_org,je_org,ks_org,ke_org
    	integer, intent(in) :: i0g,i1g,j0g,j1g,k0g,k1g
    	integer, intent(in) :: nxg,nyg,nzg
    	integer :: i0gm,i1gm,j0gm,j1gm,k0gm,k1gm
    	integer :: i,j,k,n,is,js,ks,ie,je,ke
    	integer, allocatable, dimension(:) :: im,jm,km
    	integer :: il,jl,kl
    	
		
    	!copy ranges
    	is=is_org
    	js=js_org
    	ks=ks_org
    	ie=ie_org
    	je=je_org
    	ke=ke_org

    	select case(ftype)
    		case(1)
    		    call allocateArray(im,is-1,ie)
				call allocateArray(jm,js,je)
				call allocateArray(km,ks,ke)
    		case(2)
    		    call allocateArray(im,is,ie)
				call allocateArray(jm,js-1,je)
				call allocateArray(km,ks,ke)
    		case(3)
    		    call allocateArray(im,is,ie)
				call allocateArray(jm,js,je)
				call allocateArray(km,ks-1,ke)
    		case default
    		    call allocateArray(im,is,ie)
				call allocateArray(jm,js,je)
				call allocateArray(km,ks,ke)
    	end select

    		
    	!exception periodic bc
    	if (wrap_x) then
    		do i=is+offset,ie-offset
				im(i)=modulo(i-1,nxg)+1
			end do
    	else
    		do i=is+offset,ie-offset
				im(i)=i
			end do	
    	end if
    	if (wrap_y) then
    		do j=js+offset,je-offset
				jm(j)=modulo(j-1,nyg)+1
			end do
    	else
    		do j=js+offset,je-offset
				jm(j)=j
			end do	
    	end if
    	if (wrap_z) then
    		do k=ks+offset,ke-offset
				km(k)=modulo(k-1,nzg)+1
			end do
    	else
    		do k=ks+offset,ke-offset
				km(k)=k
			end do	
    	end if
    	
    	!complete boundary indexes
    	do i=1,offset
    		im(is+offset-i)=im(is+offset)-i
    		jm(js+offset-i)=jm(js+offset)-i
    		km(ks+offset-i)=km(ks+offset)-i
    		im(ie-offset+i)=im(ie-offset)+i
    		jm(je-offset+i)=jm(je-offset)+i
    		km(ke-offset+i)=km(ke-offset)+i
    	end do

    	i0gm=i0g-offset
    	i1gm=i1g+offset
    	j0gm=j0g-offset
    	j1gm=j1g+offset
    	k0gm=k0g-offset
    	k1gm=k1g+offset
    		
    	
    	idx_buff(1,:)=-1
    	n=1
    	
    	select case(ftype)
    		case(1)		
	   			is=is-1	
				i0gm=i0gm-1
				im(is)=im(is+1)-1
			case(2)
	   			js=js-1	
				j0gm=j0gm-1
				jm(js)=jm(js+1)-1
			case(3)
	   			ks=ks-1	
				k0gm=k0gm-1
				km(ks)=km(ks+1)-1
			case default
		end select	
    	
		do k=ks,ke
			do j=js,je
				do i=is,ie

					il=im(i)-i0g+1
					jl=jm(j)-j0g+1
					kl=km(k)-k0g+1
					
					if ((im(i)<=i1gm).AND.(im(i)>=i0gm)) then
						if ((jm(j)<=j1gm).AND.(jm(j)>=j0gm)) then
							if ((km(k)<=k1gm).AND.(km(k)>=k0gm)) then
								idx_buff(1,n)=0
								idx_buff(2,n)=i
								idx_buff(3,n)=j
								idx_buff(4,n)=k
								val_buff(n)=f(il,jl,kl)
								n=n+1
							end if
						end if
					end if
					
				end do
			end do
		end do
		
				

    end subroutine   	
!========================================================================================!

!========================================================================================!
    subroutine unPackRecvFromMasterBuff(idx_buff,val_buff,f,unPack_op)
    	integer, allocatable, dimension(:,:), intent(in) :: idx_buff
    	real(DP), allocatable, dimension(:), intent(in) :: val_buff
    	real(DP), allocatable, dimension(:,:,:), intent(inout) :: f
    	integer, intent(in) :: unPack_op
    	integer :: sizeBuff,i,j,k,n
    	real(DP) :: c0,cv
    
    	sizeBuff=size(val_buff)

		select case(unPack_op)
			case(UNPACK_MAX)
    			do n=1,sizeBuff
    				if (idx_buff(1,n)==0) then
    					i=idx_buff(2,n)
    					j=idx_buff(3,n)
    					k=idx_buff(4,n)
    					c0=f(i,j,k)
    					cv=val_buff(n)
    					f(i,j,k)=max(c0,cv)
    				end if
    			end do 
    		case(UNPACK_SUM)
    			do n=1,sizeBuff
    				if (idx_buff(1,n)==0) then
    					i=idx_buff(2,n)
    					j=idx_buff(3,n)
    					k=idx_buff(4,n)
    					f(i,j,k)=f(i,j,k)+val_buff(n)
    				end if
    			end do 
    		case default
    	end select   	
    	
    end subroutine   	
!========================================================================================!

!========================================================================================!
    subroutine unPackRecvFromSlaveBuff(f,idx_buff,val_buff)
    	real(DP), allocatable, dimension(:,:,:), intent(inout) :: f 
    	integer, allocatable, dimension(:,:), intent(in) :: idx_buff
    	real(DP), allocatable, dimension(:), intent(in) :: val_buff
    	integer :: i,j,k,n,sizeBuff	
    
    	sizeBuff=size(val_buff)
    
    	do n=1,sizeBuff
    		if (idx_buff(1,n)==0) then
    			i=idx_buff(2,n)
    			j=idx_buff(3,n)
    			k=idx_buff(4,n)
    			f(i,j,k)=val_buff(n)
    		end if
    	end do


    end subroutine   	
!========================================================================================!

!========================================================================================!
    subroutine excLists(mesh)
    	type(grid), intent(in) :: mesh
    	type(mpiControl), pointer :: mpic
    	integer :: b,blk,n,nprocs,ierror,m,i,j
    	integer, allocatable, dimension(:,:) :: tmp_data,tmp_proc
    	integer :: tmp_size
    	

		mpic => mesh%ptrMPIC_
		nprocs = mpic%nProcs_
		
		call fillUpMasterList(mesh,s_blk_data,s_blk_proc)
		
		!$OMP PARALLEL DO DEFAULT(none) &
		!$OMP SHARED(nprocs,s_gbList,s_nb) &
		!$OMP PRIVATE(i,j)
		do j=1,s_nb
			do i=0,nprocs-1
				s_gbList(1:2,i,j)=-1
			end do
		end do
		!$OMP END PARALLEL DO

 		do n=0,nprocs-1
			if (n==mpic%rank_) then
				tmp_size=s_nblk			
			end if
			call MPI_BCAST(tmp_size, 1, MPI_INTEGER, n, mpic%cartComm_, ierror)
			
			if (tmp_size==0) then
 				cycle
			end if
			
			call reAllocateArray(tmp_data,1,8,1,tmp_size)
			call reAllocateArray(tmp_proc,1,tmp_size,0,nprocs-1)
			if (n==mpic%rank_) then
				tmp_data=s_blk_data
				tmp_proc=s_blk_proc
			end if
       		call MPI_BCAST(tmp_data, size(tmp_data), MPI_INTEGER, n, mpic%cartComm_, ierror)
      		call MPI_BCAST(tmp_proc, size(tmp_proc), MPI_INTEGER, n, mpic%cartComm_, ierror)
      		
        		
       		!check
       		do blk=1,tmp_size      		
    			do m=0,nprocs-1
       				if (tmp_proc(blk,m)>=0) then    				
       					b=tmp_data(1,blk)
       					s_gbList(1,m,b)=tmp_data(2,blk)
       					s_gbList(2,m,b)=m
       					s_gbList(3:,m,b)=tmp_data(3:,blk)
       				end if
       			end do
    		end do	
	
		end do

    end subroutine   	
!========================================================================================!

!========================================================================================!
    subroutine gatherLogicalExchange(mesh)
    	type(grid), intent(in) :: mesh
    	type(mpiControl), pointer :: mpic
    	logical :: s_exchange
    	integer :: ierror
    	
    	mpic => mesh%ptrMPIC_

		s_exchange = any(s_exchange_b)

		call MPI_ALLGather(s_exchange,1,MPI_LOGICAL,s_exchange_g,1,MPI_LOGICAL,&
						   mpic%cartComm_,ierror)
		
     end subroutine   	
!========================================================================================!

!========================================================================================!
    subroutine fillUpMasterList(mesh,blk_data,blk_proc)
    	type(grid), intent(in) :: mesh
    	integer, allocatable, dimension(:,:), intent(inout) :: blk_data,blk_proc
    	integer :: is,ie,js,je,ks,ke,bn
    	integer :: ism,jsm,ksm
    	type(mpiControl), pointer :: mpic
    	integer :: i0g,i1g,j0g,j1g,k0g,k1g,nxg,nyg,nzg
    	integer :: i,n,nproc,b
    	logical :: wrap_x,wrap_y,wrap_z
    	logical :: isX,isY,isZ
    	
    	
		mpic => mesh%ptrMPIC_
		nproc = mpic%nProcs_
		wrap_x = mpic%wrapAround_(1)
		wrap_y = mpic%wrapAround_(2)
		wrap_z = mpic%wrapAround_(3)
		
		nxg=mesh%nxg_
		nyg=mesh%nyg_
		nzg=mesh%nzg_

		blk_proc=-1

		do n=0,nproc-1

			call globalIndexesFromAll(mesh,i0g,i1g,j0g,j1g,k0g,k1g,n,'cl')

			do b=1,s_nblk
			
				is=vofBlocks(b)%idx(1)-1
				ie=vofBlocks(b)%idx(2)+1
				js=vofBlocks(b)%idx(3)-1
				je=vofBlocks(b)%idx(4)+1
				ks=vofBlocks(b)%idx(5)-1
				ke=vofBlocks(b)%idx(6)+1
				
				bn=vofBlocks(b)%bn
				
				isX=.FALSE.
				isY=.FALSE.
				isZ=.FALSE.
			
				!check x
				do i=is,ie
					ism=vofBlocks(b)%idx_mx(i)
					if ((ism<=i1g).AND.(ism>=i0g)) then
						isX=.TRUE.
						exit
					end if
				end do
				!check y
				do i=js,je
					jsm=vofBlocks(b)%idx_my(i)
					if ((jsm<=j1g).AND.(jsm>=j0g)) then
						isY=.TRUE.
						exit
					end if
				end do
				!check z
				do i=ks,ke
					ksm=vofBlocks(b)%idx_mz(i)
					if ((ksm<=k1g).AND.(ksm>=k0g)) then
						isZ=.TRUE.
						exit
					end if
				end do
				
				if ((isX).AND.(isY).AND.(isZ)) then
					!fill up block data
					blk_data(1,b)=bn
					blk_data(2,b)=vofBlocks(b)%master
					blk_data(3:,b)=vofBlocks(b)%idx		
					blk_proc(b,n)=1				
				end if
				
				
			end do
		
		end do

    end subroutine   	
!========================================================================================!

!========================================================================================!
    subroutine reInitBlockDistribution(mesh,gmesh)
   		type(grid), intent(in) :: mesh,gmesh
		real(DP), allocatable, dimension(:,:,:) :: c_blk,c0_blk
		integer :: b,n,nb,nprocs,bl,is,ie,js,je,ks,ke,i0g,i1g,j0g,j1g,k0g,k1g,nxg,nyg,nzg
		integer :: master,slave,tag,ierror,sizeBuff,position,ntria
    	integer, dimension(2) :: requests
    	integer, dimension(MPI_STATUS_SIZE,2) :: status
		logical :: isHere
		type(vofBlock), allocatable, dimension(:) :: new_blocks
		type(mpiControl), pointer :: mpic
		character(len=:), allocatable :: buff
    		
    		
		mpic => mesh%ptrMPIC_
		nprocs = mpic%nProcs_
		tag=0
		
		nxg = mesh%nxg_
		nyg = mesh%nyg_
		nzg = mesh%nzg_

    	do b=1,s_nb
    	  	
    		do n=0,nprocs-1
    		
    			!manage own bubbls
    			if ((mpic%rank_==s_gbList(1,n,b)).AND.(s_gbList(1,n,b)==s_gbList(2,n,b))) then
    			
    				i0g = mesh%i0g_
    				i1g = mesh%i1g_
    				j0g = mesh%j0g_
    				j1g = mesh%j1g_
    				k0g = mesh%k0g_
    				k1g = mesh%k1g_
    				
  					is = s_gbList(3,n,b)
					ie = s_gbList(4,n,b)
					js = s_gbList(5,n,b)
					je = s_gbList(6,n,b)
					ks = s_gbList(7,n,b)
					ke = s_gbList(8,n,b)
    			
					call isBubbleCenterHere(is,ie,js,je,ks,ke,i0g,i1g,j0g,j1g,k0g,k1g,&
    						      		    nxg,nyg,nzg,isHere)
    			
    				if (isHere) then
						!add new block
						call addNewBlock(new_blocks,nb)
					
						!copy old block
    					do bl=1,s_nblk
    						if (vofBlocks(bl)%bn==b) then
								new_blocks(nb)=vofBlocks(bl)
    						end if
    					end do
    				end if
    								
				
				!this proc is the a slave (recv)
				else if (mpic%rank_==s_gbList(2,n,b)) then
				
					master = s_gbList(1,n,b)
					slave = s_gbList(2,n,b)
					
    				i0g = mesh%i0g_
    				i1g = mesh%i1g_
    				j0g = mesh%j0g_
    				j1g = mesh%j1g_
    				k0g = mesh%k0g_
    				k1g = mesh%k1g_

  					is = s_gbList(3,n,b)
					ie = s_gbList(4,n,b)
					js = s_gbList(5,n,b)
					je = s_gbList(6,n,b)
					ks = s_gbList(7,n,b)
					ke = s_gbList(8,n,b)
								
					call isBubbleCenterHere(is,ie,js,je,ks,ke,i0g,i1g,j0g,j1g,k0g,k1g,&
    						      		    nxg,nyg,nzg,isHere)
    				
    				!aggiungi element
    				if (isHere) then
    				    
    				    !recv c,c0
    				    call reAllocateArray(c_blk,is,ie,js,je,ks,ke)
    				    call reAllocateArray(c0_blk,is,ie,js,je,ks,ke)
    				    sizeBuff=realDP_size*(size(c_blk)+size(c0_blk))
    				    call reAllocateArray(buff,sizeBuff)
    				    
    				    call MPI_Recv(buff,sizeBuff,MPI_CHARACTER,master,tag,&
    				    		      mpic%cartComm_,status,ierror)
    				    position = 0
    				    call MPI_UNPACK(buff, sizeBuff, position, c_blk, size(c_blk), &
								        MPI_DOUBLE_PRECISION, mpic%cartComm_, ierror)
    				    call MPI_UNPACK(buff, sizeBuff, position, c0_blk, size(c0_blk), &
								        MPI_DOUBLE_PRECISION, mpic%cartComm_, ierror)			
    				
    					!add block element
						call addNewBlock(new_blocks,nb)
					
						!init new block
						call initBlock(mesh,gmesh,new_blocks(nb),b,slave,is,ie,js,je,ks,ke,&
									   REDISTRIBUTION_BLK,c_blk,c0_blk)

    				end if      		    
    				
				!this proc is the master
				else if (mpic%rank_==s_gbList(1,n,b)) then
				
					slave = s_gbList(2,n,b)

					call globalIndexesFromAll(mesh,i0g,i1g,j0g,j1g,k0g,k1g,slave,'cl')	

  					is = s_gbList(3,n,b)
					ie = s_gbList(4,n,b)
					js = s_gbList(5,n,b)
					je = s_gbList(6,n,b)
					ks = s_gbList(7,n,b)
					ke = s_gbList(8,n,b)
								
					call isBubbleCenterHere(is,ie,js,je,ks,ke,i0g,i1g,j0g,j1g,k0g,k1g,&
    						      		    nxg,nyg,nzg,isHere)	

					if (isHere) then
    					!send c
    					call reAllocateArray(c_blk,is,ie,js,je,ks,ke)
    					call reAllocateArray(c0_blk,is,ie,js,je,ks,ke)
    					do bl=1,s_nblk
    						if (vofBlocks(bl)%bn==b) then
								c_blk=vofBlocks(bl)%c(is:ie,js:je,ks:ke)
								c0_blk=vofBlocks(bl)%c0(is:ie,js:je,ks:ke)
    						end if
    					end do    	
    					
    					sizeBuff=realDP_size*(size(c_blk)+size(c0_blk))
    					call reAllocateArray(buff,sizeBuff)
    					position = 0
    					call MPI_PACK(c_blk,size(c_blk), MPI_DOUBLE_PRECISION, buff, &
    								  sizeBuff, position, mpic%cartComm_, ierror)	
    					call MPI_PACK(c0_blk,size(c0_blk), MPI_DOUBLE_PRECISION, buff, &
    								  sizeBuff, position, mpic%cartComm_, ierror) 		
						call MPI_SEND(buff,sizeBuff,MPI_CHARACTER, slave, tag, &
						              mpic%cartComm_, ierror)
					end if		

				end if
			
			end do
			
		end do
		

 		!set new blocks
    	if (allocated(new_blocks)) then
			s_nblk = size(new_blocks)
			call copyBlocks(vofBlocks,new_blocks)
			call deallocateBlocks(new_blocks)
		else
			s_nblk=0
			if (allocated(vofBlocks)) then
				call deallocateBlocks(vofBlocks)
			end if
		end if
		
		!update exchange list
		call reAllocateArray(s_blk_data,1,8,1,s_nblk)
		call reAllocateArray(s_blk_proc,1,s_nblk,0,nprocs-1)
		call excLists(mesh)
		
		!reinit exchange logical
		if (s_nblk>0) then
			call reAllocateArray(s_exchange_b,1,s_nblk)
		else
			call reAllocateArray(s_exchange_b,1,1)
		end if
		s_exchange_g(mpic%rank_)=.FALSE.
		s_exchange_b = .FALSE.
		
		call measureBlocksDistr(mpic)
		
    end subroutine   	
!========================================================================================!

!========================================================================================!
    subroutine isBubbleCenterHere(is,ie,js,je,ks,ke,i0g,i1g,j0g,j1g,k0g,k1g,&
    						      nxg,nyg,nzg,isHere)
    	integer, intent(in) :: is,ie,js,je,ks,ke,i0g,i1g,j0g,j1g,k0g,k1g,nxg,nyg,nzg
		logical, intent(out) :: isHere
		integer :: ic,jc,kc,icm,jcm,kcm
		logical :: isX,isY,isZ
		
			ic=nint(0.5d0*(is+ie))
			jc=nint(0.5d0*(js+je))
			kc=nint(0.5d0*(ks+ke))
					
			icm=modulo(ic-1,nxg)+1
			jcm=modulo(jc-1,nyg)+1
			kcm=modulo(kc-1,nzg)+1
				
			isX=.FALSE.
			isY=.FALSE.
			isZ=.FALSE.
			
			!check x
			if ((icm<=i1g).AND.(icm>=i0g)) then
				isX=.TRUE.
			end if
			!check y
			if ((jcm<=j1g).AND.(jcm>=j0g)) then
				isY=.TRUE.
			end if
			!check z
			if ((kcm<=k1g).AND.(kcm>=k0g)) then
				isZ=.TRUE.
			end if
			
			if ((isX).AND.(isY).AND.(isZ)) then
				isHere=.TRUE.
			else
				isHere=.FALSE.
			end if
				
					
    end subroutine   	
!========================================================================================!

!========================================================================================!
    subroutine isBubbleInThisCore(ic,jc,kc,i0g,i1g,j0g,j1g,k0g,k1g,bool)
    	integer, intent(in) :: ic,jc,kc,i0g,i1g,j0g,j1g,k0g,k1g
    	logical, intent(out) :: bool
    	
    	if (((ic>=i0g).AND.(ic<=i1g)) .AND. &
    	   ((jc>=j0g).AND.(jc<=j1g))  .AND. &
    	   ((kc>=k0g).AND.(kc<=k1g))) then
    	   	bool = .TRUE.
    	else
    		bool = .FALSE.
    	end if
   	
     end subroutine   	
!========================================================================================!

!========================================================================================!
    subroutine printVOFblocks(nFolder)
    	integer, intent(in) :: nFolder
        character(len=20) :: dirName,bn
        integer :: i
        
        if (allocated(vofBlocks)) then
        
        write(dirName,s_intFormat) nFolder
        
		do i=1,s_nblk
						
			write(bn,s_intFormat) vofBlocks(i)%bn
        	open(UNIT=s_IOunitNumber,FILE=trim(adjustl(dirName))//'/b'//trim(adjustl(bn)), &
        	 	      form='UNFORMATTED',ACCESS='STREAM',STATUS='REPLACE',ACTION='WRITE')

        		write(s_IOunitNumber) vofBlocks(i)%master
        		write(s_IOunitNumber) vofBlocks(i)%bn
        		write(s_IOunitNumber) vofBlocks(i)%idx
        		write(s_IOunitNumber) vofBlocks(i)%c0
			
			close(s_IOunitNumber)	
					
		end do
		
		end if
        
     end subroutine   	
!========================================================================================!

!========================================================================================!
    subroutine bubblesLagrQ(t)
		real(DP), intent(in) :: t
		integer :: b,i,j,k,is,ie,js,je,ks,ke
		real(DP) :: c,Vol,Vol_blk,x,y,z,u,v,w
		character(len=6) :: bn
		character(len=100) :: dir
		logical :: exist

		do b=1,s_nblk
			Vol_blk=0.d0
			x=0.d0
			y=0.d0
			z=0.d0
			u=0.d0
			v=0.d0
			w=0.d0
			
			is=vofBlocks(b)%idx(1)
			ie=vofBlocks(b)%idx(2)
			js=vofBlocks(b)%idx(3)
			je=vofBlocks(b)%idx(4)
			ks=vofBlocks(b)%idx(5)
			ke=vofBlocks(b)%idx(6)
			
			do k=ks,ke
				do j=js,je
					do i=is,ie
						c=vofBlocks(b)%c(i,j,k)
						Vol=vofBlocks(b)%dxf(i)*vofBlocks(b)%dyf(j)*vofBlocks(b)%dzf(k)
						!volume
						Vol_blk=Vol_blk+c*Vol
						!center
						x=x+c*vofBlocks(b)%xc(i)*Vol
						y=y+c*vofBlocks(b)%yc(j)*Vol
						z=z+c*vofBlocks(b)%zc(k)*Vol
						!velocity
						u=u+c*0.5d0*(vofBlocks(b)%ux(i,j,k)+vofBlocks(b)%ux(i-1,j,k))*Vol
						v=v+c*0.5d0*(vofBlocks(b)%uy(i,j,k)+vofBlocks(b)%uy(i,j-1,k))*Vol
						w=w+c*0.5d0*(vofBlocks(b)%uz(i,j,k)+vofBlocks(b)%uz(i,j,k-1))*Vol				
					end do
				end do
			end do
			
			x=x/Vol_blk
			y=y/Vol_blk
			z=z/Vol_blk
			u=u/Vol_blk
			v=v/Vol_blk
			w=w/Vol_blk
			
			!write file
			write(bn,s_intFormat) vofBlocks(b)%bn
			dir="bubblesLagrQ"//'/b'//trim(adjustl(bn))
			
			inquire(file=dir, exist=exist)
			
			if (exist) then
				open(UNIT=s_IOunitNumber,FILE=dir,STATUS='OLD',POSITION="append",ACTION='WRITE')
					write(s_IOunitNumber,'(7'//s_doubleFormat(2:10)//')') t,x,y,z,u,v,w
				close(s_IOunitNumber)	
			else
				open(UNIT=s_IOunitNumber,FILE=dir,STATUS='NEW',ACTION='WRITE')
					write(s_IOunitNumber,'(7'//s_doubleFormat(2:10)//')') t,x,y,z,u,v,w
				close(s_IOunitNumber)				
			end if
			
		end do
		
			
        
    end subroutine
!========================================================================================!

!========================================================================================!
    subroutine max_box_size(mesh)
    	type(grid), intent(in) :: mesh
    	type(mpiControl), pointer :: comm
    	integer :: b,is,ie,js,je,ks,ke,sx,sy,sz,sx_g,sy_g,sz_g,ierror
 
 		comm => mesh%ptrMPIC_

		sx=0
		sy=0
		sz=0

		do b=1,s_nblk	

			is=vofBlocks(b)%idx(1)
			ie=vofBlocks(b)%idx(2)
			js=vofBlocks(b)%idx(3)
			je=vofBlocks(b)%idx(4)
			ks=vofBlocks(b)%idx(5)
			ke=vofBlocks(b)%idx(6)		
	
			sx=max(ie-is+1,sx)
			sy=max(je-js+1,sy)
			sz=max(ke-ks+1,sz)

		end do
		
		call Mpi_Reduce(sx, sx_g, 1, MPI_INTEGER, MPI_MAX, 0, comm%cartComm_, ierror)
		call Mpi_Reduce(sy, sy_g, 1, MPI_INTEGER, MPI_MAX, 0, comm%cartComm_, ierror)
		call Mpi_Reduce(sz, sz_g, 1, MPI_INTEGER, MPI_MAX, 0, comm%cartComm_, ierror)
		
		if (IS_MASTER) then
			write(*,*) 'MAX BOX SIZE X: ', sx_g
			write(*,*) 'MAX BOX SIZE Y: ', sy_g
			write(*,*) 'MAX BOX SIZE Z: ', sz_g
		end if
		
		
    end subroutine
!========================================================================================!	

end module vofBlocksMod

