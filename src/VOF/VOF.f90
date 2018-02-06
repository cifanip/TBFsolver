! ************************************************************************************** !
!    TBFsim - DNS turbulent bubbly flow simulator
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

module VOFMod

	use initialConditionsMod, only: pi
	use timeMod
	
	implicit none
	
	integer, parameter :: s_Xsweep = 1
	integer, parameter :: s_Ysweep = 2
	integer, parameter :: s_Zsweep = 3
	
	integer, protected :: s_sweepSelector = 0
	integer, dimension(3), protected :: s_sweep

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
		real(DP), allocatable, dimension(:,:,:) :: c, cs, cv, c0
		real(DP), allocatable, dimension(:,:,:) :: ux, uy, uz
		real(DP), allocatable, dimension(:,:,:) :: stx, sty, stz
		real(DP), allocatable, dimension(:,:,:) :: nx, ny, nz
		logical, allocatable, dimension(:,:,:) :: isMixed, isFull
		real(DP), allocatable, dimension(:,:,:) :: q
		real(DP), allocatable, dimension(:,:,:) :: k
		real(DP), allocatable, dimension(:,:,:) :: geoFlux, corrFlux, corrTerm
	end type
	
    interface assignment(=)
       module procedure assignBlock
    end interface
	
	type(vofBlock), allocatable, dimension(:) :: vofBlocks
	
	
	type, public :: VOF
	
		type(grid), pointer :: mesh_ => NULL()
		type(grid), pointer :: gmesh_ => NULL()
		type(time), pointer :: ptrTime_ => NULL() 
		
		!threshold vof field
		real(DP), private :: eps_ = 1.d-6
		
		!halo dim of box volume fraction field
		integer, private :: hd_
		
		!mat props
		real(DP) :: rhol_, mul_, rhog_, mug_
		real(DP) :: sigma_	
		
	end type
	
	!vertex field
	real(DP), allocatable, dimension(:,:,:), protected :: s_cv

	private :: storeOldVOField
	private :: cnh
	private :: cnp
	private :: reconstruct
	private :: implicitIntegratorAndUpdateCorrTerm
	private :: explicitIntegrator
	private :: computeCorrFluxes
	private :: computeGeoFluxes
	private :: computeGeoMixedFluxes
	private :: updateState
	private :: updateStateBlocks
	private :: reSetFullEmpty
	private :: clump
	private :: sweepCombination
	private :: info
	private :: bubbles_distribution_space
	private :: init_indexes_box
	private :: readSingleBubble
	private :: readBubblesArray
	private :: init_Bubble_vf
	private :: excLists
	private :: allocateBlockSendRecvBuffers
	private :: packSendToSlaveBuff
	private :: packSendToSlaveBuff_st
	private :: packSendToMasterBuff
	private :: unPackRecvFromMasterBuff
	private :: unPackRecvFromSlaveBuff
	private :: excVel
	private :: excVFK
	private :: excST
	private :: fillUpMasterList
	private :: updateBlock
	private :: reInitBlock
	private :: initBlock
	private :: reInitBlockDistribution
	private :: isBubbleCenterHere
	private :: addNewBlock
	private :: copyBlocks
	private :: assignBlock
	private :: resetFragments
	private :: checkBlockSize
	private :: copyBlockMesh
	private :: computeModuloIdx
	private :: isBubbleInThisCore
	private :: gatherLogicalExchange
	private :: bubblesLagrQ
	private :: computeNormal_sharp
	private :: computeNormal_smooth
	private :: computeNormal_youngs
	private :: computeNormal_hf
	private :: hfNormal
	private :: smoothVFblock
	private :: smoothVF
	private :: correctContantAngle
	private :: computeBlockST
	private :: computeCurvature
	private :: spreadCurvature
	private :: computeBlockCurvature
	private :: hfColumn
	private :: hfCurvature
	private :: searchBlockHF
	private :: independentPoints
	private :: localIntPos
	private :: buildLSqSystem
	private :: solveLSqSystem
	private :: interCentroids
	private :: parabFittedCurvature
	private :: fdCurvature
	private :: qCell
	private :: vCell
	private :: vArea
	private :: centroid
	private :: intersPoints
	private :: intersPlaneSegm
	private :: setColumnPar
	private :: setStencilPar
	private :: setDirections
	private :: n_interfaces
	private :: max_box_size

	public :: vofCTOR
	public :: initVOFblocks
	public :: solveVOF
	public :: updateMaterialProps
	public :: computeSurfaceTension
	public :: printVOFblocks
	public :: deallocateBlocks


contains

!========================================================================================!
    subroutine vofCTOR(this,gmesh,mesh,rt)
        type(VOF), intent(out) :: this
        type(time), intent(in), target :: rt
        type(grid), intent(in), target :: gmesh, mesh
        integer :: nprocs
        type(dictionary) :: dict
        integer :: nx,ny,nz
        
        this%ptrTime_ => rt
        this%mesh_ => mesh
        this%gmesh_ => gmesh
        
        nx=mesh%nx_
        ny=mesh%ny_
        nz=mesh%nz_
		
		this%hd_ = 3
		
		!read mat props
		call dictionaryCTOR(dict ,'parameters','specs')
		call readParameter(dict,this%rhol_,'rhol')
		call readParameter(dict,this%rhog_,'rhog')
		call readParameter(dict,this%mul_,'mul')
		call readParameter(dict,this%mug_,'mug')
		call readParameter(dict,this%sigma_,'sigma')

		!init blocks
		call initVOFblocks(mesh,gmesh)
		
		!store old c field
    	call storeOldVOField()

		nprocs=mesh%ptrMPIC_%nProcs_
    	call allocateArray(s_gbList,1,8,0,nprocs-1,1,s_nb)
    	
    	call reAllocateArray(s_blk_data,1,8,1,s_nblk)
		call reAllocateArray(s_blk_proc,1,s_nblk,0,nprocs-1)
    	
    	!allocate vertex field
    	call allocateArray(s_cv,0,nx,0,ny,0,nz)
    	
		call excLists(mesh)
		
		call allocateArray(s_exchange_g,0,nprocs-1)
		if (s_nblk>0) then
			call allocateArray(s_exchange_b,1,s_nblk)
		else
			call allocateArray(s_exchange_b,1,1)		
		end if
		s_exchange_g(mesh%ptrMPIC_%rank_)=.FALSE.
		s_exchange_b = .FALSE.	
        
        call updateStateBlocks(this)
		
    end subroutine
!========================================================================================!

!********************************** VOF-advection ***************************************!

!========================================================================================!
    subroutine solveVOF(this,c,u)
    	type(VOF), intent(inout) :: this
    	type(scalarField), intent(inout) :: c
    	type(vectorField), intent(in) :: u
    	integer :: i,b
    	real(DP) :: start, finish
    	
    	start = MPI_Wtime()
    	
    	!set c=c^{n+1}
    	call cnp()

		call sweepCombination()
		
		!exchange velocity field
		call excVel(this%mesh_,u)
		
		!output Lagrangian quantities bubbles
		!call bubblesLagrQ(this%ptrTime_%t_)
		

		!$OMP PARALLEL DO DEFAULT(none) &
		!$OMP SHARED(vofBlocks,this,s_nblk,s_sweep) &
		!$OMP PRIVATE(i,b)
		do b=1,s_nblk
		
			!reset corr term
			vofBlocks(b)%corrTerm=0.d0
			
			do i=1,3
		
				call reconstruct(this,vofBlocks(b))
				call computeGeoFluxes(this,vofBlocks(b),s_sweep(i))
				call computeCorrFluxes(this,vofBlocks(b),s_sweep(i))
				call implicitIntegratorAndUpdateCorrTerm(this,vofBlocks(b))

			end do
		
			call explicitIntegrator(this,vofBlocks(b))	
			call resetFragments(vofBlocks(b))
			!compute c^{n=1/2}
			call cnh(vofBlocks(b))	
			call updateBlock(this,vofBlocks(b),b)
			call updateState(this,vofBlocks(b))
			call reconstruct(this,vofBlocks(b))
		
		end do
		!$OMP END PARALLEL DO
	
		call gatherLogicalExchange(this)
		
		!update list
		if (any(s_exchange_g)) then
			call excLists(this%mesh_)
		end if
	
		!exchange volume fraction field
		call excVFK(this%mesh_,c,1)
		
		!exchange velocity for post-proc only
		call excVel(this%mesh_,u)
		
		!redistribute blocks
		if (vofBlocksRed(this%ptrTime_)) then
			if (IS_MASTER) then
				write(*,*) '	Blocks redistribution called'
			end if
			call reInitBlockDistribution(this)
		end if
		

		finish = MPI_Wtime()
		
		call info(this,finish-start,0)
    	
     end subroutine   	
!========================================================================================!

!========================================================================================!
    subroutine storeOldVOField()
    	integer :: b
    	
		do b=1,s_nblk
			vofBlocks(b)%c0=vofBlocks(b)%c
		end do  		
    
    end subroutine
!========================================================================================!

!========================================================================================!
    subroutine cnh(vofb)
    	type(vofBlock), intent(inout) :: vofb
    	real(DP), allocatable, dimension(:,:,:) :: tmp
    	integer :: lbi,ubi,lbj,ubj,lbk,ubk

		lbi = lbound(vofb%c,1)
		ubi = ubound(vofb%c,1)
		lbj = lbound(vofb%c,2)
		ubj = ubound(vofb%c,2)
		lbk = lbound(vofb%c,3)
		ubk = ubound(vofb%c,3)
		
		call allocateArray(tmp,lbi,ubi,lbj,ubj,lbk,ubk)
    	
		tmp=0.5d0*(vofb%c+vofb%c0)
		vofb%c0=vofb%c
		vofb%c=tmp

    end subroutine
!========================================================================================!

!========================================================================================!
    subroutine cnp()
    	integer :: b
    	
		do b=1,s_nblk
			vofBlocks(b)%c=vofBlocks(b)%c0
		end do 	 			
    
    end subroutine
!========================================================================================!

!========================================================================================!
    subroutine reconstruct(this,vofb)
    	type(VOF), intent(in) :: this
    	type(vofBlock), intent(inout) :: vofb
    	integer :: i, j, k
    	integer :: lbi,ubi,lbj,ubj,lbk,ubk
    	real(DP), dimension(3) :: m, delta
    	real(DP) :: cv
    	real(DP) :: qv
    	
    	
		lbi = lbound(vofb%q,1)
		ubi = ubound(vofb%q,1)
		lbj = lbound(vofb%q,2)
		ubj = ubound(vofb%q,2)
		lbk = lbound(vofb%q,3)
		ubk = ubound(vofb%q,3)
		

		!not needed
        call smoothVFblock(vofb)
    	
		call computeNormal_youngs(this,vofb)
		!call computeNormal_sharp(this,vofb)
		
    	do k=lbk,ubk
    		do j=lbj,ubj
    			do i=lbi,ubi
    				
    				if (vofb%isMixed(i,j,k)) then
    				
    					m(1) = vofb%nx(i,j,k)
    					m(2) = vofb%ny(i,j,k)
    					m(3) = vofb%nz(i,j,k)
    					
    					delta(1) = vofb%dxf(i)
    					delta(2) = vofb%dyf(j)
    					delta(3) = vofb%dzf(k)
    					
    					cv = vofb%c(i,j,k)
    					
    					call qCell(vofb,i,j,k,m,cv,delta,qv)
    					vofb%q(i,j,k) = qv
    				
    				end if
    				
    			end do
    		end do
    	end do		
    
    end subroutine
!========================================================================================!

!========================================================================================!
    subroutine implicitIntegratorAndUpdateCorrTerm(this,vofb)
    	type(VOF), intent(in) :: this
    	type(vofBlock), intent(inout) :: vofb
    	integer :: lbi,ubi,lbj,ubj,lbk,ubk
    	
		lbi = lbound(vofb%geoFlux,1)
		ubi = ubound(vofb%geoFlux,1)
		lbj = lbound(vofb%geoFlux,2)
		ubj = ubound(vofb%geoFlux,2)
		lbk = lbound(vofb%geoFlux,3)
		ubk = ubound(vofb%geoFlux,3)
    	
    	vofb%c(lbi:ubi,lbj:ubj,lbk:ubk) = (vofb%c(lbi:ubi,lbj:ubj,lbk:ubk) + vofb%geoFlux)/	&
    							          (1.d0 - vofb%corrFlux)
	     
    	vofb%corrTerm = vofb%corrTerm + vofb%c(lbi:ubi,lbj:ubj,lbk:ubk)*vofb%corrFlux

    	
    	call clump(vofb)
    	call updateState(this,vofb)
    	call reSetFullEmpty(vofb)
    	
     end subroutine   	
!========================================================================================!

!========================================================================================!
    subroutine explicitIntegrator(this,vofb)
    	type(VOF), intent(in) :: this
    	type(vofBlock), intent(inout) :: vofb
    	integer :: lbi,ubi,lbj,ubj,lbk,ubk
    	
		lbi = lbound(vofb%geoFlux,1)
		ubi = ubound(vofb%geoFlux,1)
		lbj = lbound(vofb%geoFlux,2)
		ubj = ubound(vofb%geoFlux,2)
		lbk = lbound(vofb%geoFlux,3)
		ubk = ubound(vofb%geoFlux,3)
    	
    	
    	vofb%c(lbi:ubi,lbj:ubj,lbk:ubk) = vofb%c(lbi:ubi,lbj:ubj,lbk:ubk) - vofb%corrTerm 
    	    	
    	call clump(vofb)
    	call updateState(this,vofb)
    	call reSetFullEmpty(vofb)


     end subroutine   	
!========================================================================================!

!========================================================================================!
    subroutine computeCorrFluxes(this,vofb,dir)
    	type(VOF), intent(in) :: this
    	type(vofBlock), intent(inout) :: vofb
    	integer, intent(in) :: dir
    	integer :: lbi,ubi,lbj,ubj,lbk,ubk
    	real(DP) :: vl, vr, dc
    	real(DP) :: dt,alpha,rdt
    	integer :: i, j, k

		lbi = lbound(vofb%geoFlux,1)
		ubi = ubound(vofb%geoFlux,1)
		lbj = lbound(vofb%geoFlux,2)
		ubj = ubound(vofb%geoFlux,2)
		lbk = lbound(vofb%geoFlux,3)
		ubk = ubound(vofb%geoFlux,3)
    	
    	dt = this%ptrTime_%dt_
    	alpha = alphaRKS(this%ptrTime_)
    	rdt = alpha*dt
    	
    	do k=lbk,ubk
    		do j=lbj,ubj
    			do i=lbi,ubi
 
    					!select right and left vel and delta along dir
    					select case(dir)
    						case(1)
    							vr = vofb%ux(i,j,k)
    							vl = vofb%ux(i-1,j,k)
    							dc  = vofb%dxf(i)
    						case(2)
    							vr = vofb%uy(i,j,k)
    							vl = vofb%uy(i,j-1,k)
    							dc  = vofb%dyf(j)					
    						case(3)
    							vr = vofb%uz(i,j,k)
    							vl = vofb%uz(i,j,k-1)
    							dc  = vofb%dzf(k)
    						case default   						
    					end select   
    					
    					vofb%corrFlux(i,j,k) = rdt*(vr-vl)/dc
    				
    			end do
    		end do
    	end do
    	
     end subroutine   	
!========================================================================================!

!========================================================================================!
    subroutine computeGeoFluxes(this,vofb,dir)
    	type(VOF), intent(in) :: this
    	type(vofBlock), intent(inout) :: vofb
    	integer, intent(in) :: dir
    	integer :: ip,jp,kp,im,jm,km
    	integer :: i,j,k
    	integer :: lbi,ubi,lbj,ubj,lbk,ubk
    	real(DP) :: vl, vr
    	real(DP) :: dc, dt, alpha, rdt, fr, fl
    	

		lbi = lbound(vofb%geoFlux,1)
		ubi = ubound(vofb%geoFlux,1)
		lbj = lbound(vofb%geoFlux,2)
		ubj = ubound(vofb%geoFlux,2)
		lbk = lbound(vofb%geoFlux,3)
		ubk = ubound(vofb%geoFlux,3)
    	
    	dt = this%ptrTime_%dt_
    	alpha = alphaRKS(this%ptrTime_)
    	rdt = alpha*dt
    	
    	!indexes adjacents cells
    	select case(dir)
    		case(1)
    			ip = 1
    			im = -1
    			jp = 0
    			jm = 0
    			kp = 0
    			km = 0
    		case(2)
    			ip = 0
    			im = 0
    			jp = 1
    			jm = -1
    			kp = 0
    			km = 0  						
    		case(3)
    			ip = 0
    			im = 0
    			jp = 0
    			jm = 0
    			kp = 1
    			km = -1
    		case default 					
    	end select


    	do k=lbk,ubk
    		do j=lbj,ubj
    			do i=lbi,ubi
    			
    				    !reset flux
    					vofb%geoFlux(i,j,k) = 0.d0
    			
    					!select right and left vel and delta along dir
    					select case(dir)
    						case(1)
    							vr = vofb%ux(i,j,k)
    							vl = vofb%ux(i-1,j,k)
    							dc  = vofb%dxf(i)
    						case(2)
    							vr = vofb%uy(i,j,k)
    							vl = vofb%uy(i,j-1,k)
    							dc  = vofb%dyf(j)						
    						case(3)
    							vr = vofb%uz(i,j,k)
    							vl = vofb%uz(i,j,k-1)
    							dc  = vofb%dzf(k)
    						case default    						
    					end select
     	    					
     	    			
     	    			!right face	
    					if (vr >= 0.d0) then
    					
    						if (vofb%isMixed(i,j,k)) then
    							call computeGeoMixedFluxes(this,vofb,dir,i,j,k,i,j,k,vr,fr)
    						else if (vofb%isFull(i,j,k)) then
    							fr = vr*rdt/dc
    						else
    							fr = 0.d0
    						end if
    					
    						vofb%geoFlux(i,j,k) = vofb%geoFlux(i,j,k) - fr
    					
    					else
    				
    						if (vofb%isMixed(i+ip,j+jp,k+kp)) then
    							call computeGeoMixedFluxes(this,vofb,dir,i,j,k,i+ip,j+jp,k+kp,vr,fr)
    						else if (vofb%isFull(i+ip,j+jp,k+kp)) then
    							fr = -vr*rdt/dc
    						else
    							fr = 0.d0
    						end if
    					
    						vofb%geoFlux(i,j,k) = vofb%geoFlux(i,j,k) + fr
    					
    					end if
    					
    				
    					!left face
    					if (vl >= 0.d0) then
    					
    						if (vofb%isMixed(i+im,j+jm,k+km)) then
    							call computeGeoMixedFluxes(this,vofb,dir,i,j,k,i+im,j+jm,k+km,vl,fl)
    						else if (vofb%isFull(i+im,j+jm,k+km)) then
    							fl = vl*rdt/dc
    						else
    							fl = 0.d0
    						end if
    					
    						vofb%geoFlux(i,j,k) = vofb%geoFlux(i,j,k) + fl
    					
    					else
    				
    						if (vofb%isMixed(i,j,k)) then
    							call computeGeoMixedFluxes(this,vofb,dir,i,j,k,i,j,k,vl,fl)
    						else if (vofb%isFull(i,j,k)) then
    							fl = -vl*rdt/dc
    						else
    							fl = 0.d0
    						end if
    					
    						vofb%geoFlux(i,j,k) = vofb%geoFlux(i,j,k) - fl
    				
    					end if
    					
    				
    			end do
    		end do
    	end do
    
    end subroutine
!========================================================================================!

!========================================================================================!
    subroutine computeGeoMixedFluxes(this,vofb,dir,i,j,k,iu,ju,ku,vel,flux)
    	type(VOF), intent(in) :: this
    	type(vofBlock), intent(inout) :: vofb
    	integer, intent(in) :: dir
		integer, intent(in) :: i,j,k		!indexes of the considered cell
		integer, intent(in) :: iu,ju,ku 	!indexes of the upWind cell
		real(DP), intent(in) :: vel
		real(DP), intent(out) :: flux
    	real(DP), dimension(3) :: delta, m
    	real(DP) :: q, c, dt, alpha, rdt, V
    	
    	
    	dt = this%ptrTime_%dt_
    	alpha = alphaRKS(this%ptrTime_)
    	rdt = dt*alpha
    		
    	if (vel >= 0.d0) then
    		!store cell size
    		delta(1) = vofb%dxf(iu)
    		delta(2) = vofb%dyf(ju)
    		delta(3) = vofb%dzf(ku)
    				
    		!select normal component
    		m(1) = vofb%nx(iu,ju,ku)
    		m(2) = vofb%ny(iu,ju,ku)
    		m(3) = vofb%nz(iu,ju,ku)  
    					
    		!select q
    		q = vofb%q(iu,ju,ku) 
    					
	 		!translation
			q = q + m(dir)*vel*rdt
	 		q = q - m(dir)*delta(dir)
	 		!protruded portion  
	 		delta(dir)=vel*rdt 
		    call vCell(vofb,iu,ju,ku,m,q,delta,c) 		  
		     								
    	else
    		!store cell size
    		delta(1) = vofb%dxf(iu)
    		delta(2) = vofb%dyf(ju)
    		delta(3) = vofb%dzf(ku)
    				
    		!select normal component
    		m(1) = vofb%nx(iu,ju,ku)
    		m(2) = vofb%ny(iu,ju,ku)
    		m(3) = vofb%nz(iu,ju,ku)  
    					
    		!select q
    		q = vofb%q(iu,ju,ku) 
    					
    		!protruded portion
			delta(dir)=-vel*rdt 
			call vCell(vofb,iu,ju,ku,m,q,delta,c)
    					
    	end if
    	
    	V = vofb%dxf(i)*vofb%dyf(j)*vofb%dzf(k)
    	
    	flux = c*delta(1)*delta(2)*delta(3)/V
    	

    end subroutine
!========================================================================================!

!========================================================================================!
    subroutine updateState(this,vofb)
    	type(VOF), intent(in) :: this
        type(vofBlock), intent(inout) :: vofb
        integer :: i, j, k
        integer :: lbi, ubi, lbj, ubj, lbk, ubk
        
		lbi = lbound(vofb%c,1)
		ubi = ubound(vofb%c,1)
		lbj = lbound(vofb%c,2)
		ubj = ubound(vofb%c,2)
		lbk = lbound(vofb%c,3)
		ubk = ubound(vofb%c,3)
        

		do k=lbk,ubk
			do j=lbj,ubj
				do i=lbi,ubi
					
					!set mixed logical
					if ((vofb%c(i,j,k) >= 0.d0+this%eps_) & 
					   .AND. (vofb%c(i,j,k) <= 1.d0-this%eps_)) then
					   vofb%isMixed(i,j,k) = .TRUE.
					else
						vofb%isMixed(i,j,k) = .FALSE.
					end if
					
					!set full logical
					if (vofb%c(i,j,k) >= 1.d0-this%eps_) then
						vofb%isFull(i,j,k) = .TRUE.
					else
						vofb%isFull(i,j,k) = .FALSE.
					end if
					
				end do
			end do
		end do

    end subroutine
!========================================================================================!

!========================================================================================!
    subroutine updateStateBlocks(this)
    	type(VOF), intent(in) :: this
        integer :: i, j, k, b
        integer :: lbi, ubi, lbj, ubj, lbk, ubk
        
        
        do b=1,s_nblk
        
			lbi = lbound(vofBlocks(b)%c,1)
			ubi = ubound(vofBlocks(b)%c,1)
			lbj = lbound(vofBlocks(b)%c,2)
			ubj = ubound(vofBlocks(b)%c,2)
			lbk = lbound(vofBlocks(b)%c,3)
			ubk = ubound(vofBlocks(b)%c,3)
        

			do k=lbk,ubk
				do j=lbj,ubj
					do i=lbi,ubi
					
						!set mixed logical
						if ((vofBlocks(b)%c(i,j,k) >= 0.d0+this%eps_) & 
					   		.AND. (vofBlocks(b)%c(i,j,k) <= 1.d0-this%eps_)) then
					   		vofBlocks(b)%isMixed(i,j,k) = .TRUE.
						else
							vofBlocks(b)%isMixed(i,j,k) = .FALSE.
						end if
					
						!set full logical
						if (vofBlocks(b)%c(i,j,k) >= 1.d0-this%eps_) then
							vofBlocks(b)%isFull(i,j,k) = .TRUE.
						else
							vofBlocks(b)%isFull(i,j,k) = .FALSE.
						end if
					
					end do
				end do
			end do
		
		end do

    end subroutine
!========================================================================================!

!========================================================================================!
    subroutine reSetFullEmpty(vofb)
    	type(vofBlock), intent(inout) :: vofb
    	integer :: lbi, ubi, lbj, ubj, lbk, ubk
    	integer :: i,j,k
    	
		lbi = lbound(vofb%c,1)
		ubi = ubound(vofb%c,1)
		lbj = lbound(vofb%c,2)
		ubj = ubound(vofb%c,2)
		lbk = lbound(vofb%c,3)
		ubk = ubound(vofb%c,3)

		do k=lbk,ubk
			do j=lbj,ubj
				do i=lbi,ubi
    			
    				if (vofb%isFull(i,j,k)) then
    					vofb%c(i,j,k) = 1.d0
    				else if (.NOT. (vofb%isMixed(i,j,k))) then
    					vofb%c(i,j,k) = 0.d0
    				end if
    			
    			end do
    		end do
    	end do

     end subroutine   	
!========================================================================================!

!========================================================================================!
    subroutine clump(vofb)
    	type(vofBlock), intent(inout) :: vofb
    	integer :: lbi, ubi, lbj, ubj, lbk, ubk
    	integer :: i,j,k
    	
		lbi = lbound(vofb%c,1)
		ubi = ubound(vofb%c,1)
		lbj = lbound(vofb%c,2)
		ubj = ubound(vofb%c,2)
		lbk = lbound(vofb%c,3)
		ubk = ubound(vofb%c,3)

		do k=lbk,ubk
			do j=lbj,ubj
				do i=lbi,ubi
    			
    				if (vofb%c(i,j,k) > 1.d0) then
    					vofb%c(i,j,k) = 1.d0
    				else if (vofb%c(i,j,k) < 0.d0) then
    					vofb%c(i,j,k) = 0.d0
    				end if
    			
    			end do
    		end do
    	end do

     end subroutine   	
!========================================================================================!

!========================================================================================!
    subroutine sweepCombination()

		s_sweepSelector = s_sweepSelector+1
		
		if (s_sweepSelector > 6) then
			s_sweepSelector = 1
		end if
		
		select case(s_sweepSelector)
			case(1)
				s_sweep(1)=s_Xsweep
				s_sweep(2)=s_Ysweep
				s_sweep(3)=s_Zsweep
			case(2)
				s_sweep(1)=s_Ysweep
				s_sweep(2)=s_Zsweep
				s_sweep(3)=s_Xsweep
			case(3)
				s_sweep(1)=s_Zsweep
				s_sweep(2)=s_Xsweep
				s_sweep(3)=s_Ysweep
			case(4)
				s_sweep(1)=s_Xsweep
				s_sweep(2)=s_Zsweep
				s_sweep(3)=s_Ysweep
			case(5)
				s_sweep(1)=s_Ysweep
				s_sweep(2)=s_Xsweep
				s_sweep(3)=s_Zsweep
			case(6)
				s_sweep(1)=s_Zsweep
				s_sweep(2)=s_Ysweep
				s_sweep(3)=s_Xsweep
			case default
		end select
	
    end subroutine
!========================================================================================!

!========================================================================================!
    subroutine info(this,cpuTime,sw)
    	type(VOF), intent(in) :: this
    	real(DP), intent(in) :: cpuTime
    	integer, intent(in) :: sw
    	type(mpiControl), pointer :: comm
    	real(DP) :: cpuTime_max
    	integer :: ierror
    	
    	comm => this%mesh_%ptrMPIC_
    	
    	call Mpi_Reduce(cpuTime, cpuTime_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, &
    			        comm%cartComm_, ierror)

		if (IS_MASTER) then
			if (sw==0) then
				write(*,'(A,'//s_outputFormat(2:9)//')') '	VOF Eqn: CPU time = ', cpuTime_max
			else
				write(*,'(A,'//s_outputFormat(2:9)//')') '	St upd: CPU TIME = ', cpuTime_max
			end if
		end if
		
    	
    end subroutine
!========================================================================================!

!********************************* Blocks management ************************************!

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
    subroutine readSingleBubble(mesh,dict)
    	type(grid), intent(in) :: mesh
		type(dictionary), intent(in) :: dict
		integer :: nref
		real(DP) :: x0,y0,z0,R
		
		call readParameter(dict,nref,'nref')
		call readParameter(dict,x0,'x0')
		call readParameter(dict,y0,'y0')
		call readParameter(dict,z0,'z0')
		call readParameter(dict,R,'R')
					
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
    subroutine readTwoBubbles(mesh,dict)
    	type(grid), intent(in) :: mesh
		type(dictionary), intent(in) :: dict
		integer :: nref
		real(DP) :: x0_b1,y0_b1,z0_b1,R_b1
		real(DP) :: x0_b2,y0_b2,z0_b2,R_b2
		
		call readParameter(dict,nref,'nref')
		!b1
		call readParameter(dict,x0_b1,'x0_b1')
		call readParameter(dict,y0_b1,'y0_b1')
		call readParameter(dict,z0_b1,'z0_b1')
		call readParameter(dict,R_b1,'R_b1')
		!b2
		call readParameter(dict,x0_b2,'x0_b2')
		call readParameter(dict,y0_b2,'y0_b2')
		call readParameter(dict,z0_b2,'z0_b2')
		call readParameter(dict,R_b2,'R_b2')
					
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
    subroutine readBubblesArray(mesh,dict)
    	type(grid), intent(in) :: mesh
		type(dictionary), intent(in) :: dict
		integer :: nbx,nby,nbz,nref
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
		
		call readParameter(dict,nbx,'nbx')
		call readParameter(dict,nby,'nby')
		call readParameter(dict,nbz,'nbz')
		call readParameter(dict,nref,'nref')
		call readParameter(dict,R,'R')
		call readParameter(dict,random_distr,'random_distr')

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
						x0 = x0 + (rnd(1)-0.5d0)*0.85d0*sx
						y0 = y0 + (rnd(2)-0.5d0)*0.85d0*sy
						z0 = z0 + (rnd(3)-0.5d0)*0.85d0*sz
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
		xfv=mesh%xc_(is-1:ie)
		yfv=mesh%yc_(js-1:je)
		zfv=mesh%zc_(ks-1:ke)
		
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

!========================================================================================!
    subroutine initVOFblocks(mesh,gmesh)
    	type(grid), intent(in) :: mesh,gmesh
    	type(mpiControl), pointer :: mpiCTRL
    	integer, parameter :: single_bubbles = 1, array_bubbles = 2, two_bubbles=3
    	real(DP), allocatable, dimension(:,:,:) :: cblk
    	logical, allocatable, dimension(:) :: b_proc_bool
    	type(dictionary) :: dict
    	integer :: method,nref,ierror,bi,is,js,ks,ie,je,ke,nb_tmp,n_blk_max
    	real(DP) :: x0,y0,z0,R
    	logical :: present
    	
    	
    	mpiCTRL => mesh%ptrMPIC_
    	
		if (IS_MASTER) then
			write(*,'(A)') 'INIT VOF BLOCKS'
		end if
			

		if (IS_MASTER) then		
		
			call dictionaryCTOR(dict,'initBubbles','specs')
			call readParameter(dict,method,'method')
			!number of refinements
			call readParameter(dict,nref,'nref')
		
			select case(method)
				case(single_bubbles)
					call readSingleBubble(gmesh,dict)
				case(array_bubbles)
					call readBubblesArray(gmesh,dict)
				case(two_bubbles)
					call readTwoBubbles(gmesh,dict)
				case default
			end select
		
		end if
    	
    	!broadcast total number of bubbles
    	call MPI_BCAST(s_nb, 1, MPI_INTEGER, 0, mpiCTRL%cartComm_, ierror)
    	!broadcast number of refinements
    	call MPI_BCAST(nref, 1, MPI_INTEGER, 0, mpiCTRL%cartComm_, ierror)
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
    	
    	!init block
    	do bi=1,s_nb
    	
    		present=b_proc_bool(bi)
    		
    		if (present) then
    			
    			x0=s_pos_init(1,bi)
    			y0=s_pos_init(2,bi)
    			z0=s_pos_init(3,bi)
    			R=s_pos_init(4,bi)
    			
				is=s_idx_init(1,bi)
				ie=s_idx_init(2,bi)
				js=s_idx_init(3,bi)
				je=s_idx_init(4,bi)
				ks=s_idx_init(5,bi)
				ke=s_idx_init(6,bi)
    			
    			call reAllocateArray(cblk,is,ie,js,je,ks,ke)
    			call init_Bubble_vf(gmesh,cblk,x0,y0,z0,R,nref)
    			
				call addNewBlock(vofBlocks,nb_tmp)
				
				vofBlocks(nb_tmp)%bn=bi
				vofBlocks(nb_tmp)%master=mpiCTRL%rank_
				!indexes
				vofBlocks(nb_tmp)%idx=s_idx_init(:,bi)
				call allocateArray(vofBlocks(nb_tmp)%idx_mx,is-offset_c,ie+offset_c)
				call allocateArray(vofBlocks(nb_tmp)%idx_my,js-offset_c,je+offset_c)
				call allocateArray(vofBlocks(nb_tmp)%idx_mz,ks-offset_c,ke+offset_c)
				!mesh
				call allocateArray(vofBlocks(nb_tmp)%xc,is-offset_c,ie+offset_c)
				call allocateArray(vofBlocks(nb_tmp)%yc,js-offset_c,je+offset_c)
				call allocateArray(vofBlocks(nb_tmp)%zc,ks-offset_c,ke+offset_c)
				call allocateArray(vofBlocks(nb_tmp)%xf,is-offset_c-1,ie+offset_c)
				call allocateArray(vofBlocks(nb_tmp)%yf,js-offset_c-1,je+offset_c)
				call allocateArray(vofBlocks(nb_tmp)%zf,ks-offset_c-1,ke+offset_c)
				call allocateArray(vofBlocks(nb_tmp)%dxc,is-offset_c+1,ie+offset_c)
				call allocateArray(vofBlocks(nb_tmp)%dyc,js-offset_c+1,je+offset_c)
				call allocateArray(vofBlocks(nb_tmp)%dzc,ks-offset_c+1,ke+offset_c)
				call allocateArray(vofBlocks(nb_tmp)%dxf,is-offset_c,ie+offset_c)
				call allocateArray(vofBlocks(nb_tmp)%dyf,js-offset_c,je+offset_c)
				call allocateArray(vofBlocks(nb_tmp)%dzf,ks-offset_c,ke+offset_c)
				!fields
				call allocateArray(vofBlocks(nb_tmp)%c,is-offset_c,ie+offset_c,	&
												       js-offset_c,je+offset_c,	&
												       ks-offset_c,ke+offset_c)
				call allocateArray(vofBlocks(nb_tmp)%c0,is-offset_c,ie+offset_c,	&
												        js-offset_c,je+offset_c,	&
												        ks-offset_c,ke+offset_c)
				call allocateArray(vofBlocks(nb_tmp)%ux,is-offset_u-1,ie+offset_u,	&
												        js-offset_u,je+offset_u,	&
												        ks-offset_u,ke+offset_u)
				call allocateArray(vofBlocks(nb_tmp)%uy,is-offset_u,ie+offset_u,	&
												        js-offset_u-1,je+offset_u,	&
												        ks-offset_u,ke+offset_u)
				call allocateArray(vofBlocks(nb_tmp)%uz,is-offset_u,ie+offset_u,	&
												        js-offset_u,je+offset_u,	&
												        ks-offset_u-1,ke+offset_u)
				call allocateArray(vofBlocks(nb_tmp)%stx,is-offset_st-1,ie+offset_st,	&
												         js-offset_st,je+offset_st,	&
												         ks-offset_st,ke+offset_st)
				call allocateArray(vofBlocks(nb_tmp)%sty,is-offset_st,ie+offset_st,	&
												         js-offset_st-1,je+offset_st,	&
												         ks-offset_st,ke+offset_st)
				call allocateArray(vofBlocks(nb_tmp)%stz,is-offset_st,ie+offset_st,	&
												         js-offset_st,je+offset_st,	&
												         ks-offset_st-1,ke+offset_st)		
				call allocateArray(vofBlocks(nb_tmp)%nx,is-1,ie+1,js-1,je+1,ks-1,ke+1)
				call allocateArray(vofBlocks(nb_tmp)%ny,is-1,ie+1,js-1,je+1,ks-1,ke+1)
				call allocateArray(vofBlocks(nb_tmp)%nz,is-1,ie+1,js-1,je+1,ks-1,ke+1)
				call allocateArray(vofBlocks(nb_tmp)%isMixed,is-offset_c,ie+offset_c,	&
												             js-offset_c,je+offset_c,	&
												             ks-offset_c,ke+offset_c)
				call allocateArray(vofBlocks(nb_tmp)%isFull,is-offset_c,ie+offset_c,	&
												            js-offset_c,je+offset_c,	&
												            ks-offset_c,ke+offset_c)
				call allocateArray(vofBlocks(nb_tmp)%cv,is-offset_u-1,ie+offset_u,	&
												        js-offset_u-1,je+offset_u,	&
												        ks-offset_u-1,ke+offset_u)
				call allocateArray(vofBlocks(nb_tmp)%cs,is-offset_u,ie+offset_u,	&
														js-offset_u,je+offset_u,	&
														ks-offset_u,ke+offset_u)
				call allocateArray(vofBlocks(nb_tmp)%q,is-1,ie+1,js-1,je+1,ks-1,ke+1)
				call allocateArray(vofBlocks(nb_tmp)%k,is-1,ie+1,js-1,je+1,ks-1,ke+1)
				call allocateArray(vofBlocks(nb_tmp)%geoFlux,is-offset_u,ie+offset_u,	&
															 js-offset_u,je+offset_u,	&
															 ks-offset_u,ke+offset_u)
				call allocateArray(vofBlocks(nb_tmp)%corrFlux,is-offset_u,ie+offset_u,	&
															 js-offset_u,je+offset_u,	&
															 ks-offset_u,ke+offset_u)
				call allocateArray(vofBlocks(nb_tmp)%corrTerm,is-offset_u,ie+offset_u,	&
															 js-offset_u,je+offset_u,	&
															 ks-offset_u,ke+offset_u)
								
				!copy c
				vofBlocks(nb_tmp)%c=0.d0
				vofBlocks(nb_tmp)%c(is:ie,js:je,ks:ke)=cblk
				call computeModuloIdx(mesh,vofBlocks(nb_tmp),is-offset_c,ie+offset_c,&
							   			   				     js-offset_c,je+offset_c,&
							   			   				     ks-offset_c,ke+offset_c)
				call copyBlockMesh(mesh,gmesh,vofBlocks(nb_tmp))
    			
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
		
		!check max blocks per MPI_Proc
		call Mpi_Reduce(s_nblk, n_blk_max, 1, MPI_INTEGER, MPI_MAX, 0, &
					    mpiCTRL%cartComm_, ierror)
		
		if (IS_MASTER) then
			write(*,'(A)') 'END INIT VOF BLOCKS'
			write(*,*) 'Max number of blocks per MPI proc: ', n_blk_max
		end if
			
				
    	
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
    	    		
	end subroutine
!========================================================================================!

!========================================================================================!
    subroutine packSendToSlaveBuff(idx_buff,val_buff,vofb,ftype,i0g,i1g,j0g,j1g,k0g,k1g)
    	integer, allocatable, dimension(:,:), intent(inout) :: idx_buff
    	real(DP), allocatable, dimension(:), intent(inout) :: val_buff
    	type(vofBlock), intent(in) :: vofb
    	integer, intent(in) :: i0g,i1g,j0g,j1g,k0g,k1g,ftype
    	integer :: i,j,k,im,jm,km,il,jl,kl,n  	
    
    
    	idx_buff(1,:)=-1
    	n=1
    	
    	select case(ftype)
    		case(1)
    	
				do k=vofb%idx(5),vofb%idx(6)
					do j=vofb%idx(3),vofb%idx(4)
						do i=vofb%idx(1),vofb%idx(2)

							im = vofb%idx_mx(i)
							jm = vofb%idx_my(j)
							km = vofb%idx_mz(k)
					
							il=im-i0g+1
							jl=jm-j0g+1
							kl=km-k0g+1
					
							if ((im<=i1g).AND.(im>=i0g)) then
								if ((jm<=j1g).AND.(jm>=j0g)) then
									if ((km<=k1g).AND.(km>=k0g)) then
										idx_buff(1,n)=0
										idx_buff(2,n)=il
										idx_buff(3,n)=jl
										idx_buff(4,n)=kl
										val_buff(n)=vofb%c(i,j,k)
										n=n+1
									end if
								end if
							end if
					
						end do
					end do
				end do
				
			case(2)

				do k=vofb%idx(5),vofb%idx(6)
					do j=vofb%idx(3),vofb%idx(4)
						do i=vofb%idx(1),vofb%idx(2)

							im = vofb%idx_mx(i)
							jm = vofb%idx_my(j)
							km = vofb%idx_mz(k)
					
							il=im-i0g+1
							jl=jm-j0g+1
							kl=km-k0g+1
					
							if ((im<=i1g).AND.(im>=i0g)) then
								if ((jm<=j1g).AND.(jm>=j0g)) then
									if ((km<=k1g).AND.(km>=k0g)) then
										idx_buff(1,n)=0
										idx_buff(2,n)=il
										idx_buff(3,n)=jl
										idx_buff(4,n)=kl
										val_buff(n)=vofb%k(i,j,k)
										n=n+1
									end if
								end if
							end if
					
						end do
					end do
				end do		
				
			case default
		end select			
			
    						
	end subroutine
!========================================================================================!

!========================================================================================!
    subroutine packSendToSlaveBuff_st(idx_buff,val_buff,st,ftype,i0g,i1g,j0g,j1g,k0g,k1g,offset,&
    								  is_org,ie_org,js_org,je_org,ks_org,ke_org,&
    								  nxg,nyg,nzg,wrap_x,wrap_y,wrap_z)
    	integer, allocatable, dimension(:,:), intent(inout) :: idx_buff
    	real(DP), allocatable, dimension(:), intent(inout) :: val_buff
    	real(DP), allocatable, dimension(:,:,:), intent(in) :: st
    	integer, intent(in) :: i0g,i1g,j0g,j1g,k0g,k1g,ftype,offset,nxg,nyg,nzg
    	integer, intent(in) :: is_org,ie_org,js_org,je_org,ks_org,ke_org
    	integer :: i,j,k,i0gm,j0gm,k0gm,i1gm,j1gm,k1gm,il,jl,kl,n,is,js,ks,ie,je,ke
    	integer, allocatable, dimension(:) :: im,jm,km
    	logical, intent(in) :: wrap_x,wrap_y,wrap_z
    
    
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
					
					if ((im(i)<=i1g).AND.(im(i)>=i0gm)) then
						if ((jm(j)<=j1g).AND.(jm(j)>=j0g)) then
							if ((km(k)<=k1g).AND.(km(k)>=k0g)) then
								idx_buff(1,n)=0
								idx_buff(2,n)=il
								idx_buff(3,n)=jl
								idx_buff(4,n)=kl
								val_buff(n)=st(i,j,k)
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
    subroutine unPackRecvFromMasterBuff(idx_buff,val_buff,f,ftype)
    	integer, allocatable, dimension(:,:), intent(in) :: idx_buff
    	real(DP), allocatable, dimension(:), intent(in) :: val_buff
    	real(DP), allocatable, dimension(:,:,:), intent(inout) :: f
    	integer, intent(in) :: ftype
    	integer :: sizeBuff,i,j,k,n
    	real(DP) :: c0,cv
    
    	sizeBuff=size(val_buff)

		select case(ftype)
			case(1)
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
    		case(2)
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
    subroutine unPackRecvFromSlaveBuff(vofb,ftype,idx_buff,val_buff)
    	type(vofBlock), intent(inout) :: vofb
    	integer, allocatable, dimension(:,:), intent(in) :: idx_buff
    	real(DP), allocatable, dimension(:), intent(in) :: val_buff
    	integer, intent(in) :: ftype
    	integer :: i,j,k,n,sizeBuff	
    
    	sizeBuff=size(val_buff)
    
    	select case(ftype)
    		case(1)
    		    do n=1,sizeBuff
    				if (idx_buff(1,n)==0) then
    					i=idx_buff(2,n)
    					j=idx_buff(3,n)
    					k=idx_buff(4,n)
    					vofb%ux(i,j,k)=val_buff(n)
    				end if
    			end do
    		case(2)
    		    do n=1,sizeBuff
    				if (idx_buff(1,n)==0) then
    					i=idx_buff(2,n)
    					j=idx_buff(3,n)
    					k=idx_buff(4,n)
    					vofb%uy(i,j,k)=val_buff(n)
    				end if
    			end do
    		case(3)
    		    do n=1,sizeBuff
    				if (idx_buff(1,n)==0) then
    					i=idx_buff(2,n)
    					j=idx_buff(3,n)
    					k=idx_buff(4,n)
    					vofb%uz(i,j,k)=val_buff(n)
    				end if
    			end do
    		case default
    	end select

    end subroutine   	
!========================================================================================!

!========================================================================================!
    subroutine excVel(mesh,u)
    	type(grid), intent(in) :: mesh
    	type(vectorField), intent(in) :: u
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
    									      		      
    						call unPackRecvFromSlaveBuff(vofBlocks(bl),1,ux_idx,ux_blk)
    						call unPackRecvFromSlaveBuff(vofBlocks(bl),2,uy_idx,uy_blk)
    						call unPackRecvFromSlaveBuff(vofBlocks(bl),3,uz_idx,uz_blk)
  
    						
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
    					
    						call unPackRecvFromSlaveBuff(vofBlocks(bl),1,ux_idx,ux_blk)
    						call unPackRecvFromSlaveBuff(vofBlocks(bl),2,uy_idx,uy_blk)
    						call unPackRecvFromSlaveBuff(vofBlocks(bl),3,uz_idx,uz_blk)   						
    					
    					end if
    				end do
			
    				
    			end if
    		
    		end do
    	end do


    end subroutine   	
!========================================================================================!

!========================================================================================!
    subroutine excVFK(mesh,c,ftype)
    	type(grid), intent(in) :: mesh
    	type(scalarField), intent(inout) :: c
    	integer, intent(in) :: ftype
    	type(mpiControl), pointer :: mpic
		integer :: nprocs,b,master,slave,n,bl
		integer :: i0g,j0g,k0g,i1g,j1g,k1g
		integer :: is,js,ks,ie,je,ke
		integer :: tag,ierror,position,sizeBuff
    	integer, dimension(MPI_STATUS_SIZE) :: status
    	integer, allocatable, dimension(:,:) :: c_idx
    	real(DP), allocatable, dimension(:) :: c_blk
    	character(len=:), allocatable :: buff
    	
    	

		mpic => mesh%ptrMPIC_
		nprocs = mpic%nProcs_
		tag=0
		
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
    					
    						is = vofBlocks(bl)%idx(1)
    						ie = vofBlocks(bl)%idx(2)
    						js = vofBlocks(bl)%idx(3)
    						je = vofBlocks(bl)%idx(4)
    						ks = vofBlocks(bl)%idx(5)
    						ke = vofBlocks(bl)%idx(6)
    					
    						call allocateBlockSendRecvBuffers(c_idx,c_blk,0,is,ie,js,je,ks,ke)
    					
							call packSendToSlaveBuff(c_idx,c_blk,vofBlocks(bl),ftype,i0g,i1g,j0g,j1g,k0g,k1g)
							call unPackRecvFromMasterBuff(c_idx,c_blk,c%f_,ftype)
						
						end if
						
					end do
						
    			
    			!recv from master
    			else if (mpic%rank_==s_gbList(2,n,b)) then
    			
					master = s_gbList(1,n,b)

 					is = s_gbList(3,n,b)
					ie = s_gbList(4,n,b)
					js = s_gbList(5,n,b)
					je = s_gbList(6,n,b)
					ks = s_gbList(7,n,b)
					ke = s_gbList(8,n,b) 						
    				
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
					               
					call unPackRecvFromMasterBuff(c_idx,c_blk,c%f_,ftype)
				 			
    	
    			!send to slaves
    			else if (mpic%rank_==s_gbList(1,n,b)) then
    			
    				slave = s_gbList(2,n,b)

    				call globalIndexesFromAll(mesh,i0g,i1g,j0g,j1g,k0g,k1g,slave,'cl')
    			
    				is = s_gbList(3,n,b)
    				ie = s_gbList(4,n,b)
    				js = s_gbList(5,n,b)
    				je = s_gbList(6,n,b)
    				ks = s_gbList(7,n,b)
    				ke = s_gbList(8,n,b)
    				
    				do bl=1,s_nblk
    					if (vofBlocks(bl)%bn==b) then
							
							call allocateBlockSendRecvBuffers(c_idx,c_blk,0,is,ie,js,je,ks,ke)
							call packSendToSlaveBuff(c_idx,c_blk,vofBlocks(bl),ftype,i0g,i1g,j0g,j1g,k0g,k1g)
							
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
    subroutine excST(mesh,st)
    	type(grid), intent(in) :: mesh
    	type(vectorField), intent(inout) :: st
    	type(mpiControl), pointer :: mpic
		integer :: nprocs,b,master,slave,n,bl
		integer :: i0g,j0g,k0g,i1g,j1g,k1g
		integer :: is,js,ks,ie,je,ke,nxg,nyg,nzg
		integer :: tag,ierror,position
		integer :: sizeBuff_stx,sizeBuff_sty,sizeBuff_stz
    	integer, dimension(3) :: requests
    	integer, dimension(MPI_STATUS_SIZE,3) :: status
    	integer, allocatable, dimension(:,:) :: stx_idx,sty_idx,stz_idx
    	real(DP), allocatable, dimension(:) :: stx_blk,sty_blk,stz_blk
    	character(len=:), allocatable :: buff_stx,buff_sty,buff_stz
    	logical :: wrap_x,wrap_y,wrap_z
    	

		mpic => mesh%ptrMPIC_
		nprocs = mpic%nProcs_
		
		!reset field
		call set2zero_omp(st%ux_%f_)
		call set2zero_omp(st%uy_%f_)
		call set2zero_omp(st%uz_%f_)
		
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
    					
    						is = vofBlocks(bl)%idx(1)-offset_st
    						ie = vofBlocks(bl)%idx(2)+offset_st
    						js = vofBlocks(bl)%idx(3)-offset_st
    						je = vofBlocks(bl)%idx(4)+offset_st
    						ks = vofBlocks(bl)%idx(5)-offset_st
    						ke = vofBlocks(bl)%idx(6)+offset_st
    					
    						
    						call allocateBlockSendRecvBuffers(stx_idx,stx_blk,1,is,ie,js,je,ks,ke)
    						call packSendToSlaveBuff_st(stx_idx,stx_blk,vofBlocks(bl)%stx,1,&
    												    i0g,i1g,j0g,j1g,k0g,k1g,offset_st,&
    								  				    is,ie,js,je,ks,ke,nxg,nyg,nzg,&
    								  				    wrap_x,wrap_y,wrap_z)
											            
    						call allocateBlockSendRecvBuffers(sty_idx,sty_blk,2,is,ie,js,je,ks,ke)
    						call packSendToSlaveBuff_st(sty_idx,sty_blk,vofBlocks(bl)%sty,2,&
    												    i0g,i1g,j0g,j1g,k0g,k1g,offset_st,&
    								  				    is,ie,js,je,ks,ke,nxg,nyg,nzg,&
    								  				    wrap_x,wrap_y,wrap_z)
											            
    						call allocateBlockSendRecvBuffers(stz_idx,stz_blk,3,is,ie,js,je,ks,ke)
    						call packSendToSlaveBuff_st(stz_idx,stz_blk,vofBlocks(bl)%stz,3,&
    												    i0g,i1g,j0g,j1g,k0g,k1g,offset_st,&
    								  				    is,ie,js,je,ks,ke,nxg,nyg,nzg,&
    								  				    wrap_x,wrap_y,wrap_z)
							
							call unPackRecvFromMasterBuff(stx_idx,stx_blk,st%ux_%f_,2)
							call unPackRecvFromMasterBuff(sty_idx,sty_blk,st%uy_%f_,2)
							call unPackRecvFromMasterBuff(stz_idx,stz_blk,st%uz_%f_,2)
							
						end if
						
					end do
						
    			
    			!recv from master
    			else if (mpic%rank_==s_gbList(2,n,b)) then
    			
					master = s_gbList(1,n,b)

 					is = s_gbList(3,n,b)-offset_st
					ie = s_gbList(4,n,b)+offset_st
					js = s_gbList(5,n,b)-offset_st
					je = s_gbList(6,n,b)+offset_st
					ks = s_gbList(7,n,b)-offset_st
					ke = s_gbList(8,n,b)+offset_st 						
    				
    				call allocateBlockSendRecvBuffers(stx_idx,stx_blk,1,is,ie,js,je,ks,ke)
    				call allocateBlockSendRecvBuffers(sty_idx,sty_blk,2,is,ie,js,je,ks,ke)
    				call allocateBlockSendRecvBuffers(stz_idx,stz_blk,3,is,ie,js,je,ks,ke)
    				
    				tag=0
					sizeBuff_stx = integer_size*size(stx_idx) + realDP_size*size(stx_blk)
					call reAllocateArray(buff_stx,sizeBuff_stx)
					call MPI_IRECV(buff_stx,sizeBuff_stx,MPI_CHARACTER, master, &
					               tag, mpic%cartComm_, requests(1), ierror)
					               
    				tag=1
					sizeBuff_sty = integer_size*size(sty_idx) + realDP_size*size(sty_blk)
					call reAllocateArray(buff_sty,sizeBuff_sty)
					call MPI_IRECV(buff_sty,sizeBuff_sty,MPI_CHARACTER, master, &
					               tag, mpic%cartComm_, requests(2), ierror)
					               
    				tag=2
					sizeBuff_stz = integer_size*size(stz_idx) + realDP_size*size(stz_blk)
					call reAllocateArray(buff_stz,sizeBuff_stz)
					call MPI_IRECV(buff_stz,sizeBuff_stz,MPI_CHARACTER, master, &
					               tag, mpic%cartComm_, requests(3), ierror)
					               
					call MPI_WAITALL(3, requests, status, ierror) 
					
					!unpack stx               
					position = 0
					call MPI_UNPACK(buff_stx, sizeBuff_stx, position, stx_idx, size(stx_idx), &
								    MPI_INTEGER, mpic%cartComm_, ierror)
					call MPI_UNPACK(buff_stx, sizeBuff_stx, position, stx_blk, size(stx_blk), &
								    MPI_DOUBLE_PRECISION, mpic%cartComm_, ierror)					

					!unpack sty               
					position = 0
					call MPI_UNPACK(buff_sty, sizeBuff_sty, position, sty_idx, size(sty_idx), &
								    MPI_INTEGER, mpic%cartComm_, ierror)
					call MPI_UNPACK(buff_sty, sizeBuff_sty, position, sty_blk, size(sty_blk), &
								    MPI_DOUBLE_PRECISION, mpic%cartComm_, ierror)	
								    
					!unpack stz               
					position = 0
					call MPI_UNPACK(buff_stz, sizeBuff_stz, position, stz_idx, size(stz_idx), &
								    MPI_INTEGER, mpic%cartComm_, ierror)
					call MPI_UNPACK(buff_stz, sizeBuff_stz, position, stz_blk, size(stz_blk), &
								    MPI_DOUBLE_PRECISION, mpic%cartComm_, ierror)				               
					
					call unPackRecvFromMasterBuff(stx_idx,stx_blk,st%ux_%f_,2)
					call unPackRecvFromMasterBuff(sty_idx,sty_blk,st%uy_%f_,2)
					call unPackRecvFromMasterBuff(stz_idx,stz_blk,st%uz_%f_,2)
					

    			!send to slaves
    			else if (mpic%rank_==s_gbList(1,n,b)) then
    			
    				slave = s_gbList(2,n,b)

    				call globalIndexesFromAll(mesh,i0g,i1g,j0g,j1g,k0g,k1g,slave,'cl')
    			
    				is = s_gbList(3,n,b)-offset_st
    				ie = s_gbList(4,n,b)+offset_st
    				js = s_gbList(5,n,b)-offset_st
    				je = s_gbList(6,n,b)+offset_St
    				ks = s_gbList(7,n,b)-offset_st
    				ke = s_gbList(8,n,b)+offset_st
    				
    				do bl=1,s_nblk
    					if (vofBlocks(bl)%bn==b) then
							
							call allocateBlockSendRecvBuffers(stx_idx,stx_blk,1,is,ie,js,je,ks,ke)
							call allocateBlockSendRecvBuffers(sty_idx,sty_blk,2,is,ie,js,je,ks,ke)
							call allocateBlockSendRecvBuffers(stz_idx,stz_blk,3,is,ie,js,je,ks,ke)
							
    						call packSendToSlaveBuff_st(stx_idx,stx_blk,vofBlocks(bl)%stx,1,&
    												    i0g,i1g,j0g,j1g,k0g,k1g,offset_st,&
    								  				    is,ie,js,je,ks,ke,nxg,nyg,nzg,&
    								  				    wrap_x,wrap_y,wrap_z)
    						call packSendToSlaveBuff_st(sty_idx,sty_blk,vofBlocks(bl)%sty,2,&
    												    i0g,i1g,j0g,j1g,k0g,k1g,offset_st,&
    								  				    is,ie,js,je,ks,ke,nxg,nyg,nzg,&
    								  				    wrap_x,wrap_y,wrap_z)
    						call packSendToSlaveBuff_st(stz_idx,stz_blk,vofBlocks(bl)%stz,3,&
    												    i0g,i1g,j0g,j1g,k0g,k1g,offset_st,&
    								  				    is,ie,js,je,ks,ke,nxg,nyg,nzg,&
    								  				    wrap_x,wrap_y,wrap_z)
							
    					end if
    				end do
    				
    				!pack stx
    				sizeBuff_stx = integer_size*size(stx_idx) + realDP_size*size(stx_blk)	
    				call reAllocateArray(buff_stx,sizeBuff_stx) 
					position = 0
					call MPI_PACK(stx_idx, size(stx_idx), MPI_INTEGER, buff_stx, sizeBuff_stx, position, &
							      mpic%cartComm_, ierror)
					call MPI_PACK(stx_blk, size(stx_blk), MPI_DOUBLE_PRECISION, buff_stx, sizeBuff_stx, &
								  position, mpic%cartComm_, ierror)
					
					!pack sty		  
    				sizeBuff_sty = integer_size*size(sty_idx) + realDP_size*size(sty_blk)	
    				call reAllocateArray(buff_sty,sizeBuff_sty) 
					position = 0
					call MPI_PACK(sty_idx, size(sty_idx), MPI_INTEGER, buff_sty, sizeBuff_sty, position, &
							      mpic%cartComm_, ierror)
					call MPI_PACK(sty_blk, size(sty_blk), MPI_DOUBLE_PRECISION, buff_sty, sizeBuff_sty, &
								  position, mpic%cartComm_, ierror)
					
					!pack stz			  
    				sizeBuff_stz = integer_size*size(stz_idx) + realDP_size*size(stz_blk)	
    				call reAllocateArray(buff_stz,sizeBuff_stz) 
					position = 0
					call MPI_PACK(stz_idx, size(stz_idx), MPI_INTEGER, buff_stz, sizeBuff_stz, position, &
							      mpic%cartComm_, ierror)
					call MPI_PACK(stz_blk, size(stz_blk), MPI_DOUBLE_PRECISION, buff_stz, sizeBuff_stz, &
								  position, mpic%cartComm_, ierror)
		  
					tag = 0	
					call MPI_ISEND(buff_stx, sizeBuff_stx, MPI_CHARACTER, slave, &
						       	   tag, mpic%cartComm_, requests(1), ierror) 
					tag = 1
					call MPI_ISEND(buff_sty, sizeBuff_sty, MPI_CHARACTER, slave, &
						       	   tag, mpic%cartComm_, requests(2), ierror) 
					tag = 2
					call MPI_ISEND(buff_stz, sizeBuff_stz, MPI_CHARACTER, slave, &
						       	   tag, mpic%cartComm_, requests(3), ierror) 
	           
					call MPI_WAITALL(3, requests, status, ierror)		    				   											
    				
    			end if
    		
    		end do
    	end do
    	
    	call updateBoundariesV(st)


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
			
				is=vofBlocks(b)%idx(1)
				ie=vofBlocks(b)%idx(2)
				js=vofBlocks(b)%idx(3)
				je=vofBlocks(b)%idx(4)
				ks=vofBlocks(b)%idx(5)
				ke=vofBlocks(b)%idx(6)
				
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
    subroutine updateBlock(this,vofb,b)
    	type(VOF), intent(in) :: this
    	type(mpiControl), pointer :: mpic
    	type(vofBlock), intent(inout) :: vofb
    	integer, intent(in) :: b
    	integer :: i,j,k
    	integer :: is,ie,js,je,ks,ke
    	integer :: imin,imax,jmin,jmax,kmin,kmax
    	integer :: imin_blk,imax_blk,jmin_blk,jmax_blk,kmin_blk,kmax_blk
    	
    	mpic => this%mesh_%ptrMPIC_
    	
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
			 call reInitBlock(this,vofb,imin,imax,jmin,jmax,kmin,kmax)
			 
		else
			s_exchange_b(b) = .FALSE.
		end if

    end subroutine   	
!========================================================================================!

!========================================================================================!
    subroutine reInitBlock(this,vofb,is,ie,js,je,ks,ke)
   		type(VOF), intent(in) :: this
    	type(vofBlock), intent(inout) :: vofb
		integer, intent(inout) :: is,ie,js,je,ks,ke
		real(DP), allocatable, dimension(:,:,:) :: tmp_c, tmp_c0
		integer :: nxg,nyg,nzg
		
		nxg = this%mesh_%nxg_
		nyg = this%mesh_%nyg_
		nzg = this%mesh_%nzg_
		
		!tmp copy vof field
		call allocateArray(tmp_c,is,ie,js,je,ks,ke)
		call allocateArray(tmp_c0,is,ie,js,je,ks,ke)
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
		
		vofb%idx(1)=is
		vofb%idx(2)=ie
		vofb%idx(3)=js
		vofb%idx(4)=je
		vofb%idx(5)=ks
		vofb%idx(6)=ke

		call reAllocateArray(vofb%idx_mx,is-offset_c,ie+offset_c)
		call reAllocateArray(vofb%idx_my,js-offset_c,je+offset_c)
		call reAllocateArray(vofb%idx_mz,ks-offset_c,ke+offset_c)	
		
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
				
		call reAllocateArray(vofb%c,is-offset_c,ie+offset_c,	&
								    js-offset_c,je+offset_c,	&
								    ks-offset_c,ke+offset_c)
		call reAllocateArray(vofb%c0,is-offset_c,ie+offset_c,	&
								     js-offset_c,je+offset_c,	&
								     ks-offset_c,ke+offset_c)
		call reAllocateArray(vofb%ux,is-offset_u-1,ie+offset_u,	&
								     js-offset_u,je+offset_u,	&
								     ks-offset_u,ke+offset_u)
		call reAllocateArray(vofb%uy,is-offset_u,ie+offset_u,	&
								     js-offset_u-1,je+offset_u,	&
								     ks-offset_u,ke+offset_u)
		call reAllocateArray(vofb%uz,is-offset_u,ie+offset_u,	&
								     js-offset_u,je+offset_u,	&
								     ks-offset_u-1,ke+offset_u)	
		call reAllocateArray(vofb%stx,is-offset_st-1,ie+offset_st,	&
									  js-offset_st,je+offset_st,	&
									  ks-offset_st,ke+offset_st)
		call reAllocateArray(vofb%sty,is-offset_st,ie+offset_st,	&
									  js-offset_st-1,je+offset_st,	&
									  ks-offset_st,ke+offset_st)
		call reAllocateArray(vofb%stz,is-offset_st,ie+offset_st,	&
									  js-offset_st,je+offset_st,	&
									  ks-offset_st-1,ke+offset_st)	
		call reAllocateArray(vofb%nx,is-1,ie+1,js-1,je+1,ks-1,ke+1)
		call reAllocateArray(vofb%ny,is-1,ie+1,js-1,je+1,ks-1,ke+1)
		call reAllocateArray(vofb%nz,is-1,ie+1,js-1,je+1,ks-1,ke+1)
		call reAllocateArray(vofb%isMixed,is-offset_c,ie+offset_c,	&
										  js-offset_c,je+offset_c,	&
										  ks-offset_c,ke+offset_c)
		call reAllocateArray(vofb%isFull,is-offset_c,ie+offset_c,	&
									     js-offset_c,je+offset_c,	&
									     ks-offset_c,ke+offset_c)
		call reAllocateArray(vofb%cv,is-offset_u-1,ie+offset_u,js-offset_u-1,je+offset_u, &
							         ks-offset_u-1,ke+offset_u)
		call reAllocateArray(vofb%cs,is-offset_u,ie+offset_u,js-offset_u,je+offset_u, &
								     ks-offset_u,ke+offset_u)
		call reAllocateArray(vofb%q,is-1,ie+1,js-1,je+1,ks-1,ke+1)
		call reAllocateArray(vofb%k,is-1,ie+1,js-1,je+1,ks-1,ke+1)
		call reAllocateArray(vofb%geoFlux,is-offset_u,ie+offset_u,js-offset_u,	&
										  je+offset_u,ks-offset_u,ke+offset_u)
		call reAllocateArray(vofb%corrFlux,is-offset_u,ie+offset_u,js-offset_u,	&
										   je+offset_u,ks-offset_u,ke+offset_u)
		call reAllocateArray(vofb%corrTerm,is-offset_u,ie+offset_u,js-offset_u,	&
										   je+offset_u,ks-offset_u,ke+offset_u)
		
		vofb%c=0.d0
		vofb%c(is:ie,js:je,ks:ke) = tmp_c
		vofb%c0=0.d0
		vofb%c0(is:ie,js:je,ks:ke) = tmp_c0
		
		call computeModuloIdx(this%mesh_,vofb,is-offset_c,ie+offset_c,&
								              js-offset_c,je+offset_c,&
								              ks-offset_c,ke+offset_c)
						   
		!copy mesh
		call copyBlockMesh(this%mesh_,this%gmesh_,vofb)
		
		call checkBlockSize(this,vofb)

    end subroutine   	
!========================================================================================!

!========================================================================================!
    subroutine initBlock(this,vofb,bm,bn,c,c0,is,ie,js,je,ks,ke)
   		type(VOF), intent(in) :: this
    	type(vofBlock), intent(inout) :: vofb
    	real(DP), allocatable, dimension(:,:,:) :: c,c0
		integer, intent(inout) :: is,ie,js,je,ks,ke
		integer, intent(in) :: bm,bn
		integer :: nxg,nyg,nzg
		
		nxg = this%mesh_%nxg_
		nyg = this%mesh_%nyg_
		nzg = this%mesh_%nzg_
		
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
		
		!bubble master
		vofb%master=bm
		
		!bubble number
		vofb%bn=bn
		
		vofb%idx(1)=is
		vofb%idx(2)=ie
		vofb%idx(3)=js
		vofb%idx(4)=je
		vofb%idx(5)=ks
		vofb%idx(6)=ke
		
		call allocateArray(vofb%idx_mx,is-offset_c,ie+offset_c)
		call allocateArray(vofb%idx_my,js-offset_c,je+offset_c)
		call allocateArray(vofb%idx_mz,ks-offset_c,ke+offset_c)	
		
		call allocateArray(vofb%xc,is-offset_c,ie+offset_c)
		call allocateArray(vofb%yc,js-offset_c,je+offset_c)
		call allocateArray(vofb%zc,ks-offset_c,ke+offset_c)
		call allocateArray(vofb%xf,is-offset_c-1,ie+offset_c)
		call allocateArray(vofb%yf,js-offset_c-1,je+offset_c)
		call allocateArray(vofb%zf,ks-offset_c-1,ke+offset_c)
		call allocateArray(vofb%dxc,is-offset_c+1,ie+offset_c)
		call allocateArray(vofb%dyc,js-offset_c+1,je+offset_c)
		call allocateArray(vofb%dzc,ks-offset_c+1,ke+offset_c)
		call allocateArray(vofb%dxf,is-offset_c,ie+offset_c)
		call allocateArray(vofb%dyf,js-offset_c,je+offset_c)
		call allocateArray(vofb%dzf,ks-offset_c,ke+offset_c)
				
		call allocateArray(vofb%c,is-offset_c,ie+offset_c,	&
								    js-offset_c,je+offset_c,	&
								    ks-offset_c,ke+offset_c)
		call allocateArray(vofb%c0,is-offset_c,ie+offset_c,	&
								    js-offset_c,je+offset_c,	&
								    ks-offset_c,ke+offset_c)
		call allocateArray(vofb%ux,is-offset_u-1,ie+offset_u,	&
								     js-offset_u,je+offset_u,	&
								     ks-offset_u,ke+offset_u)
		call allocateArray(vofb%uy,is-offset_u,ie+offset_u,	&
								     js-offset_u-1,je+offset_u,	&
								     ks-offset_u,ke+offset_u)
		call allocateArray(vofb%uz,is-offset_u,ie+offset_u,	&
								     js-offset_u,je+offset_u,	&
								     ks-offset_u-1,ke+offset_u)	
		call allocateArray(vofb%stx,is-offset_st-1,ie+offset_st,	&
									js-offset_st,je+offset_st,	&
									ks-offset_st,ke+offset_st)
		call allocateArray(vofb%sty,is-offset_st,ie+offset_st,	&
									js-offset_st-1,je+offset_st,	&
									ks-offset_st,ke+offset_st)
		call allocateArray(vofb%stz,is-offset_st,ie+offset_st,	&
									js-offset_st,je+offset_st,	&
									ks-offset_st-1,ke+offset_st)	
		call allocateArray(vofb%nx,is-1,ie+1,js-1,je+1,ks-1,ke+1)
		call allocateArray(vofb%ny,is-1,ie+1,js-1,je+1,ks-1,ke+1)
		call allocateArray(vofb%nz,is-1,ie+1,js-1,je+1,ks-1,ke+1)
		call allocateArray(vofb%isMixed,is-offset_c,ie+offset_c,	&
										  js-offset_c,je+offset_c,	&
										  ks-offset_c,ke+offset_c)
		call allocateArray(vofb%isFull,is-offset_c,ie+offset_c,	&
									     js-offset_c,je+offset_c,	&
									     ks-offset_c,ke+offset_c)
		call allocateArray(vofb%cv,is-offset_u-1,ie+offset_u,js-offset_u-1,je+offset_u, &
							         ks-offset_u-1,ke+offset_u)
		call allocateArray(vofb%cs,is-offset_u,ie+offset_u,js-offset_u,je+offset_u, &
								     ks-offset_u,ke+offset_u)
		call allocateArray(vofb%q,is-1,ie+1,js-1,je+1,ks-1,ke+1)
		call allocateArray(vofb%k,is-1,ie+1,js-1,je+1,ks-1,ke+1)
		call allocateArray(vofb%geoFlux,is-offset_u,ie+offset_u,js-offset_u,	&
										  je+offset_u,ks-offset_u,ke+offset_u)
		call allocateArray(vofb%corrFlux,is-offset_u,ie+offset_u,js-offset_u,	&
										   je+offset_u,ks-offset_u,ke+offset_u)
		call allocateArray(vofb%corrTerm,is-offset_u,ie+offset_u,js-offset_u,	&
										   je+offset_u,ks-offset_u,ke+offset_u)
		
		vofb%c=0.d0
		vofb%c(is:ie,js:je,ks:ke) = c
		vofb%c0=0.d0
		vofb%c0(is:ie,js:je,ks:ke) = c0
		
		call computeModuloIdx(this%mesh_,vofb,is-offset_c,ie+offset_c,&
								        js-offset_c,je+offset_c,&
								        ks-offset_c,ke+offset_c)
						   
		!copy mesh
		call copyBlockMesh(this%mesh_,this%gmesh_,vofb)
		
		call checkBlockSize(this,vofb)

    end subroutine   	
!========================================================================================!

!========================================================================================!
    subroutine reInitBlockDistribution(this)
   		type(VOF), intent(in) :: this
		real(DP), allocatable, dimension(:,:,:) :: c_blk,c0_blk
		integer :: b,n,nb,nprocs,bl,is,ie,js,je,ks,ke,i0g,i1g,j0g,j1g,k0g,k1g,nxg,nyg,nzg
		integer :: master,slave,tag,ierror,sizeBuff,position
		integer, dimension(MPI_STATUS_SIZE) :: status
		logical :: isHere
		type(vofBlock), allocatable, dimension(:) :: new_blocks
		type(mpiControl), pointer :: mpic
		character(len=:), allocatable :: buff
    		

		mpic => this%mesh_%ptrMPIC_
		nprocs = mpic%nProcs_
		tag=0
		
		nxg = this%mesh_%nxg_
		nyg = this%mesh_%nyg_
		nzg = this%mesh_%nzg_

    	do b=1,s_nb
    	  	
    		do n=0,nprocs-1
    		
    			!manage own bubbls
    			if ((mpic%rank_==s_gbList(1,n,b)).AND.(s_gbList(1,n,b)==s_gbList(2,n,b))) then
    			
    				i0g = this%mesh_%i0g_
    				i1g = this%mesh_%i1g_
    				j0g = this%mesh_%j0g_
    				j1g = this%mesh_%j1g_
    				k0g = this%mesh_%k0g_
    				k1g = this%mesh_%k1g_
    				
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
					
    				i0g = this%mesh_%i0g_
    				i1g = this%mesh_%i1g_
    				j0g = this%mesh_%j0g_
    				j1g = this%mesh_%j1g_
    				k0g = this%mesh_%k0g_
    				k1g = this%mesh_%k1g_

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
    				    
    				    !ricevi c,c0
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
    					call initBlock(this,new_blocks(nb),slave,b,c_blk,c0_blk, &
    								   is,ie,js,je,ks,ke)	
    					
    				end if      		    
    				
				
				!this proc is the master
				else if (mpic%rank_==s_gbList(1,n,b)) then
				
					slave = s_gbList(2,n,b)

					call globalIndexesFromAll(this%mesh_,i0g,i1g,j0g,j1g,k0g,k1g,slave,'cl')	

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

		!update blocks state
		call updateStateBlocks(this)
		
		!reconstruct for k computation
		do b=1,s_nblk
			call reconstruct(this,vofBlocks(b))
		end do
		
		!update exchange list
		call reAllocateArray(s_blk_data,1,8,1,s_nblk)
		call reAllocateArray(s_blk_proc,1,s_nblk,0,nprocs-1)
		call excLists(this%mesh_)
		
		!reinit exchange logical
		if (s_nblk>0) then
			call reAllocateArray(s_exchange_b,1,s_nblk)
		else
			call reAllocateArray(s_exchange_b,1,1)
		end if
		s_exchange_g(mpic%rank_)=.FALSE.
		s_exchange_b = .FALSE.


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
		call assignArray(lhs%cs,rhs%cs)
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
    subroutine resetFragments(vofb)
    	type(vofBlock), intent(inout) :: vofb
    	integer :: is,ie,js,je,ks,ke
		integer :: i,j,k,ii,jj,kk
		logical :: isFrag

		
		is=vofb%idx(1)
		ie=vofb%idx(2)
		js=vofb%idx(3)
		je=vofb%idx(4)
		ks=vofb%idx(5)
		ke=vofb%idx(6)		
    	
		do k=ks,ke
			do j=js,je
				do i=is,ie
					
					if (vofb%isMixed(i,j,k)) then
						isFrag=.TRUE.
						!search stencil
						do kk=-2,2
							do jj=-2,2
								do ii=-2,2	
								
									if (vofb%isFull(i+ii,j+jj,k+kk)) then
										isFrag=.FALSE.
									end if
									
								end do
							end do
						end do	
						
						if (isFrag) then
							vofb%c(i,j,k)=0.d0
							vofb%isMixed(i,j,k)=.FALSE.
						end if
					
					end if		
					
				end do
			end do
		end do
		
    end subroutine   	
!========================================================================================!

!========================================================================================!
    subroutine checkBlockSize(this,vofb)
    	type(VOF), intent(in) :: this
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
    	
		nxg = this%mesh_%nxg_
		nyg = this%mesh_%nyg_
		nzg = this%mesh_%nzg_
		
		wrap_x = this%mesh_%ptrMPIC_%wrapAround_(1)
		wrap_y = this%mesh_%ptrMPIC_%wrapAround_(2)
		wrap_z = this%mesh_%ptrMPIC_%wrapAround_(3)
	
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
    subroutine gatherLogicalExchange(this)
    	type(VOF), intent(in) :: this
    	type(mpiControl), pointer :: mpic
    	logical :: s_exchange
    	integer :: ierror
    	
    	mpic => this%mesh_%ptrMPIC_

		s_exchange = any(s_exchange_b)

		call MPI_ALLGather(s_exchange,1,MPI_LOGICAL,s_exchange_g,1,MPI_LOGICAL,&
						   mpic%cartComm_,ierror)
		
     end subroutine   	
!========================================================================================!

!========================================================================================!
    subroutine printVOFblocks(nFolder)
    	integer, intent(in) :: nFolder
        character(len=20) :: dirName,bn
        integer :: i,lbi,ubi,lbj,ubj,lbk,ubk
			
        
        write(dirName,s_intFormat) nFolder
        
		do i=1,s_nblk
		
			lbi=vofBlocks(i)%idx(1)
			ubi=vofBlocks(i)%idx(2)
			lbj=vofBlocks(i)%idx(3)
			ubj=vofBlocks(i)%idx(4)
			lbk=vofBlocks(i)%idx(5)
			ubk=vofBlocks(i)%idx(6)
						
			write(bn,s_intFormat) vofBlocks(i)%bn
        	open(UNIT=s_IOunitNumber,FILE=trim(adjustl(dirName))//'/b'//trim(adjustl(bn)), &
        	 	      form='UNFORMATTED',ACCESS='STREAM',STATUS='REPLACE',ACTION='WRITE')
        
        		
        		write(s_IOunitNumber) lbi
        		write(s_IOunitNumber) ubi
        		write(s_IOunitNumber) lbj
        		write(s_IOunitNumber) ubj
        		write(s_IOunitNumber) lbk
        		write(s_IOunitNumber) ubk
        		
        		write(s_IOunitNumber) vofBlocks(i)%c
        		write(s_IOunitNumber) vofBlocks(i)%k
        		
        		write(s_IOunitNumber) vofBlocks(i)%ux
        		write(s_IOunitNumber) vofBlocks(i)%uy
        		write(s_IOunitNumber) vofBlocks(i)%uz
        		
        		!write(s_IOunitNumber) vofBlocks(i)%stx
        		!write(s_IOunitNumber) vofBlocks(i)%sty
        		!write(s_IOunitNumber) vofBlocks(i)%stz

			
			close(s_IOunitNumber)			
		end do
        
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

!******************************* Material properties ***********************************!

!========================================================================================!
    subroutine updateMaterialProps(this,c,cs,rho,mu)
    	type(VOF), intent(in) :: this
    	type(scalarField), intent(in) :: c
    	type(scalarField), intent(inout) :: cs,rho,mu
        integer :: lbi, ubi, lbj, ubj, lbk, ubk
        integer :: i,j,k
       
        !all the mat props have the same dimension 
		lbi = lbound(rho%f_,1)
		ubi = ubound(rho%f_,1)
		lbj = lbound(rho%f_,2)
		ubj = ubound(rho%f_,2)
		lbk = lbound(rho%f_,3)
		ubk = ubound(rho%f_,3)
		
		call smoothVF(this,c,cs)

		!$OMP PARALLEL DO DEFAULT(none) &
		!$OMP SHARED(mu,rho,cs,this) &
		!$OMP SHARED(lbk,ubk,lbj,ubj,lbi,ubi) &
		!$OMP PRIVATE(i,j,k) 
		do k=lbk,ubk
			do j=lbj,ubj
				do i=lbi,ubi
					
					mu%f_(i,j,k) = this%mug_*cs%f_(i,j,k)+(1.d0-cs%f_(i,j,k))*this%mul_
					rho%f_(i,j,k) = this%rhog_*cs%f_(i,j,k)+(1.d0-cs%f_(i,j,k))*this%rhol_
					
				end do
			end do
		end do   
		!$OMP END PARALLEL DO
		
	
    end subroutine
!========================================================================================!

!******************************* Normal and Curvature ***********************************!

!========================================================================================!
    subroutine computeNormal_sharp(this,vofb)
    	type(VOF), intent(in) :: this
        type(vofBlock), intent(inout) :: vofb
        integer :: i, j, k
        integer :: lbi, ubi, lbj, ubj, lbk, ubk
        
		lbi = lbound(vofb%nx,1)
		ubi = ubound(vofb%nx,1)
		lbj = lbound(vofb%nx,2)
		ubj = ubound(vofb%nx,2)
		lbk = lbound(vofb%nx,3)
		ubk = ubound(vofb%nx,3)

		do k=lbk,ubk
			do j=lbj,ubj
				do i=lbi,ubi
					
					vofb%nx(i,j,k) = -(vofb%c(i+1,j,k)-vofb%c(i-1,j,k))/ &
										(vofb%dxc(i+1)+vofb%dxc(i))
					vofb%ny(i,j,k) = -(vofb%c(i,j+1,k)-vofb%c(i,j-1,k))/ &
										(vofb%dyc(j+1)+vofb%dyc(j))
					vofb%nz(i,j,k) = -(vofb%c(i,j,k+1)-vofb%c(i,j,k-1))/ &
										(vofb%dzc(k+1)+vofb%dzc(k))
					
				end do
			end do
		end do    
		
		!call correctContantAngle(this,vofb)
		
    end subroutine
!========================================================================================!

!========================================================================================!
    subroutine computeNormal_smooth(this,vofb)
    	type(VOF), intent(in) :: this
        type(vofBlock), intent(inout) :: vofb
        integer :: i, j, k
        integer :: lbi, ubi, lbj, ubj, lbk, ubk
        
		lbi = lbound(vofb%nx,1)
		ubi = ubound(vofb%nx,1)
		lbj = lbound(vofb%nx,2)
		ubj = ubound(vofb%nx,2)
		lbk = lbound(vofb%nx,3)
		ubk = ubound(vofb%nx,3)

		do k=lbk,ubk
			do j=lbj,ubj
				do i=lbi,ubi
					
					vofb%nx(i,j,k) = -(vofb%cs(i+1,j,k)-vofb%cs(i-1,j,k))/ &
										(vofb%dxc(i+1)+vofb%dxc(i))
					vofb%ny(i,j,k) = -(vofb%cs(i,j+1,k)-vofb%cs(i,j-1,k))/ &
										(vofb%dyc(j+1)+vofb%dyc(j))
					vofb%nz(i,j,k) = -(vofb%cs(i,j,k+1)-vofb%cs(i,j,k-1))/ &
										(vofb%dzc(k+1)+vofb%dzc(k))
					
				end do
			end do
		end do    
		
		!call correctContantAngle(this,vofb)
		
    end subroutine
!========================================================================================!

!========================================================================================!
    subroutine computeNormal_youngs(this,vofb)
    	type(VOF), intent(in) :: this
        type(vofBlock), intent(inout) :: vofb
         integer :: i, j, k, q, iv,jv,kv
        integer :: lbi, ubi, lbj, ubj, lbk, ubk
        integer, dimension(3,8) :: ind
        real(DP), dimension(3,8) :: nv
        real(DP) :: cbar1,cbar2
        
		lbi = lbound(vofb%nx,1)
		ubi = ubound(vofb%nx,1)
		lbj = lbound(vofb%nx,2)
		ubj = ubound(vofb%nx,2)
		lbk = lbound(vofb%nx,3)
		ubk = ubound(vofb%nx,3)
		

		do k=lbk,ubk
			do j=lbj,ubj
				do i=lbi,ubi
					
					!v1
					ind(1,1)=i
					ind(2,1)=j
					ind(3,1)=k
					!v2
					ind(1,2)=i+1
					ind(2,2)=j
					ind(3,2)=k
					!v3
					ind(1,3)=i
					ind(2,3)=j
					ind(3,3)=k+1
					!v4								
					ind(1,4)=i+1
					ind(2,4)=j
					ind(3,4)=k+1
					!v5
					ind(1,5)=i
					ind(2,5)=j+1
					ind(3,5)=k+1
					!v6
					ind(1,6)=i+1
					ind(2,6)=j+1
					ind(3,6)=k+1
					!v7
					ind(1,7)=i+1
					ind(2,7)=j+1
					ind(3,7)=k
					!v8								
					ind(1,8)=i
					ind(2,8)=j+1
					ind(3,8)=k
					
			
					do q=1,8
					
						iv=ind(1,q)
						jv=ind(2,q)
						kv=ind(3,q)
						
						cbar1=0.25d0*(vofb%c(iv-1,jv,kv)+vofb%c(iv-1,jv,kv-1)+&
						  			  vofb%c(iv-1,jv-1,kv)+vofb%c(iv-1,jv-1,kv-1))
						cbar2=0.25d0*(vofb%c(iv,jv,kv)+vofb%c(iv,jv,kv-1)+&
						  			  vofb%c(iv,jv-1,kv)+vofb%c(iv,jv-1,kv-1))
						nv(1,q)=(cbar2-cbar1)/vofb%dxc(iv)
					
						cbar1=0.25d0*(vofb%c(iv,jv-1,kv)+vofb%c(iv-1,jv-1,kv)+&
						  			  vofb%c(iv,jv-1,kv-1)+vofb%c(iv-1,jv-1,kv-1))
						cbar2=0.25d0*(vofb%c(iv,jv,kv)+vofb%c(iv-1,jv,kv)+&
						  			  vofb%c(iv,jv,kv-1)+vofb%c(iv-1,jv,kv-1))
						nv(2,q)=(cbar2-cbar1)/vofb%dyc(jv)
					
						cbar1=0.25d0*(vofb%c(iv,jv,kv-1)+vofb%c(iv-1,jv,kv-1)+&
						  			  vofb%c(iv,jv-1,kv-1)+vofb%c(iv-1,jv-1,kv-1))
						cbar2=0.25d0*(vofb%c(iv,jv,kv)+vofb%c(iv-1,jv,kv)+&
						  			  vofb%c(iv,jv-1,kv)+vofb%c(iv-1,jv-1,kv))
						nv(3,q)=(cbar2-cbar1)/vofb%dzc(kv)
						
					end do
					
					vofb%nx(i,j,k)=-0.125d0*sum(nv(1,:))
					vofb%ny(i,j,k)=-0.125d0*sum(nv(2,:))
					vofb%nz(i,j,k)=-0.125d0*sum(nv(3,:))
					
				end do
			end do
		end do  
		
		!call correctContantAngle(this,vofb)
		
    end subroutine
!========================================================================================!

!========================================================================================!
    subroutine computeNormal_hf(this,vofb)
    	type(VOF), intent(in) :: this
        type(vofBlock), intent(inout) :: vofb
        integer :: i, j, k,cn
        integer :: lbi, ubi, lbj, ubj, lbk, ubk
        real(DP) :: mx, my, mz
        real(DP) :: mmax
        integer :: mloc
        real(DP), dimension(3) :: mv
        real(DP), dimension(9) :: h
        integer, dimension(9) :: lb
        logical, dimension(9) :: isValid
        integer, dimension(3,9) :: d
        real(DP), dimension(3,9) :: pos_hf
        logical :: hfk
        
		lbi = lbound(vofb%nx,1)
		ubi = ubound(vofb%nx,1)
		lbj = lbound(vofb%nx,2)
		ubj = ubound(vofb%nx,2)
		lbk = lbound(vofb%nx,3)
		ubk = ubound(vofb%nx,3)
		
		call computeNormal_youngs(this,vofb)

		do k=lbk,ubk
			do j=lbj,ubj
				do i=lbi,ubi

					hfk = .TRUE.
				
					if (vofb%isMixed(i,j,k)) then
					
						mx = vofb%nx(i,j,k)
						my = vofb%ny(i,j,k)
						mz = vofb%nz(i,j,k)
					
						mv = (/ mx, my, mz /)
						mloc = maxloc( abs(mv),1 )
						mmax = mv(mloc)
					
						!chiama procedura colonna
						call setStencilPar(mloc,d)
						do cn=1,9
							call hfColumn(this,vofb,mmax,mloc,i+d(1,cn),j+d(2,cn),k+d(3,cn),cn,h,lb,pos_hf,isValid)
							if (.NOT. (isValid(cn))) then
								hfk = .FALSE.
								!goto 100
							end if						
						end do
				
						!compute normal
						if (hfk) then
							call hfNormal(vofb,mloc,mmax,i,j,k,h,lb)		
						end if	
						
					end if				
					
				end do
			end do
		end do    
		
		!call correctContantAngle(this,vofb)
		
    end subroutine
!========================================================================================!

!========================================================================================!
    subroutine hfNormal(vofb,dir,m,i,j,k,h,lb)
        type(vofBlock), intent(inout) :: vofb
        integer, intent(in) :: dir, i, j, k
        real(DP), intent(in) :: m
        real(DP), dimension(9), intent(inout) :: h
        integer, dimension(9), intent(inout) :: lb
		real(DP) :: dhx, dhy
		integer :: lmax
		integer :: nc,n1,n2,dir1,dir2
		integer :: q
		real(DP), pointer, dimension(:) :: df => NULL()
        real(DP), pointer, dimension(:) :: dc1 => NULL()
        real(DP), pointer, dimension(:) :: dc2 => NULL()
        !$omp threadprivate(df,dc1,dc2)
		
		select case(dir)
			case(1)
				df => vofb%dxf
			case(2)
				df => vofb%dyf
			case(3)
				df => vofb%dzf
			case default
		end select
		
		lmax = maxval( lb,1 )
		
		call setDirections(dir,i,j,k,nc,n1,n2,dir1,dir2)
		
		select case(dir1)
			case(1)
				dc1 => vofb%dxc
			case(2)
				dc1 => vofb%dyc
			case(3)
				dc1 => vofb%dzc
			case default
		end select
		
		select case(dir2)
			case(1)
				dc2 => vofb%dxc
			case(2)
				dc2 => vofb%dyc
			case(3)
				dc2 => vofb%dzc
			case default
		end select
		
		!make reference equal
		if (m < 0.d0) then
			do q=1,9
				h(q) = h(q) + sum(df(nc+lb(q)+1:nc+lmax))
			end do
		else
			do q=1,9
				h(q) = h(q) + sum(df(nc-lmax:nc-lb(q)-1))	
			end do	
		end if
		
		dhx = (h(6) - h(4))/(dc1(n1+1)+dc1(n1))
		dhy = (h(8) - h(2))/(dc2(n2+1)+dc2(n2))
		
		select case(dir)
			case(1)
				vofb%nx(i,j,k)=sign(1.d0,vofb%nx(i,j,k))
				vofb%ny(i,j,k)=-dhx
				vofb%nz(i,j,k)=-dhy
			case(2)
				vofb%nx(i,j,k)=-dhx
				vofb%ny(i,j,k)=sign(1.d0,vofb%ny(i,j,k))
				vofb%nz(i,j,k)=-dhy
			case(3)
				vofb%nx(i,j,k)=-dhx
				vofb%ny(i,j,k)=-dhy
				vofb%nz(i,j,k)=sign(1.d0,vofb%nz(i,j,k))
			case default
		end select		

		
    end subroutine
!========================================================================================!

!========================================================================================!
    subroutine smoothVFblock(vofb)
        type(vofBlock), intent(inout) :: vofb
		integer :: lbi, ubi,lbj, ubj, lbk, ubk
		real(DP), allocatable, dimension(:,:,:) :: ctmp
		integer :: i
		
		lbi = lbound(vofb%c,1)
		ubi = ubound(vofb%c,1)
		lbj = lbound(vofb%c,2)
		ubj = ubound(vofb%c,2)
		lbk = lbound(vofb%c,3)
		ubk = ubound(vofb%c,3)
        
        
        !step 1
        call allocateArray(ctmp,lbi,ubi,lbj,ubj,lbk,ubk)
        ctmp = vofb%c
	    call cellToVertexBlock(ctmp,vofb%cv,vofb%xc,vofb%yc,vofb%zc,	&
	     					   vofb%xf,vofb%yf,vofb%zf)					  
		call vertexToCellBlock(vofb%cs,vofb%cv,vofb%xc,vofb%yc,vofb%zc, &
		 					   vofb%xf,vofb%yf,vofb%zf)
		!step 2	
		
		do i=1,1
		
		ctmp(lbi+1:ubi-1,lbj+1:ubj-1,lbk+1:ubk-1) = vofb%cs				   
	    call cellToVertexBlock(ctmp,vofb%cv,vofb%xc,vofb%yc,vofb%zc,	&
	     					   vofb%xf,vofb%yf,vofb%zf)					  
		call vertexToCellBlock(vofb%cs,vofb%cv,vofb%xc,vofb%yc,vofb%zc, &
		 					   vofb%xf,vofb%yf,vofb%zf)
		 					   
		end do
		 					   
		
    end subroutine
!========================================================================================!

!========================================================================================!
    subroutine smoothVF(this,c,cs)
        type(VOF), intent(in) :: this
        type(scalarField), intent(in) :: c
        type(scalarField), intent(inout) :: cs
        
        !step 1
        call cellToVertex(c%f_,s_cv,							        &	
        				  this%mesh_%xc_,this%mesh_%yc_,this%mesh_%zc_,	&
        				  this%mesh_%xf_,this%mesh_%yf_,this%mesh_%zf_)			  
        call vertexToCell(cs%f_,s_cv,						            &
        				  this%mesh_%xc_,this%mesh_%yc_,this%mesh_%zc_,	&
        				  this%mesh_%xf_,this%mesh_%yf_,this%mesh_%zf_)				  
        call updateBoundaries(cs)
        				  
		!step 2
        call cellToVertex(cs%f_,s_cv,							        &	
        				  this%mesh_%xc_,this%mesh_%yc_,this%mesh_%zc_,	&
        				  this%mesh_%xf_,this%mesh_%yf_,this%mesh_%zf_)				  
        call vertexToCell(cs%f_,s_cv,						            &
        				  this%mesh_%xc_,this%mesh_%yc_,this%mesh_%zc_,	&
        				  this%mesh_%xf_,this%mesh_%yf_,this%mesh_%zf_)
		call updateBoundaries(cs)
		
    end subroutine
!========================================================================================!

!========================================================================================!
    subroutine correctContantAngle(this,vofb)
        type(VOF), intent(in) :: this
        type(vofBlock), intent(inout) :: vofb
        integer :: i, j, k, is, ks
        integer :: lbi, ubi,lbj, ubj, lbk, ubk, ny, jw
        
        ny=this%mesh_%nyg_
        
		lbi = lbound(vofb%nx,1)
		ubi = ubound(vofb%nx,1)
		lbj = lbound(vofb%nx,2)
		ubj = ubound(vofb%nx,2)
		lbk = lbound(vofb%nx,3)
		ubk = ubound(vofb%nx,3)
		
		do j=lbj,ubj
			if (vofb%idx_my(j)==1) then
				jw=j
        		!bottom boundary
				do k=lbk,ubk
					do i=lbi,ubi
				
						if (vofb%isMixed(i,jw,k)) then
							vofb%ny(i,jw,k) = -(vofb%c(i,jw+1,k)-vofb%c(i,jw,k))/ &
										  	   (vofb%dyc(jw+1))					
							do ks=-1,1
							 do is=-1,1
								if (vofb%c(i+is,jw,k+ks)==0.d0) then
									vofb%ny(i,jw,k) = -1.d0	
				 				end if
				 			end do
							end do
						end if
								
					end do
				end do
			end if
		end do
		
 
		
        !top boundary
		do j=lbj,ubj
			if (vofb%idx_my(j)==ny) then
				jw=j
				!top boundary
				do k=lbk,ubk
					do i=lbi,ubi
						if (vofb%isMixed(i,jw,k)) then
							vofb%ny(i,jw,k) = -(vofb%c(i,jw,k)-vofb%c(i,jw-1,k))/ &
									   	  	   (vofb%dyc(jw))					
							do ks=-1,1
							 do is=-1,1
					 			if (vofb%c(i+is,jw,k+ks)==0.d0) then
					 				vofb%ny(i,jw,k) = 1.d0	
					 			end if
					 		end do
							end do
						end if
								
					end do
				end do 
			end if 
		end do 
		
    end subroutine
!========================================================================================!

!========================================================================================!
    subroutine computeSurfaceTension(this,st,k)
    	type(VOF), intent(in) :: this
    	type(vectorField), intent(inout) :: st
    	type(scalarField), intent(inout) :: k
    	integer :: b
    	real(DP) :: t_S, t_E
    	
    	t_S = MPI_Wtime()
    	
		call computeCurvature(this,k)

		!$OMP PARALLEL DO DEFAULT(none) &
		!$OMP SHARED(vofBlocks,this,s_nblk) &
		!$OMP PRIVATE(b)  		
    	do b=1,s_nblk
    		call computeBlockST(this,vofBlocks(b))
    	end do	
    	!$OMP END PARALLEL DO	
		
		call excST(this%mesh_,st)
		
		t_E = MPI_Wtime()
		
		call info(this,t_E-t_S,1)

    end subroutine
!========================================================================================!

!========================================================================================!
    subroutine computeBlockST(this,vofb)
    	type(VOF), intent(in) :: this
		type(vofBlock), intent(inout) :: vofb
		integer :: lbi,ubi,lbj,ubj,lbk,ubk,i,j,k
		real(DP) :: nf,kf
    	
    	!stx
		lbi = lbound(vofb%stx,1)
		ubi = ubound(vofb%stx,1)
		lbj = lbound(vofb%stx,2)
		ubj = ubound(vofb%stx,2)
		lbk = lbound(vofb%stx,3)
		ubk = ubound(vofb%stx,3)
		
		do k=lbk,ubk
			do j=lbj,ubj
				do i=lbi,ubi
				
					nf = (vofb%c(i+1,j,k)-vofb%c(i,j,k))/vofb%dxc(i+1)

					if ((vofb%isMixed(i+1,j,k)) .AND. (vofb%isMixed(i,j,k))) then
						kf = 0.5d0*(vofb%k(i+1,j,k)+vofb%k(i,j,k))
					else if (vofb%isMixed(i+1,j,k)) then
						kf = vofb%k(i+1,j,k)
					else if (vofb%isMixed(i,j,k)) then
						kf = vofb%k(i,j,k)
					else
						kf = 0.d0
					end if
					
					vofb%stx(i,j,k)=this%sigma_*kf*nf	

				end do
			end do
		end do	
		
    	!sty
		lbi = lbound(vofb%sty,1)
		ubi = ubound(vofb%sty,1)
		lbj = lbound(vofb%sty,2)
		ubj = ubound(vofb%sty,2)
		lbk = lbound(vofb%sty,3)
		ubk = ubound(vofb%sty,3)
		
		do k=lbk,ubk
			do j=lbj,ubj
				do i=lbi,ubi
				
					nf = (vofb%c(i,j+1,k)-vofb%c(i,j,k))/vofb%dyc(j+1)
					
					if ((vofb%isMixed(i,j+1,k)) .AND. (vofb%isMixed(i,j,k))) then
						kf = 0.5d0*(vofb%k(i,j+1,k)+vofb%k(i,j,k))
					else if (vofb%isMixed(i,j+1,k)) then
						kf = vofb%k(i,j+1,k)
					else if (vofb%isMixed(i,j,k)) then
						kf = vofb%k(i,j,k)
					else
						kf = 0.d0
					end if
					
					vofb%sty(i,j,k)=this%sigma_*kf*nf	

				end do
			end do
		end do	
		
    	!stz
		lbi = lbound(vofb%stz,1)
		ubi = ubound(vofb%stz,1)
		lbj = lbound(vofb%stz,2)
		ubj = ubound(vofb%stz,2)
		lbk = lbound(vofb%stz,3)
		ubk = ubound(vofb%stz,3)
		
		do k=lbk,ubk
			do j=lbj,ubj
				do i=lbi,ubi
				
					nf = (vofb%c(i,j,k+1)-vofb%c(i,j,k))/vofb%dzc(k+1)
					
					if ((vofb%isMixed(i,j,k+1)) .AND. (vofb%isMixed(i,j,k))) then
						kf = 0.5d0*(vofb%k(i,j,k+1)+vofb%k(i,j,k))
					else if (vofb%isMixed(i,j,k+1)) then
						kf = vofb%k(i,j,k+1)	
					else if (vofb%isMixed(i,j,k)) then
						kf = vofb%k(i,j,k)	
					else
						kf = 0.d0
					end if

					vofb%stz(i,j,k)=this%sigma_*kf*nf	

				end do
			end do
		end do

    end subroutine
!========================================================================================!


!========================================================================================!
    subroutine computeCurvature(this,k)
    	type(VOF), intent(in) :: this
    	type(scalarField), intent(inout) :: k
    	integer :: b

		!$OMP PARALLEL DO DEFAULT(none) &
		!$OMP SHARED(vofBlocks,this,s_nblk) &
		!$OMP PRIVATE(b)    	
    	do b=1,s_nblk
    		call computeBlockCurvature(this,vofBlocks(b))
    		!call spreadCurvature(vofBlocks(b))
    	end do
    	!$OMP END PARALLEL DO
    
    	call excVFK(this%mesh_,k,2)

    end subroutine
!========================================================================================!

!========================================================================================!
    subroutine spreadCurvature(vofb)
    	type(vofBlock), intent(inout) :: vofb
    	integer :: lbi,ubi,lbj,ubj,lbk,ubk
    	integer :: i,j,k,is,js,ks,count
    	real(DP) :: kav
		
		lbi = lbound(vofb%k,1)
		ubi = ubound(vofb%k,1)
		lbj = lbound(vofb%k,2)
		ubj = ubound(vofb%k,2)
		lbk = lbound(vofb%k,3)
		ubk = ubound(vofb%k,3)
		
		do k=lbk+1,ubk-1
			do j=lbj+1,ubj-1
				do i=lbi+1,ubi-1
				
					if (vofb%k(i,j,k)==0.d0) then
					
						count = 0
						kav=0.d0
						
						!search the stencil			
						do ks=-1,1
							do js=-1,1
								do is=-1,1
									if (vofb%isMixed(i+is,j+js,k+ks)) then
										kav = kav+vofb%k(i+is,j+js,k+ks)
										count = count+1
									end if
								end do
							end do
						end do	
						
						if (count>0) then
							vofb%k(i,j,k) = kav/count
						end if
						
					end if
					
				end do
			end do
		end do				
    
    end subroutine
!========================================================================================!

!========================================================================================!
    subroutine computeBlockCurvature(this,vofb)
! ************************************************************************************** !
! This subroutine is an implementation based on the GHF method developed in:
! Popinet, Stéphane. "An accurate adaptive solver for surface-tension-driven interfacial 
! flows." Journal of Computational Physics 228.16 (2009): 5838-5866.
! ************************************************************************************** !
        type(VOF), intent(in) :: this
        type(vofBlock), intent(inout) :: vofb
        integer :: i, j, k
        real(DP) :: mx, my, mz
        real(DP) :: mmax
        integer :: mloc
        real(DP), dimension(3) :: mv
        real(DP), dimension(9) :: h
        integer, dimension(9) :: lb
        logical, dimension(9) :: isValid
        logical, dimension(27) :: isBlockValid
        integer, dimension(3,9) :: d
        real(DP), dimension(3,9) :: pos_hf
        real(DP), dimension(3,27) :: pos_block
        logical :: hfk,failed
        integer :: cn
    	integer :: lbi,ubi,lbj,ubj,lbk,ubk
    	
		lbi = lbound(vofb%k,1)
		ubi = ubound(vofb%k,1)
		lbj = lbound(vofb%k,2)
		ubj = ubound(vofb%k,2)
		lbk = lbound(vofb%k,3)
		ubk = ubound(vofb%k,3)

		do k=lbk,ubk
			do j=lbj,ubj
				do i=lbi,ubi
				
					!reset curvature
					vofb%k(i,j,k) = 0.d0
					hfk = .TRUE.
				
					if (vofb%isMixed(i,j,k)) then
					
					mx = vofb%nx(i,j,k)
					my = vofb%ny(i,j,k)
					mz = vofb%nz(i,j,k)
					
					mv = (/ mx, my, mz /)
					mloc = maxloc( abs(mv),1 )
					mmax = mv(mloc)
					
					!chiama procedura colonna
					call setStencilPar(mloc,d)
					do cn=1,9
						call hfColumn(this,vofb,mmax,mloc,i+d(1,cn),j+d(2,cn),k+d(3,cn),cn,h,lb,pos_hf,isValid)
						if (.NOT. (isValid(cn))) then
							hfk = .FALSE.
						end if						
					end do
				
					!compute curvature	
					if (hfk) then
						call hfCurvature(vofb,mloc,mmax,i,j,k,h,lb)		
					else
						call searchBlockHF(this,vofb,i,j,k,mv,mloc,isValid,pos_hf,isBlockValid,pos_block)
						call parabFittedCurvature(vofb,i,j,k,pos_block,isBlockValid,.TRUE.,failed)	
						if (failed) then
							call interCentroids(vofb,i,j,k,pos_block,isBlockValid)
							call parabFittedCurvature(vofb,i,j,k,pos_block,isBlockValid,.FALSE.,failed)	
						end if		
					!	call fdCurvature(vofb,i,j,k)			
					end if		

	
					end if
					
				end do
			end do
		end do      
		
		
    end subroutine
!========================================================================================!

!========================================================================================!
    subroutine hfColumn(this,vofb,m,dir,i,j,k,cn,h,lb,posi,isValid)
    	type(VOF), intent(in) :: this
        type(vofBlock), intent(in) :: vofb
        real(DP),intent(in) :: m
        integer, intent(in) :: dir, i, j ,k, cn
        real(DP), dimension(9), intent(inout) :: h
        integer, dimension(9), intent(inout) :: lb
        logical, dimension(9), intent(inout) :: isValid
        real(DP), dimension(3,9), intent(inout) :: posi
        integer :: swi,swj,swk,it,jt,kt,ib,jb,kb
        integer :: ii, jj, kk, l
        real(DP), pointer, dimension(:) :: df => NULL()
        real(DP), pointer, dimension(:) :: posf => NULL()
        !$omp threadprivate(df,posf)


		!chiama set par      
		call setColumnPar(m,dir,swi,swj,swk,it,jt,kt,ib,jb,kb)
		
		select case(dir)
			case(1)
				df => vofb%dxf
				posf => vofb%xf
			case(2)
				df => vofb%dyf
				posf => vofb%yf
			case(3)
				df => vofb%dzf
				posf => vofb%zf
			case default
		end select
		
		h(cn) = vofb%c(i,j,k)*df(swi*i+swj*j+swk*k)
		
		!search the top column
		isValid(cn) = .FALSE.
		do l=1,this%hd_
			ii = i+l*it
			jj = j+l*jt
			kk = k+l*kt
			h(cn) = h(cn) + vofb%c(ii,jj,kk)*df(swi*ii+swj*jj+swk*kk)
			
			if ((.NOT. (vofb%isFull(ii,jj,kk))) .AND. (.NOT. (vofb%isMixed(ii,jj,kk)))) then
				isValid(cn) = .TRUE.
				exit
			end if
			
		end do
		
		if (.NOT. (isValid(cn))) then
			posi(:,cn)=0.d0
			return;
		end if
		
		!search the bottom column
		isValid(cn) = .FALSE.
		do l=1,this%hd_
			ii = i+l*ib
			jj = j+l*jb
			kk = k+l*kb
			h(cn) = h(cn) + vofb%c(ii,jj,kk)*df(swi*ii+swj*jj+swk*kk)
			
			if (vofb%isFull(ii,jj,kk)) then
				isValid(cn) = .TRUE.
				lb(cn) = l
				exit
			end if
			
		end do
		
		if (isValid(cn)) then
			!store interface position
			posi(1,cn)=vofb%xc(i)
			posi(2,cn)=vofb%yc(j)
			posi(3,cn)=vofb%zc(k)
			if (m >= 0.d0) then
				posi(dir,cn)=posf(swi*i+swj*j+swk*k-lb(cn)-1)+h(cn)
			else
				posi(dir,cn)=posf(swi*i+swj*j+swk*k+lb(cn))-h(cn)
			end if
		else
			posi(:,cn)=0.d0
		end if
		
    end subroutine
!========================================================================================!

!========================================================================================!
    subroutine hfCurvature(vofb,dir,m,i,j,k,h,lb)
        type(vofBlock), intent(inout) :: vofb
        integer, intent(in) :: dir, i, j, k
        real(DP), intent(in) :: m
        real(DP), dimension(9), intent(inout) :: h
        integer, dimension(9), intent(inout) :: lb
		real(DP) :: dhxx, dhyy, dhxy, dhx, dhy, denom, num
		integer :: lmax
		integer :: nc,n1,n2,dir1,dir2
		integer :: q
        real(DP), pointer, dimension(:) :: df => NULL()
        real(DP), pointer, dimension(:) :: df1 => NULL()
        real(DP), pointer, dimension(:) :: df2 => NULL()
        real(DP), pointer, dimension(:) :: dc1 => NULL()
        real(DP), pointer, dimension(:) :: dc2 => NULL()
        !$omp threadprivate(df,df1,df2,dc1,dc2)
		
		select case(dir)
			case(1)
				df => vofb%dxf
			case(2)
				df => vofb%dyf
			case(3)
				df => vofb%dzf
			case default
		end select
		
		lmax = maxval( lb,1 )
		
		call setDirections(dir,i,j,k,nc,n1,n2,dir1,dir2)
		
		select case(dir1)
			case(1)
				df1 => vofb%dxf
				dc1 => vofb%dxc
			case(2)
				df1 => vofb%dyf
				dc1 => vofb%dyc
			case(3)
				df1 => vofb%dzf
				dc1 => vofb%dzc
			case default
		end select
		
		select case(dir2)
			case(1)
				df2 => vofb%dxf
				dc2 => vofb%dxc
			case(2)
				df2 => vofb%dyf
				dc2 => vofb%dyc
			case(3)
				df2 => vofb%dzf
				dc2 => vofb%dzc
			case default
		end select
		
		!make reference equal
		if (m < 0.d0) then
			do q=1,9
				h(q) = h(q) + sum(df(nc+lb(q)+1:nc+lmax))
			end do
		else
			do q=1,9
				h(q) = h(q) + sum(df(nc-lmax:nc-lb(q)-1))	
			end do	
		end if
		
		dhx = (h(6) - h(4))/(dc1(n1+1)+dc1(n1))
		dhy = (h(8) - h(2))/(dc2(n2+1)+dc2(n2))
		dhxx = ((h(6)-h(5))/(dc1(n1+1))-(h(5)-h(4))/(dc1(n1))) / & 
			   df1(n1)
		dhyy = ((h(8)-h(5))/(dc2(n2+1))-(h(5)-h(2))/(dc2(n2))) / &
				df2(n2)
		dhxy = ((h(9)-h(7)-h(3)+h(1))/(dc1(n1+1)+dc1(n1))) / &
				(dc2(n2+1)+dc2(n2))
		
		!3D curvature
		num = dhxx*(1.d0+dhy*dhy)+dhyy*(1.d0+dhx*dhx)-2.d0*dhxy*dhx*dhy
		denom = (1.d0+dhx*dhx+dhy*dhy)*sqrt(1.d0+dhx*dhx+dhy*dhy)
		vofb%k(i,j,k) = -num/denom

		
    end subroutine
!========================================================================================!

!========================================================================================!
    subroutine searchBlockHF(this,vofb,i,j,k,m,mloc,isValid,pos_hf,isBlockValid,pos_block)
    	type(VOF), intent(in) :: this
    	type(vofBlock), intent(in) :: vofb
    	integer, intent(in) :: i,j,k
    	real(DP), dimension(3), intent(in) :: m
    	integer, intent(in) :: mloc
    	real(DP), dimension(3,9), intent(inout) :: pos_hf
    	logical, dimension(9), intent(inout) :: isValid
    	logical, dimension(27), intent(out) :: isBlockValid
    	real(DP), dimension(3,27), intent(out) :: pos_block
    	integer, dimension(3,9) :: d
    	integer :: cn,mloc_2,mloc_3
        real(DP), dimension(9) :: h
        integer, dimension(9) :: lb
        real(DP) :: m_2,m_3
		

		!store search dir 1
		isBlockValid(1:9)=isValid
		pos_block(:,1:9)=pos_hf
		
		select case(mloc)
			case(1)
				mloc_2=2
				m_2=m(2)
				mloc_3=3
				m_3=m(3)
			case(2)
				mloc_2=1
				m_2=m(1)
				mloc_3=3
				m_3=m(3)
			case(3)
				mloc_2=1
				m_2=m(1)
				mloc_3=2
				m_3=m(2)
			case default
		end select		
		
		!search dir 2
		call setStencilPar(mloc_2,d)
		do cn=1,9
			call hfColumn(this,vofb,m_2,mloc_2,i+d(1,cn),j+d(2,cn),k+d(3,cn),cn,h,lb,pos_hf,isValid)
		end do	
		isBlockValid(10:18)=isValid
		pos_block(:,10:18)=pos_hf			
		
		!search dir 3
		call setStencilPar(mloc_3,d)
		do cn=1,9
			call hfColumn(this,vofb,m_3,mloc_3,i+d(1,cn),j+d(2,cn),k+d(3,cn),cn,h,lb,pos_hf,isValid)
		end do	
		isBlockValid(19:27)=isValid
		pos_block(:,19:27)=pos_hf		
		
		
    end subroutine
!========================================================================================!

!========================================================================================!
    subroutine independentPoints(ds,vin,isValid,fromHF,vout,np,lab_c0,enoughPoints)
    	real(DP), dimension(3), intent(in) :: ds
    	real(DP), dimension(3,27), intent(inout) :: vin
    	real(DP), dimension(3,27), intent(out) :: vout
    	logical, dimension(27), intent(in) :: isValid
    	logical, intent(in) :: fromHF
    	integer, intent(out) :: np, lab_c0
    	logical, intent(out) :: enoughPoints
    	integer :: ii,jj,kk
    	real(DP) :: d,dx,dy,dz,dmin,delta
    	
    	!check for valid central point first
    	if (fromHF) then
    		if (.NOT.isValid(5)) then
    			if (.NOT.isValid(14)) then
    				if (.NOT.isValid(23)) then
    					enoughPoints=.FALSE.
					return
    				else
    					lab_c0 = 23
    				end if
    			else
    				lab_c0 = 14
    			end if
    		else
    			lab_c0 = 5
    		end if
    	else
    		if (.NOT.isValid(14)) then
    			enoughPoints=.FALSE.
				return
    		else
    			lab_c0 = 14
    		end if
    	end if

    	np = 0
    	vout = 0.d0
    	dmin = 1.d-2*min(ds(1),ds(2),ds(3))
    	!dmin = 0.5d0*min(ds(1),ds(2),ds(3))
    	delta = dmin*dmin 
    	!delta = 0.d0
 	
    	!search for independent points
    	do ii=1,size(vin,2)
    		if (isValid(ii)) then
				vout(:,1)=vin(:,ii)
				np=1
				exit
    		end if
    	end do
    	do jj=ii+1,size(vin,2)
    		if (isValid(jj)) then
    			do kk=1,np
    	    		dx = vout(1,kk)-vin(1,jj)
    				dy = vout(2,kk)-vin(2,jj)
    				dz = vout(3,kk)-vin(3,jj)
    				d = dx*dx+dy*dy+dz*dz
    				if (d<=delta) then
    					goto 10
    				end if
    			end do
    			vout(:,np+1)=vin(:,jj)
    			np=np+1
    		end if
10    	end do

		if (np<6) then
			enoughPoints=.FALSE.
			return
		end if
		
		enoughPoints=.TRUE.
		
    end subroutine
!========================================================================================!

!========================================================================================!
    subroutine localIntPos(vg,vl,np,m,p0)
    	real(DP), dimension(3,27), intent(in) :: vg
    	real(DP), dimension(3,27), intent(out) :: vl
    	real(DP), dimension(3), intent(in) :: m,p0
    	integer, intent(in) :: np
    	real(DP), dimension(3) :: i1,i2,i3
    	real(DP) :: s1,s2
    	integer :: j,mloc
    	
    	vl = 0.d0
    	
    	s2 = sqrt(m(1)*m(1)+m(2)*m(2)+m(3)*m(3))
    	i2=m/s2
    	
    	mloc = maxloc( abs(m),1 )
    	
		select case(mloc)
    		case(1)
        		mloc=maxloc(abs((/m(2),m(3)/)),1)
        		select case(mloc)
            		case(1)
               			i1(1)=-i2(2)
               			i1(2)=i2(1)
               			i1(3)=0.d0
               			s1=sqrt(i1(1)*i1(1)+i1(2)*i1(2))
            		case(2)
               			i1(1)=-i2(3)
               			i1(2)=0.d0
               			i1(3)=i2(1)
               			s1=sqrt(i1(1)*i1(1)+i1(3)*i1(3))
            		case default
        		end select
    		case(2)
        		mloc=maxloc(abs((/m(1),m(3)/)),1)
        		select case(mloc)
            		case(1)
               			i1(1)=-i2(2)
               			i1(2)=i2(1)
               			i1(3)=0.d0
               			s1=sqrt(i1(1)*i1(1)+i1(2)*i1(2))
            		case(2)
               			i1(1)=0.d0
               			i1(2)=-i2(3)
               			i1(3)=i2(2)
               			s1=sqrt(i1(2)*i1(2)+i1(3)*i1(3))
            		case default
        		end select
     		case(3)
        		mloc=maxloc(abs((/m(1),m(2)/)),1)
        		select case(mloc)
            		case(1)
               			i1(1)=-i2(3)
               			i1(2)=0.d0
               			i1(3)=i2(1)
               			s1=sqrt(i1(1)*i1(1)+i1(3)*i1(3))
            		case(2)
               			i1(1)=0.d0
               			i1(2)=-i2(3)
               			i1(3)=i2(2)
               			s1=sqrt(i1(2)*i1(2)+i1(3)*i1(3))
            		case default
        		end select
    		case default
		end select
		 	
    	i1=i1/s1

		i3(1) = i1(2)*i2(3)-i1(3)*i2(2)
		i3(2) = i1(3)*i2(1)-i1(1)*i2(3)
		i3(3) = i1(1)*i2(2)-i1(2)*i2(1)	

		do j=1,np
    		vl(1,j) = i1(1)*(vg(1,j)-p0(1)) + i1(2)*(vg(2,j)-p0(2)) + i1(3)*(vg(3,j)-p0(3))
			vl(2,j) = i2(1)*(vg(1,j)-p0(1)) + i2(2)*(vg(2,j)-p0(2)) + i2(3)*(vg(3,j)-p0(3))
			vl(3,j) = i3(1)*(vg(1,j)-p0(1)) + i3(2)*(vg(2,j)-p0(2)) + i3(3)*(vg(3,j)-p0(3))	
		end do

		
    end subroutine
!========================================================================================!

!========================================================================================!
    subroutine buildLSqSystem(v,np,A,b)
		real(DP), dimension(3,27), intent(in) :: v
		integer, intent(in) :: np
		real(DP), dimension(6,6), intent(out) :: A
		real(DP), dimension(6), intent(out) :: b
		real(DP) :: q
		integer :: i
		
		A=0.d0
		b=0.d0

		!col 1
		do i=1,np
			q=v(1,i)*v(1,i)
			A(1,1)=A(1,1)+q*v(1,i)*v(1,i)
			A(2,1)=A(2,1)+q*v(3,i)*v(3,i)
			A(3,1)=A(3,1)+q*v(1,i)*v(3,i)
			A(4,1)=A(4,1)+q*v(1,i)
			A(5,1)=A(5,1)+q*v(3,i)
			A(6,1)=A(6,1)+q
		end do
		
		!col 2
		do i=1,np
			q=v(3,i)*v(3,i)
			A(1,2)=A(1,2)+q*v(1,i)*v(1,i)
			A(2,2)=A(2,2)+q*v(3,i)*v(3,i)
			A(3,2)=A(3,2)+q*v(1,i)*v(3,i)
			A(4,2)=A(4,2)+q*v(1,i)
			A(5,2)=A(5,2)+q*v(3,i)
			A(6,2)=A(6,2)+q
		end do
		
		!col 3
		do i=1,np
			q=v(1,i)*v(3,i)
			A(1,3)=A(1,3)+q*v(1,i)*v(1,i)
			A(2,3)=A(2,3)+q*v(3,i)*v(3,i)
			A(3,3)=A(3,3)+q*v(1,i)*v(3,i)
			A(4,3)=A(4,3)+q*v(1,i)
			A(5,3)=A(5,3)+q*v(3,i)
			A(6,3)=A(6,3)+q
		end do
		
		!col 4
		do i=1,np
			q=v(1,i)
			A(1,4)=A(1,4)+q*v(1,i)*v(1,i)
			A(2,4)=A(2,4)+q*v(3,i)*v(3,i)
			A(3,4)=A(3,4)+q*v(1,i)*v(3,i)
			A(4,4)=A(4,4)+q*v(1,i)
			A(5,4)=A(5,4)+q*v(3,i)
			A(6,4)=A(6,4)+q
		end do
		
		!col 5
		do i=1,np
			q=v(3,i)
			A(1,5)=A(1,5)+q*v(1,i)*v(1,i)
			A(2,5)=A(2,5)+q*v(3,i)*v(3,i)
			A(3,5)=A(3,5)+q*v(1,i)*v(3,i)
			A(4,5)=A(4,5)+q*v(1,i)
			A(5,5)=A(5,5)+q*v(3,i)
			A(6,5)=A(6,5)+q
		end do
		
		!col 6
		do i=1,np
			q=1.d0
			A(1,6)=A(1,6)+q*v(1,i)*v(1,i)
			A(2,6)=A(2,6)+q*v(3,i)*v(3,i)
			A(3,6)=A(3,6)+q*v(1,i)*v(3,i)
			A(4,6)=A(4,6)+q*v(1,i)
			A(5,6)=A(5,6)+q*v(3,i)
			A(6,6)=A(6,6)+q
		end do
		
		!rhs
		do i=1,np
			b(1)=b(1)+v(2,i)*(v(1,i)*v(1,i))
			b(2)=b(2)+v(2,i)*(v(3,i)*v(3,i))
			b(3)=b(3)+v(2,i)*(v(1,i)*v(3,i))
			b(4)=b(4)+v(2,i)*v(1,i)
			b(5)=b(5)+v(2,i)*v(3,i)
			b(6)=b(6)+v(2,i)
		end do	

		
    end subroutine
!========================================================================================!

!========================================================================================!
    subroutine solveLSqSystem(a,b,singular)
		real(DP), dimension(6,6), intent(inout) :: a
		real(DP), dimension(6), intent(inout) :: b
		logical, intent(out) :: singular
		integer :: n,nrhs,lda,ldb,info
		integer, dimension(6) :: ipiv

		n=6
		nrhs=1
		lda=n
		ldb=n 
		
		call DGESV(n,nrhs,a,lda,ipiv,b,ldb,info)
		
		if (.NOT.(info==0)) then
			singular = .TRUE.
		else
			singular = .FALSE.
		end if
		
		
    end subroutine
!========================================================================================!

!========================================================================================!
    subroutine interCentroids(vofb,i,j,k,cts,isValid)
    	type(vofBlock), intent(in) :: vofb
    	integer, intent(in) :: i,j,k
    	integer :: is,js,ks,ii,jj,kk
    	real(DP), dimension(3,27), intent(out) :: cts
    	logical, dimension(27), intent(out) :: isValid
    	real(DP), dimension(3) :: ds,c0,m,ct
    	real(DP) :: q
    	integer :: count
    	logical :: failed
    	
    	cts = 0.d0
    	isValid = .FALSE.
    	
    	count = 0
    	do ks=-1,1
    		do js=-1,1
    			do is=-1,1
    				ii = i+is
    				jj= j+js
    				kk= k+ks
    				count = count+1
    				if (vofb%isMixed(ii,jj,kk)) then
    					ds(1)=vofb%dxf(ii)
    					ds(2)=vofb%dyf(jj)
    					ds(3)=vofb%dzf(kk)
    					c0(1)=vofb%xf(ii-1)
    					c0(2)=vofb%yf(jj-1)
    					c0(3)=vofb%zf(kk-1)
						m(1)=vofb%nx(ii,jj,kk)
						m(2)=vofb%ny(ii,jj,kk)
						m(3)=vofb%nz(ii,jj,kk)   
						q=vofb%q(ii,jj,kk) 					
    					call centroid(ds,c0,m,q,ct,failed)
    					if (.NOT.failed) then
    						cts(:,count)=ct
    						isValid(count)=.TRUE.
    					end if
    				end if
    			end do
    		end do
    	end do
		
		
    end subroutine
!========================================================================================!

!========================================================================================!
    subroutine parabFittedCurvature(vofb,i,j,k,posi,isValid,fromHF,failed)
    	type(vofBlock), intent(inout) :: vofb
    	integer, intent(in) :: i,j,k
    	real(DP), dimension(3,27), intent(inout) :: posi
    	logical, dimension(27), intent(inout) :: isValid
    	logical, intent(in) :: fromHF
    	logical, intent(out) :: failed
    	real(DP), dimension(3,27) :: posi_ind, ploc
    	real(DP), dimension(3) :: m,ds,c0
		real(DP), dimension(6,6) :: A
		real(DP), dimension(6) :: b
    	real(DP) :: num,denom,r
    	integer :: np, lab_c0
    	logical :: enoughPoints,singular
    	
    	ds(1)=vofb%dxf(i)
    	ds(2)=vofb%dyf(j)
    	ds(3)=vofb%dzf(k)
		m(1)=vofb%nx(i,j,k)
		m(2)=vofb%ny(i,j,k)
		m(3)=vofb%nz(i,j,k) 

		call independentPoints(ds,posi,isValid,fromHF,posi_ind,np,lab_c0,enoughPoints)
		if (enoughPoints) then
			c0 = posi(:,lab_c0)
			call localIntPos(posi_ind,ploc,np,m,c0)
			call buildLSqSystem(ploc,np,A,b)
			call solveLSqSystem(A,b,singular)
		
			if (singular) then
				failed = .TRUE.
				return
			else
				!calcola curv
				num=b(1)*(1.d0+b(5)*b(5))+b(2)*(1.d0+b(4)*b(4))-b(3)*b(4)*b(5)
				r=1.d0+b(4)*b(4)+b(5)*b(5)
				denom=r*sqrt(r)
				vofb%k(i,j,k)=-2.d0*num/denom
			end if
		else
			failed = .TRUE.
			return
		end if
		
		failed = .FALSE.
		
    end subroutine
!========================================================================================!

!========================================================================================!
    subroutine fdCurvature(vofb,i,j,k)
        type(vofBlock), intent(inout) :: vofb
        integer, intent(in) :: i,j,k
        real(DP), dimension(3) :: nihp,nihm,njhp,njhm,nkhp,nkhm
        real(DP) :: uihp,uihm,ujhp,ujhm,ukhp,ukhm
        real(DP) :: dx, dy, dz

		!normal vector
		!face i+1/2
		nihp(1) = (vofb%cs(i+1,j,k)-vofb%cs(i,j,k))/vofb%dxc(i+1)
		nihp(2) = 0.5d0*(										&
						vofb%cs(i,j+1,k)-vofb%cs(i,j-1,k)+		&
						vofb%cs(i+1,j+1,k)-vofb%cs(i+1,j-1,k)	&
						)/(vofb%dyc(j+1)+vofb%dyc(j))
		nihp(3) = 0.5d0*(										&
						vofb%cs(i,j,k+1)-vofb%cs(i,j,k-1)+		&
						vofb%cs(i+1,j,k+1)-vofb%cs(i+1,j,k-1)	&
						)/(vofb%dzc(k+1)+vofb%dzc(k))
		
		!face i-1/2
		nihm(1) = (vofb%cs(i,j,k)-vofb%cs(i-1,j,k))/vofb%dxc(i)
		nihm(2) = 0.5d0*(										&
						vofb%cs(i,j+1,k)-vofb%cs(i,j-1,k)+		&
						vofb%cs(i-1,j+1,k)-vofb%cs(i-1,j-1,k)	&
						)/(vofb%dyc(j+1)+vofb%dyc(j))
		nihm(3) = 0.5d0*(										&
						vofb%cs(i,j,k+1)-vofb%cs(i,j,k-1)+		&
						vofb%cs(i-1,j,k+1)-vofb%cs(i-1,j,k-1)	&
						)/(vofb%dzc(k+1)+vofb%dzc(k))
						
		!face j+1/2
		njhp(1) = 0.5d0*(										&
						vofb%cs(i+1,j,k)-vofb%cs(i-1,j,k)+		&
						vofb%cs(i+1,j+1,k)-vofb%cs(i-1,j+1,k)	&
						)/(vofb%dxc(i+1)+vofb%dxc(i))
		
		njhp(2) = (vofb%cs(i,j+1,k)-vofb%cs(i,j,k))/vofb%dyc(j+1)
		njhp(3) = 0.5d0*(										&
						vofb%cs(i,j,k+1)-vofb%cs(i,j,k-1)+		&
						vofb%cs(i,j+1,k+1)-vofb%cs(i,j+1,k-1)	&
						)/(vofb%dzc(k+1)+vofb%dzc(k))

		!face j-1/2
		njhm(1) = 0.5d0*(										&
						vofb%cs(i+1,j,k)-vofb%cs(i-1,j,k)+		&
						vofb%cs(i+1,j-1,k)-vofb%cs(i-1,j-1,k)	&
						)/(vofb%dxc(i+1)+vofb%dxc(i))
		
		njhm(2) = (vofb%cs(i,j,k)-vofb%cs(i,j-1,k))/vofb%dyc(j)
		njhm(3) = 0.5d0*(										&
						vofb%cs(i,j,k+1)-vofb%cs(i,j,k-1)+		&
						vofb%cs(i,j-1,k+1)-vofb%cs(i,j-1,k-1)	&
						)/(vofb%dzc(k+1)+vofb%dzc(k))
						
		!face k+1/2
		nkhp(1) = 0.5d0*(										&
						vofb%cs(i+1,j,k)-vofb%cs(i-1,j,k)+		&
						vofb%cs(i+1,j,k+1)-vofb%cs(i-1,j,k+1)	&
						)/(vofb%dxc(i+1)+vofb%dxc(i))
		
		nkhp(2) = 0.5d0*(										&
						vofb%cs(i,j+1,k)-vofb%cs(i,j-1,k)+		&
						vofb%cs(i,j+1,k+1)-vofb%cs(i,j-1,k+1)	&
						)/(vofb%dyc(j+1)+vofb%dyc(j)) 
		nkhp(3) = (vofb%cs(i,j,k+1)-vofb%cs(i,j,k))/vofb%dzc(k+1)
		
		!face k-1/2
		nkhm(1) = 0.5d0*(										&
						vofb%cs(i+1,j,k)-vofb%cs(i-1,j,k)+		&
						vofb%cs(i+1,j,k-1)-vofb%cs(i-1,j,k-1)	&
						)/(vofb%dxc(i+1)+vofb%dxc(i))
		
		nkhm(2) = 0.5d0*(										&
						vofb%cs(i,j+1,k)-vofb%cs(i,j-1,k)+		&
						vofb%cs(i,j+1,k-1)-vofb%cs(i,j-1,k-1)	&
						)/(vofb%dyc(j+1)+vofb%dyc(j)) 
		nkhm(3) = (vofb%cs(i,j,k)-vofb%cs(i,j,k-1))/vofb%dzc(k)
		
		!unit normal components
		!face i+1/2
		uihp = nihp(1)/sqrt(nihp(1)*nihp(1)+nihp(2)*nihp(2)+nihp(3)*nihp(3))
		!face i-1/2	
		uihm = nihm(1)/sqrt(nihm(1)*nihm(1)+nihm(2)*nihm(2)+nihm(3)*nihm(3))
		!face j+1/2
		ujhp = njhp(2)/sqrt(njhp(1)*njhp(1)+njhp(2)*njhp(2)+njhp(3)*njhp(3))
		!face j-1/2	
		ujhm = njhm(2)/sqrt(njhm(1)*njhm(1)+njhm(2)*njhm(2)+njhm(3)*njhm(3))
		!face k+1/2
		ukhp = nkhp(3)/sqrt(nkhp(1)*nkhp(1)+nkhp(2)*nkhp(2)+nkhp(3)*nkhp(3))
		!face k-1/2	
		ukhm = nkhm(3)/sqrt(nkhm(1)*nkhm(1)+nkhm(2)*nkhm(2)+nkhm(3)*nkhm(3))
		
		dx = vofb%dxf(i)
		dy = vofb%dyf(j)
		dz = vofb%dzf(k)
		
		vofb%k(i,j,k) = -( (uihp-uihm)/dx + (ujhp-ujhm)/dy + (ukhp-ukhm)/dz )
		
		
		
    end subroutine
!========================================================================================!

!******************************* Interface Geom. props **********************************!

!========================================================================================!
    subroutine qCell(vofb,i,j,k,m0,v,delta,q)
! ************************************************************************************** !
! This subroutine is an implementation of the reconstruction routines developed in:  
! Scardovelli, Ruben, and Stephane Zaleski. "Analytical relations connecting linear 
! interfaces and volume fractions in rectangular grids." 
! Journal of Computational Physics 164.1 (2000): 228-237.
! ************************************************************************************** !
    	type(vofBlock), intent(inout) :: vofb
    	integer, intent(in) :: i,j,k
    	real(DP), dimension(3), intent(in) :: m0, delta
    	real(DP), intent(inout) :: v
    	real(DP), intent(out) :: q
    	real(DP), dimension(3) :: m, ms
    	real(DP) :: m_par, m1, m2, m3
    	real(DP) :: m12,V1,V2,V3,vh
    	real(DP) :: m1m2m3,m1p2,m2p2,m3p2,m1p3,m2p3,m3p3
    	real(DP) :: p0,q0,rp0,theta,sn,cs
    	integer :: l_min,l_max
    
		!scale size
		ms(1)=m0(1)*delta(1)
		ms(2)=m0(2)*delta(2)
		ms(3)=m0(3)*delta(3)
	
		!mirror around h/2
		m(1)=abs(ms(1))
		m(2)=abs(ms(2))
		m(3)=abs(ms(3))
	
		!scaled plane equation
		m_par=m(1)+m(2)+m(3)
	
		!check normal interface cell
		if (m_par == 0.d0) then
			if (v - 0.5d0 > 0.d0) then
!				write(*,'(A,3'//s_intFormat(2:3)//',A,'//s_outputFormat(2:9)//',A)') &
!							"Warning: failed normal cell ", i,j,k, &
!			 				" --> c = ", v, &
!			 				" reset to 1"
				vofb%c(i,j,k) = 1.d0
				vofb%isFull(i,j,k) = .TRUE.
				q = 0.d0
				return 
			end if
		
			if (v - 0.5d0 < 0.d0) then
!				write(*,'(A,3'//s_intFormat(2:3)//',A,'//s_outputFormat(2:9)//',A)') &
!							"Warning: failed normal cell ", i,j,k, &
!			 				" --> c = ", v, &
!			 				" reset to 0"
				vofb%c(i,j,k) = 0.d0
				vofb%isFull(i,j,k) = .FALSE.
				vofb%isMixed(i,j,k) = .FALSE.
				q = 0.d0
				return 
			end if
		end if


		m = m/m_par
	
		!permutation
  		l_min=minloc(m,1)
  		l_max=maxloc(m,1)
  		m1=m(l_min)
  		m3=m(l_max)
  		m(l_min)=-1.d0
  		m(l_max)=-1.d0
  		m2=maxval(m)

		!calculates bounds
		V1=(m1*m1)/max(6.d0*m2*m3,tiny(0.d0))
		V2=(m2-m1)/(2*m3)+V1
		m12=m1+m2
		m1m2m3=max(m1*m2*m3,tiny(0.d0))
		if (m3<m12) then
   		 	V3=(m3*m3*(3.d0*m12-m3)+m1*m1*(m1-3.d0*m3)+m2*m2*(m2-3.d0*m3))/ &
   		 	    (6.d0*m1m2m3)
		else
   		 	V3=m12/(2.d0*m3)
		end if

		!check bounds
		vh=min(v,1.d0-v)		
		if (vh<=V1) then
    		q=(6.d0*m1m2m3*vh)**(1.d0/3.d0)
		else if (vh<=V2) then
   			q=0.5d0*(m1+sqrt(m1*m1+8.d0*m2*m3*(vh-V1)))
		else if (vh<=V3) then
    		p0= -2.d0*m1*m2
    		q0=3.d0*m1*m2*(0.5d0*m12-m3*vh)
    		rp0 = sqrt(-p0)
    		theta=acos(q0/(rp0*(-p0)))/3.d0
    		sn = sin(theta)
    		cs = sqrt(1.d0-sn*sn)
   			q=rp0*(sqrt(3.d0)*sn-cs)+m12
		else if (m3<m12) then

			m1p2=m1*m1
			m2p2=m2*m2
			m3p2=m3*m3
			m1p3=m1*m1p2
			m2p3=m2*m2p2
			m3p3=m3*m3p2

    		p0=0.5d0*(m1p2+m2p2+m3p2-0.5d0)
    		q0=(-3.d0*(m1p2+m2p2+m3p2)+2.d0*(m1p3+m2p3+m3p3)- &
    			12.d0*(m1m2m3*vh)+1.d0)/8.d0
    		rp0 = sqrt(-p0)
    		theta=acos(q0/(rp0*(-p0)))/3.d0
    		sn = sin(theta)
    		cs = sqrt(1.d0-sn*sn)
    		q=rp0*(sqrt(3.d0)*sn-cs)+0.5d0
		else
		    q=m3*vh+0.5d0*m12
    	end if

		!symmetry
		if (v>0.5d0) then
    		q=1-q
    	end if

		q=q*m_par
		
	
		!get back to the coord. system

		if (m0(1) < 0.d0) then
  	  		q = q + ms(1)
  	  	end if
		if (m0(2) < 0.d0) then
 	   		q = q + ms(2)
 	   	end if
		if (m0(3) < 0.d0) then
  	  		q = q + ms(3)
  	  	end if

    end subroutine
!========================================================================================!

!========================================================================================!
    subroutine vCell(vofb,i,j,k,m0,q0,delta,v)
! ************************************************************************************** !
! This subroutine is an implementation of the reconstruction routines developed in:  
! Scardovelli, Ruben, and Stephane Zaleski. "Analytical relations connecting linear 
! interfaces and volume fractions in rectangular grids." 
! Journal of Computational Physics 164.1 (2000): 228-237.
! ************************************************************************************** !
    	type(vofBlock), intent(inout) :: vofb
    	integer, intent(in) :: i,j,k
    	real(DP), dimension(3), intent(in) :: m0, delta
    	real(DP), intent(in) :: q0
    	real(DP), intent(out) :: v
    	real(DP), dimension(3) :: ms, m
    	real(DP) :: m_par, m1, m2, m3, q, qH, m12, m4th, m1m2m3, V1
    	integer :: l_min,l_max
	
		q = q0
	
		!scale size
		ms(1)=m0(1)*delta(1)
		ms(2)=m0(2)*delta(2)
		ms(3)=m0(3)*delta(3)
	
		!mirror around h/2
		m(1)=abs(ms(1))
		m(2)=abs(ms(2))
		m(3)=abs(ms(3))
	
		!scaled plane equation
		m_par=m(1)+m(2)+m(3)
	
		!check normal interface cell
		if (m_par == 0.d0) then
!			write(*,'(A,3'//s_intFormat(2:3)//',A,'//s_outputFormat(2:9)//',A)') &
!						"Warning: failed normal cell ", i,j,k, &
!			 			" --> c = ", v, &
!			 			" reset to 0"
				vofb%c(i,j,k) = 0.d0
				vofb%isFull(i,j,k) = .FALSE.
				vofb%isMixed(i,j,k) = .FALSE.
				v = 0.d0
				return 
		end if
	
		m = m/m_par
	
		!permutation
  		l_min=minloc(m,1)
  		l_max=maxloc(m,1)
  		m1=m(l_min)
  		m3=m(l_max)
  		m(l_min)=-1.d0
  		m(l_max)=-1.d0
  		m2=maxval(m)
	
		q = q+max(0.d0,-ms(1))+max(0.d0,-ms(2))+max(0.d0,-ms(3))
		q=q/m_par
  	
  		if (q > 1.d0) then
   			v = 1.d0
   			return
   		end if
		if (q < 0.d0) then
  		 	v = 0.d0
  		 	return
  		end if
  	
  	
  		qH = min(q,1.d0-q)
		m12=m1+m2
		m4th=min(m12,m3)
		m1m2m3=max(m1*m2*m3,tiny(0.d0))
		V1=(m1*m1)/max(6.d0*m2*m3,tiny(0.d0))

		if (qH<=m1) then
    		v = (qH*qH*qH)/(6.d0*m1m2m3)
		else if (qH<=m2) then
    		v=0.5d0*qH*(qH - m1)/(m2*m3)+V1
		else if (qH<=m4th) then
    		v=(qH*qH*(3.d0*m12-qH)+m1*m1*(m1-3.d0*qH)+m2*m2*(m2-3.d0*qH))/(6.d0*m1m2m3)
		else if (m3<m12) then
        	v=(qH*qH*(3.d0-2.d0*qH)+m1*m1*(m1-3.d0*qH)+m2*m2*(m2-3.d0*qH)+ & 
        	   m3*m3*(m3-3.d0*qH))/(6.d0*m1m2m3)
		else
    		v = (qH-0.5d0*m4th)/m3
		end if

		!symmetry
		if (q>0.5d0) then
    		v=1.d0-v
    	end if


    end subroutine
!========================================================================================!

!========================================================================================!
    subroutine vArea(m0,q0,delta,area)
    	real(DP), dimension(3), intent(in) :: m0, delta
    	real(DP), intent(in) :: q0
    	real(DP), intent(out) :: area
    	real(DP) :: q,m_par,m1,m2,m3,d1,d2,d3,tmp,qH,m12,m4th,m1d,m2d,m3d
    	real(DP), dimension(3) :: ms,m,i1,i2,i3,p1,p2,p3,p4,p5,p6
    	real(DP), dimension(3) :: p1p,p2p,p3p,p4p,p5p,p6p
    	
    	q=q0
	
		!scale size
		ms(1)=m0(1)*delta(1)
		ms(2)=m0(2)*delta(2)
		ms(3)=m0(3)*delta(3)
	
		!mirror around h/2
		m(1)=abs(ms(1))
		m(2)=abs(ms(2))
		m(3)=abs(ms(3))
	
		!scaled plane equation
		m_par=m(1)+m(2)+m(3)
	
		!check normal interface cell 
		if (m_par==0.d0) then
			call mpiAbort('invalid normal found in vArea')
		end if
		m = m/m_par
	
		!permutation
  		!keep track of cell size
  		if (m(1) >= m(2)) then
    		m3 = m(1)
    		m1 = m(2)
    		d3 = delta(1)
    		d1 = delta(2)
		else
    		m3 = m(2)
    		m1 = m(1)
    		d3 = delta(2)
    		d1 = delta(1)
    	end if
		m2 = m(3)
    	d2 = delta(3)
 
 		if (m2 < m1) then
    		tmp = m1
    		m1 = m2
    		m2 = tmp
    	
        	tmp = d1
        	d1 = d2
        	d2 = tmp    	
  		else if (m2 > m3) then
    		tmp = m3
    		m3 = m2
    		m2 = tmp
    	
        	tmp = d3
    		d3 = d2
    		d2 = tmp
  		end if
  	
  	
		q = q+max(0.d0,-ms(1))+max(0.d0,-ms(2))+max(0.d0,-ms(3))
		q=q/m_par
  	
  		!check valid values of q
  		if ((q<0.d0).OR.(q>1.d0)) then
  			call mpiAbort('invalid value of q found in vArea')
  		end if
  	
  		qH = min(q,1.d0-q)
		m12=m1+m2
		m4th=min(m12,m3)

		!coordinate transformation + scaling back
		m1d = m1/d1
		m2d = m2/d2
		m3d = m3/d3
		i1 = (/-m2d,m1d,0.d0/)*(1.d0/(sqrt(m1d*m1d+m2d*m2d)))
		i2 = (/m1d,m2d,m3d/)*(1.d0/sqrt(m1d*m1d+m2d*m2d+m3d*m3d))
		i3(1) = i1(2)*i2(3)-i1(3)*i2(2)
		i3(2) = i1(3)*i2(1)-i1(1)*i2(3)
		i3(3) = i1(1)*i2(2)-i1(2)*i2(1)

		!check for particular case: m1=m2=0
		if ( (m1==0.d0).AND.(m2==0.d0) ) then
    		area = d1*d2
    		return
    	end if

		if (qH<=m1) then
		
    		p1=(/qH/m1d,0.d0,0.d0/)
    		p2=(/0.d0,qH/m2d,0.d0/)
    		p3=(/0.d0,0.d0,qH/m3d/)
    	
			p1p(1) = i1(1)*p1(1) + i1(2)*p1(2) + i1(3)*p1(3)
			p1p(2) = i2(1)*p1(1) + i2(2)*p1(2) + i2(3)*p1(3)
			p1p(3) = i3(1)*p1(1) + i3(2)*p1(2) + i3(3)*p1(3)
		
			p2p(1) = i1(1)*p2(1) + i1(2)*p2(2) + i1(3)*p2(3)
			p2p(2) = i2(1)*p2(1) + i2(2)*p2(2) + i2(3)*p2(3)
			p2p(3) = i3(1)*p2(1) + i3(2)*p2(2) + i3(3)*p2(3)
		
			p3p(1) = i1(1)*p3(1) + i1(2)*p3(2) + i1(3)*p3(3)
			p3p(2) = i2(1)*p3(1) + i2(2)*p3(2) + i2(3)*p3(3)
			p3p(3) = i3(1)*p3(1) + i3(2)*p3(2) + i3(3)*p3(3)
		
    		area = 0.5d0*abs(p1p(1)*p2p(3)-p2p(1)*p1p(3)	&
                      	+ p2p(1)*p3p(3)-p3p(1)*p2p(3)		&
                      	+ p3p(1)*p1p(3)-p1p(1)*p3p(3))

		else if (qH<=m2) then

    		p1=(/d1,(qH-m1)/m2d,0.d0/)
    		p2=(/0.d0,qH/m2d,0.d0/)
    		p3=(/0.d0,0.d0,qH/m3d/)
    		p4=(/d1,0.d0,(qH-m1)/m3d/)
    	
			p1p(1) = i1(1)*p1(1) + i1(2)*p1(2) + i1(3)*p1(3)
			p1p(2) = i2(1)*p1(1) + i2(2)*p1(2) + i2(3)*p1(3)
			p1p(3) = i3(1)*p1(1) + i3(2)*p1(2) + i3(3)*p1(3)
		
			p2p(1) = i1(1)*p2(1) + i1(2)*p2(2) + i1(3)*p2(3)
			p2p(2) = i2(1)*p2(1) + i2(2)*p2(2) + i2(3)*p2(3)
			p2p(3) = i3(1)*p2(1) + i3(2)*p2(2) + i3(3)*p2(3)
		
			p3p(1) = i1(1)*p3(1) + i1(2)*p3(2) + i1(3)*p3(3)
			p3p(2) = i2(1)*p3(1) + i2(2)*p3(2) + i2(3)*p3(3)
			p3p(3) = i3(1)*p3(1) + i3(2)*p3(2) + i3(3)*p3(3)
		
			p4p(1) = i1(1)*p4(1) + i1(2)*p4(2) + i1(3)*p4(3)
			p4p(2) = i2(1)*p4(1) + i2(2)*p4(2) + i2(3)*p4(3)
			p4p(3) = i3(1)*p4(1) + i3(2)*p4(2) + i3(3)*p4(3)
		
    		area = 0.5d0*abs(p1p(1)*p2p(3)-p2p(1)*p1p(3) &
                   	  + p2p(1)*p3p(3)-p3p(1)*p2p(3)		 &
                      + p3p(1)*p4p(3)-p4p(1)*p3p(3)      &
                      + p4p(1)*p1p(3)-p1p(1)*p4p(3))   
                                  
		else if (qH<=m4th) then

    		p1=(/d1,(qH-m1)/m2d,0.d0/)
    		p2=(/(qH-m2)/m1d,d2,0.d0/)
    		p3=(/0.d0,d2,(qH-m2)/m3d/)
    		p4=(/0.d0,0.d0,qH/m3d/)
    		p5=(/d1,0.d0,(qH-m1)/m3d/)

			p1p(1) = i1(1)*p1(1) + i1(2)*p1(2) + i1(3)*p1(3)
			p1p(2) = i2(1)*p1(1) + i2(2)*p1(2) + i2(3)*p1(3)
			p1p(3) = i3(1)*p1(1) + i3(2)*p1(2) + i3(3)*p1(3)
		
			p2p(1) = i1(1)*p2(1) + i1(2)*p2(2) + i1(3)*p2(3)
			p2p(2) = i2(1)*p2(1) + i2(2)*p2(2) + i2(3)*p2(3)
			p2p(3) = i3(1)*p2(1) + i3(2)*p2(2) + i3(3)*p2(3)
		
			p3p(1) = i1(1)*p3(1) + i1(2)*p3(2) + i1(3)*p3(3)
			p3p(2) = i2(1)*p3(1) + i2(2)*p3(2) + i2(3)*p3(3)
			p3p(3) = i3(1)*p3(1) + i3(2)*p3(2) + i3(3)*p3(3)
		
			p4p(1) = i1(1)*p4(1) + i1(2)*p4(2) + i1(3)*p4(3)
			p4p(2) = i2(1)*p4(1) + i2(2)*p4(2) + i2(3)*p4(3)
			p4p(3) = i3(1)*p4(1) + i3(2)*p4(2) + i3(3)*p4(3)
		
			p5p(1) = i1(1)*p5(1) + i1(2)*p5(2) + i1(3)*p5(3)
			p5p(2) = i2(1)*p5(1) + i2(2)*p5(2) + i2(3)*p5(3)
			p5p(3) = i3(1)*p5(1) + i3(2)*p5(2) + i3(3)*p5(3)
		
    		area = 0.5d0*abs(p1p(1)*p2p(3)-p2p(1)*p1p(3)	&
                   			+ p2p(1)*p3p(3)-p3p(1)*p2p(3)	&
                   			+ p3p(1)*p4p(3)-p4p(1)*p3p(3)	&
                   			+ p4p(1)*p5p(3)-p5p(1)*p4p(3)	&
                   			+ p5p(1)*p1p(3)-p1p(1)*p5p(3)) 
                   			
		else if (m12<m3) then
		
    		p1=(/d1,d2,(qH-m1-m2)/m3d/)
    		p2=(/0.d0,d2,(qH-m2)/m3d/)
    		p3=(/0.d0,0.d0,qH/m3d/)
    		p4=(/d1,0.d0,(qH-m1)/m3d/)
		
			p1p(1) = i1(1)*p1(1) + i1(2)*p1(2) + i1(3)*p1(3)
			p1p(2) = i2(1)*p1(1) + i2(2)*p1(2) + i2(3)*p1(3)
			p1p(3) = i3(1)*p1(1) + i3(2)*p1(2) + i3(3)*p1(3)
		
			p2p(1) = i1(1)*p2(1) + i1(2)*p2(2) + i1(3)*p2(3)
			p2p(2) = i2(1)*p2(1) + i2(2)*p2(2) + i2(3)*p2(3)
			p2p(3) = i3(1)*p2(1) + i3(2)*p2(2) + i3(3)*p2(3)
		
			p3p(1) = i1(1)*p3(1) + i1(2)*p3(2) + i1(3)*p3(3)
			p3p(2) = i2(1)*p3(1) + i2(2)*p3(2) + i2(3)*p3(3)
			p3p(3) = i3(1)*p3(1) + i3(2)*p3(2) + i3(3)*p3(3)
		
			p4p(1) = i1(1)*p4(1) + i1(2)*p4(2) + i1(3)*p4(3)
			p4p(2) = i2(1)*p4(1) + i2(2)*p4(2) + i2(3)*p4(3)
			p4p(3) = i3(1)*p4(1) + i3(2)*p4(2) + i3(3)*p4(3)
		
   			area = 0.5d0*abs(p1p(1)*p2p(3)-p2p(1)*p1p(3)	&
                   			+ p2p(1)*p3p(3)-p3p(1)*p2p(3)	&
                   			+ p3p(1)*p4p(3)-p4p(1)*p3p(3)	&
                   			+ p4p(1)*p1p(3)-p1p(1)*p4p(3))
		else

    		p1=(/d1,(qH-m1)/m2d,0.d0/)
    		p2=(/(qH-m2)/m1d,d2,0.d0/)
    		p3=(/0.d0,d2,(qH-m2)/m3d/)
    		p4=(/0.d0,(qH-m3)/m2d,d3/)
    		p5=(/(qH-m3)/m1d,0.d0,d3/)
    		p6=(/d1,0.d0,(qH-m1)/m3d/)
    	
			p1p(1) = i1(1)*p1(1) + i1(2)*p1(2) + i1(3)*p1(3)
			p1p(2) = i2(1)*p1(1) + i2(2)*p1(2) + i2(3)*p1(3)
			p1p(3) = i3(1)*p1(1) + i3(2)*p1(2) + i3(3)*p1(3)
		
			p2p(1) = i1(1)*p2(1) + i1(2)*p2(2) + i1(3)*p2(3)
			p2p(2) = i2(1)*p2(1) + i2(2)*p2(2) + i2(3)*p2(3)
			p2p(3) = i3(1)*p2(1) + i3(2)*p2(2) + i3(3)*p2(3)
		
			p3p(1) = i1(1)*p3(1) + i1(2)*p3(2) + i1(3)*p3(3)
			p3p(2) = i2(1)*p3(1) + i2(2)*p3(2) + i2(3)*p3(3)
			p3p(3) = i3(1)*p3(1) + i3(2)*p3(2) + i3(3)*p3(3)
		
			p4p(1) = i1(1)*p4(1) + i1(2)*p4(2) + i1(3)*p4(3)
			p4p(2) = i2(1)*p4(1) + i2(2)*p4(2) + i2(3)*p4(3)
			p4p(3) = i3(1)*p4(1) + i3(2)*p4(2) + i3(3)*p4(3)
		
			p5p(1) = i1(1)*p5(1) + i1(2)*p5(2) + i1(3)*p5(3)
			p5p(2) = i2(1)*p5(1) + i2(2)*p5(2) + i2(3)*p5(3)
			p5p(3) = i3(1)*p5(1) + i3(2)*p5(2) + i3(3)*p5(3)

			p6p(1) = i1(1)*p6(1) + i1(2)*p6(2) + i1(3)*p6(3)
			p6p(2) = i2(1)*p6(1) + i2(2)*p6(2) + i2(3)*p6(3)
			p6p(3) = i3(1)*p6(1) + i3(2)*p6(2) + i3(3)*p6(3)
		
    		area = 0.5d0*abs(p1p(1)*p2p(3)-p2p(1)*p1p(3)	&
                   			+ p2p(1)*p3p(3)-p3p(1)*p2p(3)	&
                   			+ p3p(1)*p4p(3)-p4p(1)*p3p(3)	&
                   			+ p4p(1)*p5p(3)-p5p(1)*p4p(3)	&
                   			+ p5p(1)*p6p(3)-p6p(1)*p5p(3) 	&
                   			+ p6p(1)*p1p(3)-p1p(1)*p6p(3)) 
		end if
  		 
    end subroutine
!========================================================================================!

!========================================================================================!
    subroutine centroid(d,c,n,q,ct,failed)
    	real(DP), dimension(3), intent(in) :: d,c
        real(DP), dimension(3), intent(in) :: n
        real(DP), intent(in) :: q
        real(DP), dimension(3), intent(out) :: ct
        logical, intent(out) :: failed
        real(DP), dimension(3,12) :: v
        integer :: nv
        
        !note: the centroid is computed by arithmetic average of the vertices
        call intersPoints(d,c,n,q,v,nv)
        
        if (nv < 3) then
        	!call mpiABORT('Valid centroid not found ')
        	ct(1) = 0.d0
        	ct(2) = 0.d0
        	ct(3) = 0.d0
        	failed = .TRUE.
        	return
        end if
        
        ct(1) = sum(v(1,1:nv))/nv
        ct(2) = sum(v(2,1:nv))/nv
        ct(3) = sum(v(3,1:nv))/nv
        
        failed = .FALSE.
	
    end subroutine
!========================================================================================!

!========================================================================================!
    subroutine intersPoints(d,c,n,q,v,nv)
    	real(DP), dimension(3), intent(in) :: d,c
        real(DP), dimension(3), intent(in) :: n
        real(DP), intent(in) :: q
        real(DP), dimension(3,12), intent(out) :: v
        integer, intent(out) :: nv
        real(DP), dimension(3) :: p
        real(DP) :: dx, dy, dz
        real(DP), dimension(3) :: p1,p2,p3,p4,p5,p6,p7,p8
        integer :: l,m
        
        v = 0.d0
        p = 0.d0
        
        dx = d(1)
        dy = d(2)
        dz = d(3)
		
		p1 = (/0.d0,0.d0,0.d0/)
		p2 = (/dx,0.d0,0.d0/)
		p3 = (/dx,dy,0.d0/)
		p4 = (/0.d0,dy,0.d0/)
		p5 = (/0.d0,0.d0,dz/)
		p6 = (/dx,0.d0,dz/)
		p7 = (/dx,dy,dz/)
		p8 = (/0.d0,dy,dz/)
	
		l = 1
		if (intersPlaneSegm(n,q,p1,p2,p,1)) then
			v(:,l) = p
			l = l+1
		end if
		if (intersPlaneSegm(n,q,p2,p3,p,2)) then
			v(:,l) = p
			l = l+1
		end if
		if (intersPlaneSegm(n,q,p3,p4,p,1)) then
			v(:,l) = p
			l = l+1
		end if
		if (intersPlaneSegm(n,q,p4,p1,p,2)) then
			v(:,l) = p
			l = l+1
		end if
		if (intersPlaneSegm(n,q,p5,p6,p,1)) then
			v(:,l) = p
			l = l+1
		end if
		if (intersPlaneSegm(n,q,p6,p7,p,2)) then
			v(:,l) = p
			l = l+1
		end if
		if (intersPlaneSegm(n,q,p7,p8,p,1)) then
			v(:,l) = p
			l = l+1
		end if
		if (intersPlaneSegm(n,q,p8,p5,p,2)) then
			v(:,l) = p
			l = l+1
		end if
		if (intersPlaneSegm(n,q,p1,p5,p,3)) then
			v(:,l) = p
			l = l+1
		end if
		if (intersPlaneSegm(n,q,p2,p6,p,3)) then
			v(:,l) = p
			l = l+1
		end if
		if (intersPlaneSegm(n,q,p3,p7,p,3)) then
			v(:,l) = p
			l = l+1
		end if
		if (intersPlaneSegm(n,q,p4,p8,p,3)) then
			v(:,l) = p
			l = l+1	
		end if

	!global coord system
	do m=1,l-1
		v(1,m) = v(1,m) + c(1)
		v(2,m) = v(2,m) + c(2)
		v(3,m) = v(3,m) + c(3)
	end do
	
	nv = l-1

		
    end subroutine
!========================================================================================!

!========================================================================================!
    function intersPlaneSegm(m,q,p1,p2,pv,sdir) result(isInters)
        real(DP), dimension(3), intent(in) :: m
        real(DP), intent(in) :: q
        real(DP), dimension(3), intent(in) :: p1, p2
        integer, intent(in) :: sdir
        real(DP), dimension(3), intent(out) :: pv
        real(DP), dimension(3) :: t
        real(DP) :: s
        logical :: isInters
        

		t = 0.d0
		pv = 0.d0
		
		if (m(sdir) == 0.d0) then
			isInters = .FALSE.
			return
		else
			s = (q-m(1)*p1(1)-m(2)*p1(2)-m(3)*p1(3))/(m(sdir)*(p2(sdir)-p1(sdir)))
    		t(sdir)=s
    		if ( (s < 0.d0) .OR. ( s > 1.d0 ) ) then
    			isInters = .FALSE.
       			return
    		else 
    			pv(1)=p1(1)*(1.d0-t(1))+p2(1)*t(1)
    			pv(2)=p1(2)*(1.d0-t(2))+p2(2)*t(2)
    			pv(3)=p1(3)*(1.d0-t(3))+p2(3)*t(3)
    		end if
    	end if
	
		isInters = .TRUE.
		
    end function
!========================================================================================!

!************************************* helpers ******************************************!

!========================================================================================!
    subroutine setColumnPar(m,dir,swi,swj,swk,it,jt,kt,ib,jb,kb)
        integer, intent(in) :: dir
        real(DP),intent(in) :: m
        integer, intent(out) :: swi,swj,swk,it,jt,kt,ib,jb,kb
        
        !switch for delta cells
		select case(dir)
        	case(1)
        		swi = 1
        		swj = 0
        		swk = 0
        	case(2)
        		swi = 0
        		swj = 1
        		swk = 0
        	case(3)
        		swi = 0
        		swj = 0
        		swk = 1
        	case default
        end select
        
        !set search indexes 
        if (m >= 0.d0) then
        	select case(dir)
        		case(1)
        			it = 1
        			jt = 0
        			kt = 0 
        			ib = -1
        			jb = 0
        			kb = 0 
        		case(2)
        			it = 0
        			jt = 1
        			kt = 0 
        			ib = 0
        			jb = -1
        			kb = 0 
        		case(3)
        			it = 0
        			jt = 0
        			kt = 1 
        			ib = 0
        			jb = 0
        			kb = -1 
        		case default
        	end select
        else
        	select case(dir)
        		case(1)
        			it = -1
        			jt = 0
        			kt = 0 
        			ib = 1
        			jb = 0
        			kb = 0 
        		case(2)
        			it = 0
        			jt = -1
        			kt = 0
        			ib = 0
        			jb = 1
        			kb = 0  
        		case(3)
        			it = 0
        			jt = 0
        			kt = -1 
        			ib = 0
        			jb = 0
        			kb = 1 
        		case default
        	end select
        end if

		     
		
		
    end subroutine
!========================================================================================!

!========================================================================================!
    subroutine setStencilPar(dir,d)
        integer, intent(in) :: dir
        integer, intent(out), dimension(3,9) :: d
        
        d = 0
        
		select case(dir)
        	case(1)
        		!delta j
        		d(2,1) = -1
        		d(2,3) =  1
        		d(2,4) = -1
        		d(2,6) =  1
        		d(2,7) = -1
        		d(2,9) = 1
        		!delta k
        		d(3,1:3) = -1
        		d(3,7:9) = 1		
        	case(2)
        		!delta i
        		d(1,1) = -1
        		d(1,3) =  1
        		d(1,4) = -1
        		d(1,6) =  1
        		d(1,7) = -1
        		d(1,9) = 1
        		!delta k
        		d(3,1:3) = -1
        		d(3,7:9) = 1
        	case(3)
        		!delta i
        		d(1,1) = -1
        		d(1,3) =  1
        		d(1,4) = -1
        		d(1,6) =  1
        		d(1,7) = -1
        		d(1,9) = 1
        		!delta j
        		d(2,1:3) = -1
        		d(2,7:9) = 1
        	case default
        end select
		     
		
		
    end subroutine
!========================================================================================!

!========================================================================================!
    subroutine setDirections(dir,i,j,k,nc,n1,n2,dir1,dir2)
        integer, intent(in) :: dir, i, j, k
        integer, intent(out) :: nc,n1,n2,dir1,dir2
		
		select case(dir)
			case(1)
				nc = i
				dir1 = 2
				dir2 = 3
				n1 = j
				n2 = k
			case(2)
				nc = j
				dir1 = 1
				dir2 = 3
				n1 = i
				n2 = k
			case(3)
				nc = k
				dir1 = 1
				dir2 = 2
				n1 = i
				n2 = j
			case default
		end select
		

    end subroutine
!========================================================================================!

!========================================================================================!
    subroutine n_interfaces(this)
    	type(VOF), intent(in) :: this
    	type(mpiControl), pointer :: comm
    	integer :: b,c_l,c_g,i,j,k,lbi,ubi,lbj,ubj,lbk,ubk,ierror
 
 		comm => this%mesh_%ptrMPIC_
    	c_l=0

		do b=1,s_nblk	

			lbi = lbound(vofBlocks(b)%c,1)
			ubi = ubound(vofBlocks(b)%c,1)
			lbj = lbound(vofBlocks(b)%c,2)
			ubj = ubound(vofBlocks(b)%c,2)
			lbk = lbound(vofBlocks(b)%c,3)
			ubk = ubound(vofBlocks(b)%c,3)	
					
				do k=lbk,ubk
					do j=lbj,ubj
						do i=lbi,ubi
						
							if (vofBlocks(b)%isMixed(i,j,k)) then
							
								c_l=c_l+1
							
							end if
						
						end do
					end do
				end do	
		end do
		
		call Mpi_Reduce(c_l, c_g, 1, MPI_INTEGER, MPI_SUM, 0, comm%cartComm_, ierror)
		
		if (IS_MASTER) then
			write(*,*) 'Total INTERFACES: ', c_g
		end if
		
		
    end subroutine
!========================================================================================!	

!========================================================================================!
    subroutine max_box_size(this)
    	type(VOF), intent(in) :: this
    	type(mpiControl), pointer :: comm
    	integer :: b,is,ie,js,je,ks,ke,sx,sy,sz,sx_g,sy_g,sz_g,ierror
 
 		comm => this%mesh_%ptrMPIC_

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

end module VOFMod


	



