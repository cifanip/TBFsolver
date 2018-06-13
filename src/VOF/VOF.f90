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

module vofMOD

	use vofBlocksMod
	use timeMod
	
	integer, parameter :: s_Xsweep = 1
	integer, parameter :: s_Ysweep = 2
	integer, parameter :: s_Zsweep = 3
	
	integer, protected :: s_sweepSelector = 0
	integer, dimension(3), protected :: s_sweep
	
	type, public :: VOF
	
		type(grid), pointer :: mesh_ => NULL()
		type(grid), pointer :: gmesh_ => NULL()
		type(time), pointer :: ptrTime_ => NULL() 
		
		!threshold vof field
		real(DP), private :: eps_ = 1.d-8
		
		!halo dim of box volume fraction field
		integer, private :: hd_
		
		!mat props
		real(DP) :: rhol_, mul_, rhog_, mug_
		real(DP) :: sigma_	
		
	end type
	
	!vertex field
	real(DP), allocatable, dimension(:,:,:), protected :: s_cv
	

	public :: vofCTOR
	public :: updateMaterialProps
	public :: solveVOF
	public :: computeSurfaceTension
	
	private :: smoothVF
	private :: cnh
	private :: cnp
	private :: storeOldVOField
	private :: reconstruct
	private :: reconstruct_blocks
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
	private :: resetFragments
	private :: computeNormal_sharp
	private :: computeNormal_smooth
	private :: computeNormal_youngs
	private :: computeNormal_hf
	private :: hfNormal
	private :: smoothVFblock
	private :: correctContantAngle
	private :: computeBlockST
	private :: computeCurvature
	private :: spreadCurvature
	private :: computeBlockCurvature
	private :: hfColumn
	private :: hfCurvature
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
	
contains

!========================================================================================!
    subroutine vofCTOR(this,gmesh,mesh,rt,flow_solver,tpf)
        type(VOF), intent(out) :: this
        type(time), intent(in), target :: rt
        type(grid), intent(in), target :: gmesh, mesh
        integer, intent(in) :: flow_solver,tpf
        integer :: nprocs
        type(parFile) :: pfile
        integer :: nx,ny,nz
        
        this%ptrTime_ => rt
        this%mesh_ => mesh
        this%gmesh_ => gmesh
        
        nx=mesh%nx_
        ny=mesh%ny_
        nz=mesh%nz_
		
		this%hd_ = 3
		
		!read mat props
		call parFileCTOR(pfile ,'parameters','specs')
		call readParameter(pfile,this%rhol_,'rhol')
		call readParameter(pfile,this%rhog_,'rhog')
		call readParameter(pfile,this%mul_,'mul')
		call readParameter(pfile,this%mug_,'mug')
		call readParameter(pfile,this%sigma_,'sigma')
		
		if (flow_solver==tpf) then
		
			!init boxes
			call vofBlocksCTOR(mesh,gmesh)
		
    		!allocate vertex field
    		call allocateArray(s_cv,0,nx,0,ny,0,nz)
    	
			!store old c field
    		call storeOldVOField()

			call updateStateBlocks(this)
			
		end if
		
    end subroutine
!========================================================================================!

!========================================================================================!
    subroutine updateMaterialProps(this,c,cs,rho,mu)
    	type(VOF), intent(in) :: this
    	type(field), intent(in) :: c
    	type(field), intent(inout) :: cs,rho,mu
        integer :: lbi, ubi, lbj, ubj, lbk, ubk
        integer :: i,j,k
       
        !all the mat props have the same box size 
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

!========================================================================================!
    subroutine smoothVF(this,c,cs)
        type(VOF), intent(in) :: this
        type(field), intent(in) :: c
        type(field), intent(inout) :: cs
        
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
    subroutine storeOldVOField()
    	integer :: b
    	
		do b=1,s_nblk
			vofBlocks(b)%c0=vofBlocks(b)%c
		end do  		
    
    end subroutine
!========================================================================================!

!========================================================================================!
    subroutine solveVOF(this,c,u)
    	type(VOF), intent(inout) :: this
    	type(field), intent(inout) :: c
    	type(vfield), intent(in) :: u
    	integer :: i,b
    	real(DP) :: start, finish
    	
    	start = MPI_Wtime()
    	
    	!set c=c^{n+1}
    	call cnp()

		call sweepCombination()
		
		!exchange velocity field
		call grid_2_boxes_u(this%mesh_,u)

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
			call updateBlock(this%mesh_,this%gmesh_,vofBlocks(b),b)
			call updateState(this,vofBlocks(b))
			call reconstruct(this,vofBlocks(b))
		
		end do
		!$OMP END PARALLEL DO
	
		call gatherLogicalExchange(this%mesh_)
		
		!update list
		if (any(s_exchange_g)) then
			call excLists(this%mesh_)
		end if
	
		!exchange volume fraction field
		call boxes_2_grid_f(this%mesh_,c,PACK_BOX_VF,UNPACK_MAX)
		
		!redistribute blocks
		if (vofBlocksRed(this%ptrTime_)) then
			if (IS_MASTER) then
				write(*,*) '	Blocks redistribution called'
			end if
			call reInitBlockDistribution(this%mesh_,this%gmesh_)	
			call updateStateBlocks(this)
			call reconstruct_blocks(this)
		end if	

		finish = MPI_Wtime()
		
		call info(this,finish-start,0)
    	
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
    subroutine reconstruct_blocks(this)
    	type(VOF), intent(in) :: this
    	integer :: b
    	
    	do b=1,s_nblk
    		call reconstruct(this,vofBlocks(b))
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
	
	vofb%corrTerm = 0.d0
    	    	
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
        integer :: b
        
        do b=1,s_nblk
        
        	call updateState(this,vofBlocks(b))
		
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
							call hfColumn(vofb,mmax,mloc,i+d(1,cn),j+d(2,cn),k+d(3,cn),cn,h,lb,pos_hf,isValid)
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
    	type(vfield), intent(inout) :: st
    	type(field), intent(inout) :: k
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
		
		call boxes_2_grid_vf(this%mesh_,st,PACK_BOX_ST,UNPACK_SUM)
		
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
    	type(field), intent(inout) :: k
    	integer :: b

		!$OMP PARALLEL DO DEFAULT(none) &
		!$OMP SHARED(vofBlocks,s_nblk) &
		!$OMP PRIVATE(b)    	
    	do b=1,s_nblk
    		call computeBlockCurvature(vofBlocks(b))
    		!call spreadCurvature(vofBlocks(b))
    	end do
    	!$OMP END PARALLEL DO
    
    	call boxes_2_grid_f(this%mesh_,k,PACK_BOX_K,UNPACK_SUM)

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
    subroutine computeBlockCurvature(vofb)
! ************************************************************************************** !
! This subroutine is an implementation of the GHF method developed by:
! Popinet, S. (2009). An accurate adaptive solver for surface-tension-driven 
! interfacial flows. Journal of Computational Physics, 228(16), 5838-5866.
! ************************************************************************************** !
        type(vofBlock), intent(inout) :: vofb
        integer :: i, j, k
        real(DP) :: mx, my, mz
        real(DP) :: m_1,m_2,m_3,small
        integer :: mloc_1,mloc_2,mloc_3
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
				
					if (vofb%isMixed(i,j,k)) then
					
					mx = vofb%nx(i,j,k)
					my = vofb%ny(i,j,k)
					mz = vofb%nz(i,j,k)
					
					mv = (/ mx, my, mz /)
					small=tiny(1.d0)
					!sort
					mloc_1 = maxloc( abs(mv),1 )
					m_1 = mv(mloc_1)
					mv(mloc_1)=small
					mloc_2 = maxloc( abs(mv),1 )
					m_2 = mv(mloc_2)
					mv(mloc_2)=small
					mloc_3 = maxloc( abs(mv),1 )
					m_3 = mv(mloc_3)
					
					!search dir 1
					call setStencilPar(mloc_1,d)
					hfk = .TRUE.
					do cn=1,9
						call hfColumn(vofb,m_1,mloc_1,i+d(1,cn),j+d(2,cn),k+d(3,cn),cn,h,lb,pos_hf,isValid)
						if (.NOT. (isValid(cn))) then
							hfk = .FALSE.
						end if						
					end do	
					isBlockValid(1:9)=isValid
					pos_block(:,1:9)=pos_hf		
					if (hfk) then
						call hfCurvature(vofb,mloc_1,m_1,i,j,k,h,lb)	
						cycle	
					end if	
					
					!search dir 2
					call setStencilPar(mloc_2,d)
					hfk = .TRUE.
					do cn=1,9
						call hfColumn(vofb,m_2,mloc_2,i+d(1,cn),j+d(2,cn),k+d(3,cn),cn,h,lb,pos_hf,isValid)
						if (.NOT. (isValid(cn))) then
							hfk = .FALSE.
						end if						
					end do	
					isBlockValid(10:18)=isValid
					pos_block(:,10:18)=pos_hf	
					if (hfk) then
						call hfCurvature(vofb,mloc_2,m_2,i,j,k,h,lb)
						cycle	
					end if
					
					!search dir 3
					call setStencilPar(mloc_3,d)
					hfk = .TRUE.
					do cn=1,9
						call hfColumn(vofb,m_3,mloc_3,i+d(1,cn),j+d(2,cn),k+d(3,cn),cn,h,lb,pos_hf,isValid)
						if (.NOT. (isValid(cn))) then
							hfk = .FALSE.
						end if						
					end do	
					isBlockValid(19:27)=isValid
					pos_block(:,19:27)=pos_hf	
					if (hfk) then
						call hfCurvature(vofb,mloc_3,m_3,i,j,k,h,lb)	
						cycle	
					end if						

					!standard HF failed (execute parabolic fit)
					
					call parabFittedCurvature(vofb,i,j,k,pos_block,isBlockValid,.TRUE.,failed)	
					
					if (failed) then
						call interCentroids(vofb,i,j,k,pos_block,isBlockValid)
						call parabFittedCurvature(vofb,i,j,k,pos_block,isBlockValid,.FALSE.,failed)						
					end if
							
					!	call fdCurvature(vofb,i,j,k)	
						
					end if
	
					
				end do
			end do
		end do      
		
		
    end subroutine
!========================================================================================!

!========================================================================================!
    subroutine hfColumn(vofb,m,dir,i,j,k,cn,h,lb,posi,isValid)
        type(vofBlock), intent(in) :: vofb
        real(DP),intent(in) :: m
        integer, intent(in) :: dir, i, j ,k, cn
        real(DP), dimension(9), intent(inout) :: h
        integer, dimension(9), intent(inout) :: lb
        logical, dimension(9), intent(inout) :: isValid
        real(DP), dimension(3,9), intent(inout) :: posi
        integer :: swi,swj,swk,it,jt,kt,ib,jb,kb
        integer :: ii, jj, kk, l, idx_lb, idx_ub
        logical :: search,isEmpty,isMixed,isFull
        real(DP), pointer, dimension(:) :: df => NULL()
        real(DP), pointer, dimension(:) :: posf => NULL()
        !$omp threadprivate(df,posf)


		!chiama set par      
		call setColumnPar(m,dir,swi,swj,swk,it,jt,kt,ib,jb,kb)
		
		select case(dir)
			case(1)
				df => vofb%dxf
				posf => vofb%xf
				idx_ub=vofb%idx(2)+offset_c-1
				idx_lb=vofb%idx(1)-offset_c+1
			case(2)
				df => vofb%dyf
				posf => vofb%yf
				idx_ub=vofb%idx(4)+offset_c-1
				idx_lb=vofb%idx(3)-offset_c+1
			case(3)
				df => vofb%dzf
				posf => vofb%zf
				idx_ub=vofb%idx(6)+offset_c-1
				idx_lb=vofb%idx(5)-offset_c+1
			case default
		end select
		
		h(cn) = vofb%c(i,j,k)*df(swi*i+swj*j+swk*k)
		
		!search the top column
		isEmpty=((.NOT.(vofb%isFull(i,j,k))).AND.&
				 (.NOT.(vofb%isMixed(i,j,k))))
		isMixed=vofb%isMixed(i,j,k)
		isFull=vofb%isFull(i,j,k)
		if (.not.isFull) then
			search=.FALSE.
		else
			search=.TRUE.
		end if
		l=1
		do while (search.OR.isMixed)
			ii = i+l*it
			jj = j+l*jt
			kk = k+l*kt
			h(cn) = h(cn) + vofb%c(ii,jj,kk)*df(swi*ii+swj*jj+swk*kk)
			
			if (vofb%isMixed(ii,jj,kk)) then
				search=.FALSE.
				isMixed=.TRUE.
			else
				isMixed=.FALSE.
			end if
			
			isEmpty=((.NOT.(vofb%isFull(ii,jj,kk))).AND.&
				     (.NOT.(vofb%isMixed(ii,jj,kk))))
			l=l+1
			
			select case(dir)
				case(1)
					if ((ii>=idx_ub).OR.(ii<=idx_lb)) then
						isValid(cn) =.FALSE.
						posi(:,cn)=0.d0
						return							
					end if
				case(2)
					if ((jj>=idx_ub).OR.(jj<=idx_lb)) then
						isValid(cn) =.FALSE.
						posi(:,cn)=0.d0
						return							
					end if
				case(3)
					if ((kk>=idx_ub).OR.(kk<=idx_lb)) then
						isValid(cn) =.FALSE.
						posi(:,cn)=0.d0
						return							
					end if
				case default
			end select
			
		end do

		if (isEmpty) then
			isValid(cn) = .TRUE.
		else
			isValid(cn) =.FALSE.
			posi(:,cn)=0.d0
			return		
		end if
		

		!search the bottom column
		isEmpty=((.NOT.(vofb%isFull(i,j,k))).AND.&
				 (.NOT.(vofb%isMixed(i,j,k))))
		isMixed=vofb%isMixed(i,j,k)
		isFull=vofb%isFull(i,j,k)
		if (.not.isEmpty) then
			search=.FALSE.
		else
			search=.TRUE.
		end if
		l=1
		do while (search.OR.isMixed)
			ii = i+l*ib
			jj = j+l*jb
			kk = k+l*kb
			h(cn) = h(cn) + vofb%c(ii,jj,kk)*df(swi*ii+swj*jj+swk*kk)
			
			if (vofb%isMixed(ii,jj,kk)) then
				search=.FALSE.
				isMixed=.TRUE.
			else
				isMixed=.FALSE.
			end if
			
			isFull=vofb%isFull(ii,jj,kk)
			
			l=l+1
				
			select case(dir)
				case(1)
					if ((ii>=idx_ub).OR.(ii<=idx_lb)) then
						isValid(cn) =.FALSE.
						posi(:,cn)=0.d0
						return							
					end if
				case(2)
					if ((jj>=idx_ub).OR.(jj<=idx_lb)) then
						isValid(cn) =.FALSE.
						posi(:,cn)=0.d0
						return							
					end if
				case(3)
					if ((kk>=idx_ub).OR.(kk<=idx_lb)) then
						isValid(cn) =.FALSE.
						posi(:,cn)=0.d0
						return							
					end if
				case default
			end select

		end do
		
		if (isFull) then
			isValid(cn) = .TRUE.
			lb(cn) = l-1
		else
			isValid(cn) =.FALSE.
			posi(:,cn)=0.d0
			return		
		end if

		!store interface position
		posi(1,cn)=vofb%xc(i)
		posi(2,cn)=vofb%yc(j)
		posi(3,cn)=vofb%zc(k)
		if (m >= 0.d0) then
			posi(dir,cn)=posf(swi*i+swj*j+swk*k-lb(cn)-1)+h(cn)
		else
			posi(dir,cn)=posf(swi*i+swj*j+swk*k+lb(cn))-h(cn)
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
! This subroutine is an implementation of the analytical relations developed by:
! Scardovelli, R., & Zaleski, S. (2000). Analytical relations connecting linear 
! interfaces and volume fractions in rectangular grids. 
! Journal of Computational Physics, 164(1), 228-237.
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
! This subroutine is an implementation of the analytical relations developed by:
! Scardovelli, R., & Zaleski, S. (2000). Analytical relations connecting linear 
! interfaces and volume fractions in rectangular grids. 
! Journal of Computational Physics, 164(1), 228-237.
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

end module VOFMod








