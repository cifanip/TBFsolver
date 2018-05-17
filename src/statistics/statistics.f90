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

module statisticsMod

	use vfieldMod
	
	implicit none
	
	integer, parameter :: s_whole_region = 0
	integer, parameter :: s_gas_region = 1
	integer, parameter :: s_liquid_region = 2
	
	integer, parameter :: stat_scalar = 0
	integer, parameter :: stat_vector = 1
	integer, parameter :: stat_symtensor = 2
	integer, parameter :: stat_tensor = 3
	integer, parameter :: stat_tensor_sqr = 4
	
	!averages
	!volume fraction
	real(DP), allocatable, dimension(:,:) :: cm_
	
	!whole field
	real(DP), allocatable, dimension(:,:) :: pm_,ppm_
	real(DP), allocatable, dimension(:,:) :: um_,uum_,wm_,wwm_
	
	!gas field
	real(DP), allocatable, dimension(:,:) :: pmg_,ppmg_
	real(DP), allocatable, dimension(:,:) :: umg_,uumg_,wmg_,wwmg_
	
	!liquid field
	real(DP), allocatable, dimension(:,:) :: pml_,ppml_
	real(DP), allocatable, dimension(:,:) :: uml_,uuml_,wml_,wwml_
	
	!turbulent dissipation rate
	real(DP), allocatable, dimension(:,:) :: dum_,dudum_
	real(DP), allocatable, dimension(:,:) :: dumg_,dudumg_
	real(DP), allocatable, dimension(:,:) :: duml_,duduml_
	
	!time average total wall shear stress
	real(DP) :: tauw_tav
	
	type, public :: statistics
	
		real(DP), private :: Ts_
		
		type(grid), pointer :: gmesh_ => NULL()
		type(vfield), pointer :: u_ => NULL()
		type(vfield), pointer :: w_ => NULL()
		type(field), pointer :: p_ => NULL()
		type(field), pointer :: c_ => NULL()
		type(field), pointer :: mu_ => NULL()
		
		!slice volume
		real(DP), allocatable, dimension(:) :: Vs_
		
		real(DP) :: nu_
		
		logical, private :: isTimeAverage_, hasTimeAvStarted_

		contains	
	end type
	

	private :: updateStatsVF
	private :: updateStatsFields
	private :: updateStatsGradU
	private :: updateTimeAvVF
	private :: updateTimeAvFields
	private :: updateTimeAvGradU
	private :: reduceStat
	private :: computeSliceVolume
	private :: writeStats_VF
	private :: writeStats_region
	private :: shearStress
	private :: checkTimeAverage
	private :: info
		
	public :: statisticsCTOR
	public :: updateStats
	public :: writeStats

	
contains


!========================================================================================!
	subroutine statisticsCTOR(this,u,w,p,c,mu,nu,gmesh)
		type(statistics), intent(out) :: this
		type(vfield), target, intent(in) :: u,w
		type(field), target, intent(in) :: p,c,mu
		real(DP), intent(in) :: nu
		type(grid), target, intent(in) :: gmesh
		type(parFile) :: pfile
		integer :: nyg
		
		this%p_ =>  p
		this%u_ =>  u
		this%w_ =>  w
		this%c_ =>  c
		this%mu_ => mu
		
		this%nu_ = nu
		
		if (IS_MASTER) then
			this%gmesh_ => gmesh
		end if

		call parFileCTOR(pfile,'timeControl','specs')
		call readParameter(pfile,this%Ts_,'Ts')

		nyg = p%ptrMesh_%nyg_
		
		if (IS_MASTER) then
			!temporal averages
			!volume fraction
			call allocateArray(cm_,1,1,1,nyg)
			!whole
			call allocateArray(pm_,1,1,1,nyg)
			call allocateArray(um_,1,3,1,nyg)
			call allocateArray(wm_,1,3,1,nyg)
			call allocateArray(ppm_,1,1,1,nyg)
			call allocateArray(uum_,1,6,1,nyg)
			call allocateArray(wwm_,1,6,1,nyg)
			!gas
			call allocateArray(pmg_,1,1,1,nyg)
			call allocateArray(umg_,1,3,1,nyg)
			call allocateArray(wmg_,1,3,1,nyg)
			call allocateArray(ppmg_,1,1,1,nyg)
			call allocateArray(uumg_,1,6,1,nyg)
			call allocateArray(wwmg_,1,6,1,nyg)
			!liquid
			call allocateArray(pml_,1,1,1,nyg)
			call allocateArray(uml_,1,3,1,nyg)
			call allocateArray(wml_,1,3,1,nyg)
			call allocateArray(ppml_,1,1,1,nyg)
			call allocateArray(uuml_,1,6,1,nyg)
			call allocateArray(wwml_,1,6,1,nyg)
			
			!grad U tensor
			!whole
			call allocateArray(dum_,1,9,1,nyg)
			call allocateArray(dudum_,1,12,1,nyg)
			!gas
			call allocateArray(dumg_,1,9,1,nyg)
			call allocateArray(dudumg_,1,12,1,nyg)
			!liquid
			call allocateArray(duml_,1,9,1,nyg)
			call allocateArray(duduml_,1,12,1,nyg)		

			!init to zero
			cm_ = 0.d0
		
			pm_ = 0.d0
			um_ = 0.d0
			wm_ = 0.d0
			uum_ = 0.d0
			wwm_ = 0.d0
		
			pmg_ = 0.d0
			umg_ = 0.d0
			wmg_ = 0.d0
			uumg_ = 0.d0
			wwmg_ = 0.d0
		
			pml_ = 0.d0
			uml_ = 0.d0
			wml_ = 0.d0
			uuml_ = 0.d0
			wwml_ = 0.d0
		
			dum_ = 0.d0
			dudum_ = 0.d0
		
			dumg_ = 0.d0
			dudumg_ = 0.d0
		
			duml_ = 0.d0
			duduml_ = 0.d0
			
			tauw_tav = 0.d0
			
		end if
		
		!allocate space average volume
		call allocateArray(this%Vs_,1,nyg)
		this%Vs_ = 0.d0
		call computeSliceVolume(this)
		
		this%hasTimeAvStarted_=.FALSE.

	end subroutine
!========================================================================================!

!========================================================================================!
	subroutine updateStats(this,t,dt)
		type(statistics), intent(inout) :: this
		real(DP), intent(in) :: t, dt
		real(DP) :: start,finish
		
		start = MPI_Wtime()  
		
		if (.not.this%hasTimeAvStarted_) then
			call checkTimeAverage(this,t)
		end if
		
		if (this%isTimeAverage_) then
			
			call updateStatsVF(this,t,dt)
	
			call updateStatsFields(this,t,dt,s_whole_region)
			call updateStatsFields(this,t,dt,s_gas_region)
			call updateStatsFields(this,t,dt,s_liquid_region)
		
			call updateStatsGradU(this,t,dt,s_whole_region)
			call updateStatsGradU(this,t,dt,s_gas_region)
			call updateStatsGradU(this,t,dt,s_liquid_region)
        
        end if
        
        call shearStress(this,t,dt)
        
        finish = MPI_Wtime()
        
        call info(finish-start)
	
	end subroutine
!========================================================================================!

!========================================================================================!
	subroutine updateStatsVF(this,t,dt)
		type(statistics), intent(in) :: this
		real(DP), intent(in) :: t, dt
		type(grid), pointer :: mesh
		type(mpiControl), pointer :: mpic
		real(DP), allocatable, dimension(:,:) :: cm,scm
		real(DP) :: V
		integer :: i,j,k,nx,ny,nz

		mesh => this%p_%ptrMesh_
		mpic => mesh%ptrMPIC_
		
		nx = this%p_%nx_
		ny = this%p_%ny_
		nz = this%p_%nz_

		!allocate tmp fields
		call allocateArray(cm,1,1,1,ny)
		
		!reset to zero
		cm = 0.d0
		
		!mean along wall normal planes
		!$OMP PARALLEL DO DEFAULT(none) &
		!$OMP SHARED(this,mesh) &
		!$OMP SHARED(nx,ny,nz) &
		!$OMP PRIVATE(V) &
		!$OMP PRIVATE(i,j,k) &
		!$OMP REDUCTION(+:cm)
		do k=1,nz
			do j=1,ny
				do i=1,nx
				
					V=mesh%V_(i,j,k)
					cm(1,j)=cm(1,j)+this%c_%f_(i,j,k)*V					
										
				end do
			end do
		end do
		!$OMP END PARALLEL DO

		
		!reduce mean
		call reduceStat(cm,scm,this%Vs_,mpic,stat_scalar)
		
		!update time averages
		if (IS_MASTER) then
			call updateTimeAvVF(this,t,dt,scm)
		end if

	end subroutine
!========================================================================================!

!========================================================================================!
	subroutine updateStatsFields(this,t,dt,region)
		type(statistics), intent(in) :: this
		real(DP), intent(in) :: t, dt
		integer, intent(in) :: region
		type(grid), pointer :: mesh
		type(mpiControl), pointer :: mpic
		real(DP), allocatable, dimension(:,:) :: pm,ppm,spm,sppm
		real(DP), allocatable, dimension(:,:) :: um,uum,sum,suum
		real(DP), allocatable, dimension(:,:) :: wm,wwm,swm,swwm
		real(DP) :: pv,V,cv
		real(DP) :: uxv,uyv,uzv
		real(DP) :: wxv,wyv,wzv
		integer :: i,j,k,nx,ny,nz

		mesh => this%p_%ptrMesh_
		mpic => mesh%ptrMPIC_
		
		nx = this%p_%nx_
		ny = this%p_%ny_
		nz = this%p_%nz_
		
		!allocate tmp fields
		call allocateArray(pm,1,1,1,ny)
		call allocateArray(ppm,1,1,1,ny)
		call allocateArray(um,1,3,1,ny)
		call allocateArray(uum,1,6,1,ny)
		call allocateArray(wm,1,3,1,ny)
		call allocateArray(wwm,1,6,1,ny)

		pm = 0.d0
		um = 0.d0
		wm = 0.d0
		ppm = 0.d0
		uum = 0.d0
		wwm = 0.d0
		
		!mean along wall normal planes
		!$OMP PARALLEL DO DEFAULT(none) &
		!$OMP SHARED(this,mesh,region) &
		!$OMP SHARED(nx,ny,nz) &
		!$OMP PRIVATE(V,cv,pv,uxv,uyv,uzv,wxv,wyv,wzv) &
		!$OMP PRIVATE(i,j,k) &
		!$OMP REDUCTION(+:pm,ppm,um,uum,wm,wwm)
		do k=1,nz
			do j=1,ny
				do i=1,nx
				
					V=mesh%V_(i,j,k)
				
					!volume fraction
					select case(region)
						case(s_whole_region)
							cv=1.d0
						case(s_gas_region)
							cv = this%c_%f_(i,j,k)
						case(s_liquid_region)
							cv = 1.d0-this%c_%f_(i,j,k)
						case default
					end select

					!pressure
					pv = this%p_%f_(i,j,k)
					pm(1,j) = pm(1,j) + cv*pv*V
					ppm(1,j) = ppm(1,j) + cv*pv*pv*V
					
					!velocity
					uxv = 0.5d0*(this%u_%ux_%f_(i,j,k)+this%u_%ux_%f_(i-1,j,k))
					uyv = 0.5d0*(this%u_%uy_%f_(i,j,k)+this%u_%uy_%f_(i,j-1,k))
					uzv = 0.5d0*(this%u_%uz_%f_(i,j,k)+this%u_%uz_%f_(i,j,k-1))
					um(1,j) = um(1,j) + cv*uxv*V
					um(2,j) = um(2,j) + cv*uyv*V
					um(3,j) = um(3,j) + cv*uzv*V
					uum(1,j) = uum(1,j) + cv*uxv*uxv*V
					uum(2,j) = uum(2,j) + cv*uxv*uyv*V
					uum(3,j) = uum(3,j) + cv*uxv*uzv*V
					uum(4,j) = uum(4,j) + cv*uyv*uyv*V
					uum(5,j) = uum(5,j) + cv*uyv*uzv*V
					uum(6,j) = uum(6,j) + cv*uzv*uzv*V
					
					!vorticity
					wxv = this%w_%ux_%f_(i,j,k)
					wyv = this%w_%uy_%f_(i,j,k)
					wzv = this%w_%uz_%f_(i,j,k)
					wm(1,j) = wm(1,j) + cv*wxv*V
					wm(2,j) = wm(2,j) + cv*wyv*V
					wm(3,j) = wm(3,j) + cv*wzv*V
					wwm(1,j) = wwm(1,j) + cv*wxv*wxv*V
					wwm(2,j) = wwm(2,j) + cv*wxv*wyv*V
					wwm(3,j) = wwm(3,j) + cv*wxv*wzv*V
					wwm(4,j) = wwm(4,j) + cv*wyv*wyv*V
					wwm(5,j) = wwm(5,j) + cv*wyv*wzv*V
					wwm(6,j) = wwm(6,j) + cv*wzv*wzv*V				
										
				end do
			end do
		end do
		!$OMP END PARALLEL DO  

		call reduceStat(pm,spm,this%Vs_,mpic,stat_scalar)
		call reduceStat(um,sum,this%Vs_,mpic,stat_vector)
		call reduceStat(wm,swm,this%Vs_,mpic,stat_vector)
		
		call reduceStat(ppm,sppm,this%Vs_,mpic,stat_scalar)
		call reduceStat(uum,suum,this%Vs_,mpic,stat_symtensor)
		call reduceStat(wwm,swwm,this%Vs_,mpic,stat_symtensor)
		
		!output time average shear stress liquid
		if (region==s_liquid_region) then
			!call shearStress_liquid(this,sum,t)
		end if
		
		!update time averages
		if (IS_MASTER) then
			call updateTimeAvFields(this,t,dt,spm,sppm,sum,suum,swm,swwm,region)
		end if

	end subroutine
!========================================================================================!

!========================================================================================!
	subroutine updateStatsGradU(this,t,dt,region)
		type(statistics), intent(in) :: this
		real(DP), intent(in) :: t, dt
		integer, intent(in) :: region
		type(grid), pointer :: mesh
		type(mpiControl), pointer :: mpic
		real(DP), allocatable, dimension(:,:) :: dum,dudum,sdum,sdudum
		real(DP) :: V,cv,duxdx,duxdy,duxdz,duydx,duydy,duydz,duzdx,duzdy,duzdz
		integer :: i,j,k,nx,ny,nz

		mesh => this%p_%ptrMesh_
		mpic => mesh%ptrMPIC_
		
		nx = this%p_%nx_
		ny = this%p_%ny_
		nz = this%p_%nz_
		
		!allocate tmp fields
		call allocateArray(dum,1,9,1,ny)
		call allocateArray(dudum,1,12,1,ny)
		
		!reset to zero
		dum = 0.d0
		dudum = 0.d0
		
		!mean along wall normal planes
		!$OMP PARALLEL DO DEFAULT(none) &
		!$OMP SHARED(this,mesh,region) &
		!$OMP SHARED(nx,ny,nz) &
		!$OMP PRIVATE(V,cv,duxdx,duxdy,duxdz,duydx,duydy,duydz,duzdx,duzdy,duzdz) &
		!$OMP PRIVATE(i,j,k) &
		!$OMP REDUCTION(+:dum,dudum)
		do k=1,nz
			do j=1,ny
				do i=1,nx

					V=mesh%V_(i,j,k)
					
					!volume fraction
					select case(region)
						case(s_whole_region)
							cv=1.d0
						case(s_gas_region)
							cv = this%c_%f_(i,j,k)
						case(s_liquid_region)
							cv = 1.d0-this%c_%f_(i,j,k)
						case default
					end select					
					
					!dux
					duxdx=(this%u_%ux_%f_(i,j,k)-this%u_%ux_%f_(i-1,j,k))/mesh%dxf_(i)	
					duxdy=(0.5d0*(this%u_%ux_%f_(i,j+1,k)+this%u_%ux_%f_(i-1,j+1,k))-  &
						   0.5d0*(this%u_%ux_%f_(i,j-1,k)+this%u_%ux_%f_(i-1,j-1,k)))/ &
						   (mesh%dyc_(j+1)+mesh%dyc_(j))
					duxdz=(0.5d0*(this%u_%ux_%f_(i,j,k+1)+this%u_%ux_%f_(i-1,j,k+1))-  &
						   0.5d0*(this%u_%ux_%f_(i,j,k-1)+this%u_%ux_%f_(i-1,j,k-1)))/ &
						   (mesh%dzc_(k+1)+mesh%dzc_(k))
					!duy
					duydx=(0.5d0*(this%u_%uy_%f_(i+1,j,k)+this%u_%uy_%f_(i+1,j-1,k))-  &
						   0.5d0*(this%u_%uy_%f_(i-1,j,k)+this%u_%uy_%f_(i-1,j-1,k)))/ &
						   (mesh%dxc_(i+1)+mesh%dxc_(i))
					duydy=(this%u_%uy_%f_(i,j,k)-this%u_%uy_%f_(i,j-1,k))/mesh%dyf_(j)
					duydz=(0.5d0*(this%u_%uy_%f_(i,j,k+1)+this%u_%uy_%f_(i,j-1,k+1))-  &
						   0.5d0*(this%u_%uy_%f_(i,j,k-1)+this%u_%uy_%f_(i,j-1,k-1)))/ &
						   (mesh%dzc_(k+1)+mesh%dzc_(k))
					
					!duz
					duzdx=(0.5d0*(this%u_%uz_%f_(i+1,j,k)+this%u_%uz_%f_(i+1,j,k-1))-  &
						   0.5d0*(this%u_%uz_%f_(i-1,j,k)+this%u_%uz_%f_(i-1,j,k-1)))/ &
						   (mesh%dxc_(i+1)+mesh%dxc_(i))
					duzdy=(0.5d0*(this%u_%uz_%f_(i,j+1,k)+this%u_%uz_%f_(i,j+1,k-1))-  &
						   0.5d0*(this%u_%uz_%f_(i,j-1,k)+this%u_%uz_%f_(i,j-1,k-1)))/ &
						   (mesh%dyc_(j+1)+mesh%dyc_(j))
					duzdz=(this%u_%uz_%f_(i,j,k)-this%u_%uz_%f_(i,j,k-1))/mesh%dzf_(k)	
					
					dum(1,j) = dum(1,j) + cv*duxdx*V
					dum(2,j) = dum(2,j) + cv*duxdy*V
					dum(3,j) = dum(3,j) + cv*duxdz*V
					dum(4,j) = dum(4,j) + cv*duydx*V
					dum(5,j) = dum(5,j) + cv*duydy*V
					dum(6,j) = dum(6,j) + cv*duydz*V
					dum(7,j) = dum(7,j) + cv*duzdx*V
					dum(8,j) = dum(8,j) + cv*duzdy*V
					dum(9,j) = dum(9,j) + cv*duzdz*V	
					
					dudum(1,j) = dudum(1,j) + cv*duxdx*duxdx*V
					dudum(2,j) = dudum(2,j) + cv*duxdy*duxdy*V
					dudum(3,j) = dudum(3,j) + cv*duxdz*duxdz*V
					dudum(4,j) = dudum(4,j) + cv*duydx*duydx*V
					dudum(5,j) = dudum(5,j) + cv*duydy*duydy*V
					dudum(6,j) = dudum(6,j) + cv*duydz*duydz*V
					dudum(7,j) = dudum(7,j) + cv*duzdx*duzdx*V
					dudum(8,j) = dudum(8,j) + cv*duzdy*duzdy*V
					dudum(9,j) = dudum(9,j) + cv*duzdz*duzdz*V	
					!mixed products
					dudum(10,j) = dudum(10,j) + cv*duxdy*duydx*V
					dudum(11,j) = dudum(11,j) + cv*duxdz*duzdx*V
					dudum(12,j) = dudum(12,j) + cv*duydz*duzdy*V
														
				end do
			end do
		end do
		!$OMP END PARALLEL DO

		call reduceStat(dum,sdum,this%Vs_,mpic,stat_tensor)
		call reduceStat(dudum,sdudum,this%Vs_,mpic,stat_tensor_sqr)
	
		!update time averages
		if (IS_MASTER) then
			call updateTimeAvGradU(this,t,dt,sdum,sdudum,region)
		end if

	end subroutine
!========================================================================================!

!========================================================================================!
	subroutine updateTimeAvVF(this,t,dt,scm)
		type(statistics), intent(in) :: this
		real(DP), intent(in) :: t, dt
		real(DP), allocatable, dimension(:,:), intent(in) :: scm
		integer :: j,nyg
		
		nyg = this%p_%ptrMesh_%nyg_

		do j=1,nyg							
			cm_(1,j) = ( cm_(1,j)*(t-this%Ts_)+scm(1,j)*dt ) / (t-this%Ts_+dt)
		end do

	end subroutine
!========================================================================================!

!========================================================================================!
	subroutine updateTimeAvFields(this,t,dt,spm,sppm,sum,suum,swm,swwm,region)
		type(statistics), intent(in) :: this
		real(DP), intent(in) :: t, dt
		real(DP), allocatable, dimension(:,:), intent(in) :: spm,sppm,sum,suum,swm,swwm
		integer, intent(in) :: region
		integer :: j,nyg

		nyg = this%p_%ptrMesh_%nyg_

		select case(region)
			case(s_whole_region)
				do j=1,nyg
					pm_(1,j) = ( pm_(1,j)*(t-this%Ts_)+spm(1,j)*dt ) / (t-this%Ts_+dt)
					um_(:,j) = ( um_(:,j)*(t-this%Ts_)+sum(:,j)*dt ) / (t-this%Ts_+dt)
					wm_(:,j) = ( wm_(:,j)*(t-this%Ts_)+swm(:,j)*dt ) / (t-this%Ts_+dt)
					ppm_(1,j) = ( ppm_(1,j)*(t-this%Ts_)+sppm(1,j)*dt ) / (t-this%Ts_+dt)				
					uum_(:,j) = ( uum_(:,j)*(t-this%Ts_)+suum(:,j)*dt ) / (t-this%Ts_+dt)
					wwm_(:,j) = ( wwm_(:,j)*(t-this%Ts_)+swwm(:,j)*dt ) / (t-this%Ts_+dt)
				end do				
			case(s_gas_region)
				do j=1,nyg		
					pmg_(1,j) = ( pmg_(1,j)*(t-this%Ts_)+spm(1,j)*dt ) / (t-this%Ts_+dt)
					umg_(:,j) = ( umg_(:,j)*(t-this%Ts_)+sum(:,j)*dt ) / (t-this%Ts_+dt)
					wmg_(:,j) = ( wmg_(:,j)*(t-this%Ts_)+swm(:,j)*dt ) / (t-this%Ts_+dt)
					ppmg_(1,j) = ( ppmg_(1,j)*(t-this%Ts_)+sppm(1,j)*dt ) / (t-this%Ts_+dt)				
					uumg_(:,j) = ( uumg_(:,j)*(t-this%Ts_)+suum(:,j)*dt ) / (t-this%Ts_+dt)
					wwmg_(:,j) = ( wwmg_(:,j)*(t-this%Ts_)+swwm(:,j)*dt ) / (t-this%Ts_+dt)	
				end do
			case(s_liquid_region)
				do j=1,nyg			
					pml_(1,j) = ( pml_(1,j)*(t-this%Ts_)+spm(1,j)*dt ) / (t-this%Ts_+dt)
					uml_(:,j) = ( uml_(:,j)*(t-this%Ts_)+sum(:,j)*dt ) / (t-this%Ts_+dt)
					wml_(:,j) = ( wml_(:,j)*(t-this%Ts_)+swm(:,j)*dt ) / (t-this%Ts_+dt)
					ppml_(1,j) = ( ppml_(1,j)*(t-this%Ts_)+sppm(1,j)*dt ) / (t-this%Ts_+dt)				
					uuml_(:,j) = ( uuml_(:,j)*(t-this%Ts_)+suum(:,j)*dt ) / (t-this%Ts_+dt)
					wwml_(:,j) = ( wwml_(:,j)*(t-this%Ts_)+swwm(:,j)*dt ) / (t-this%Ts_+dt)	
				end do
			case default
		end select
		
	end subroutine
!========================================================================================!

!========================================================================================!
	subroutine updateTimeAvGradU(this,t,dt,sdum,sdudum,region)
		type(statistics), intent(in) :: this
		real(DP), intent(in) :: t, dt
		real(DP), allocatable, dimension(:,:), intent(in) :: sdum,sdudum
		integer, intent(in) :: region
		integer :: j,nyg
		
		nyg = this%p_%ptrMesh_%nyg_

		select case(region)
			case(s_whole_region)
				do j=1,nyg
					dum_(:,j) = ( dum_(:,j)*(t-this%Ts_)+sdum(:,j)*dt ) / (t-this%Ts_+dt)
					dudum_(:,j) = ( dudum_(:,j)*(t-this%Ts_)+sdudum(:,j)*dt ) / (t-this%Ts_+dt)
				end do				
			case(s_gas_region)
				do j=1,nyg		
					dumg_(:,j) = ( dumg_(:,j)*(t-this%Ts_)+sdum(:,j)*dt ) / (t-this%Ts_+dt)
					dudumg_(:,j) = ( dudumg_(:,j)*(t-this%Ts_)+sdudum(:,j)*dt ) / (t-this%Ts_+dt)
				end do
			case(s_liquid_region)
				do j=1,nyg			
					duml_(:,j) = ( duml_(:,j)*(t-this%Ts_)+sdum(:,j)*dt ) / (t-this%Ts_+dt)
					duduml_(:,j) = ( duduml_(:,j)*(t-this%Ts_)+sdudum(:,j)*dt ) / (t-this%Ts_+dt)
				end do
			case default
		end select

	end subroutine
!========================================================================================!

!========================================================================================!
	subroutine reduceStat(ql,qg,V,mpic,stat_type)
		real(DP), allocatable, dimension(:,:), intent(in) :: ql
		real(DP), allocatable, dimension(:,:), intent(out) :: qg
		real(DP), allocatable, dimension(:), intent(in) :: V
		integer, intent(in) :: stat_type
		type(mpiControl) :: mpic
		real(DP), allocatable, dimension(:,:) :: tmp
		integer :: j,l,n,ny,nyg,j0g,j1g,dim,nproc,ierror
		
		
		ny = ubound(ql,2)
		nyg = ubound(V,1)
		nproc = mpic%nProcs_
		
		select case(stat_type)
			case(stat_scalar)
				dim = 1
			case(stat_vector)
				dim = 3
			case(stat_symtensor)
				dim = 6
			case(stat_tensor)
				dim = 9
			case(stat_tensor_sqr)
				dim = 12
			case default
		end select
		
		if (IS_MASTER) then
			call allocateArray(qg,1,dim,1,nyg)
		end if
		call allocateArray(tmp,1,dim,1,ny*nproc)
		
		call MPI_Gather(ql,dim*ny,MPI_DOUBLE_PRECISION,tmp,dim*ny,MPI_DOUBLE_PRECISION,&
						0,mpic%cartComm_,ierror)
		
		if (IS_MASTER) then
		
			qg = 0.d0
		
			do n=0,nproc-1
				j0g = mpic%gCoords_(2,n)*ny+1
				j1g = j0g + ny - 1
				
				l=ny*n+1
				do j=j0g,j1g	
					qg(:,j)=qg(:,j)+tmp(:,l)
					l=l+1
				end do
				
			end do
			
			do j=1,nyg
				qg(:,j) = qg(:,j)/V(j)
			end do
			
		end if
		
	end subroutine
!========================================================================================!

!========================================================================================!
	subroutine computeSliceVolume(this)
		type(statistics), intent(inout) :: this
		type(mpiControl), pointer :: mpic
		real(DP) :: tmp
		integer :: i,j,k
		integer :: nx, ny, nz, nyg, j0g, j1g
		integer :: ierror
		
		!global indexes wall normal
		j0g = this%p_%ptrMesh_%j0g_
		j1g = this%p_%ptrMesh_%j1g_
		
		nx = this%p_%nx_
		ny = this%p_%ny_
		nz = this%p_%nz_
		
		nyg =  this%p_%ptrMesh_%nyg_
		
		mpic => this%p_%ptrMesh_%ptrMPIC_
		
		!mean along wall normal planes
		do k=1,nz
			do j=j0g,j1g
				do i=1,nx		
				
					!mean Volume
					this%Vs_(j) = this%Vs_(j) + this%p_%ptrMesh_%V_(i,j-j0g+1,k)
					
				end do
			end do
		end do
		
		!reduce to the master node		
		do j=1,nyg
			call Mpi_reduce(this%Vs_(j), tmp, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, mpic%cartComm_, ierror)
			if (IS_MASTER) then
				this%Vs_(j) = tmp
			end if
		end do
		
		
	end subroutine
!========================================================================================!

!========================================================================================!
	subroutine writeStats(this,nf)
		type(statistics), intent(in) :: this
		integer, intent(in) :: nf
		
		call writeStats_VF(this,nf,cm_)
		
		call writeStats_region(this,nf,pm_,ppm_,um_,uum_,wm_,wwm_,dum_,dudum_,s_whole_region)
		call writeStats_region(this,nf,pmg_,ppmg_,umg_,uumg_,wmg_,wwmg_,dumg_,dudumg_,s_gas_region)
		call writeStats_region(this,nf,pml_,ppml_,uml_,uuml_,wml_,wwml_,duml_,duduml_,s_liquid_region)
		
	end subroutine
!========================================================================================!

!========================================================================================!
	subroutine writeStats_VF(this,nf,cm)
		type(statistics), intent(in) :: this
		integer, intent(in) :: nf
		real(DP), allocatable, dimension(:,:) ::cm
		integer :: j, ny
		character(len=20) :: dirName
		
		ny = this%p_%ptrMesh_%nyg_
		
		if (IS_MASTER) then
        
        	write(dirName,s_intFormat) nf
        	
			!volume fraction
			open(UNIT=s_IOunitNumber,FILE=adjustl(trim(dirName)//'/stats_c'),&
				 STATUS='REPLACE',ACTION='WRITE')
			do j=1,ny
				write(s_IOunitNumber,'(1'//s_doubleFormat(2:10)//')') cm(1,j)
			end do
			close(s_IOunitNumber)		
		
		end if
		
		
	end subroutine
!========================================================================================!

!========================================================================================!
	subroutine writeStats_region(this,nf,pm,ppm,um,uum,wm,wwm,dum,dudum,region)
		type(statistics), intent(in) :: this
		integer, intent(in) :: nf, region
		real(DP), allocatable, dimension(:,:) :: pm,ppm,um,uum,wm,wwm,dum,dudum
		integer :: j, ny
		character(len=20) :: dirName
		character(len=2) :: char_reg
		
		ny = this%p_%ptrMesh_%nyg_
		
		if (IS_MASTER) then
        
        	write(dirName,s_intFormat) nf
        	
        	select case(region)
        		case(s_whole_region)
        			char_reg = '_W'
        		case(s_gas_region)
        			char_reg = '_G'
        		case(s_liquid_region)
        			char_reg = '_L'
        		case default
        	end select
		
			!pressure
			open(UNIT=s_IOunitNumber,FILE=adjustl(trim(dirName)//'/stats_p'//char_reg),&
				 STATUS='REPLACE',ACTION='WRITE')
			do j=1,ny
				write(s_IOunitNumber,'(2'//s_doubleFormat(2:10)//')') pm(1,j),ppm(1,j)
			end do
			close(s_IOunitNumber)
			
			!velocity
			open(UNIT=s_IOunitNumber,FILE=adjustl(trim(dirName)//'/stats_u'//char_reg),&
				 STATUS='REPLACE',ACTION='WRITE')
			do j=1,ny	
				write(s_IOunitNumber,'(9'//s_doubleFormat(2:10)//')') &
					  um(1,j),um(2,j),um(3,j),	  &
					  uum(1,j),uum(2,j),uum(3,j),  &
					  uum(4,j),uum(5,j),uum(6,j)
			end do
			close(s_IOunitNumber)
			
			!vorticity
			open(UNIT=s_IOunitNumber,FILE=adjustl(trim(dirName)//'/stats_w'//char_reg),&
				 STATUS='REPLACE',ACTION='WRITE')
			do j=1,ny	
				write(s_IOunitNumber,'(9'//s_doubleFormat(2:10)//')') &
					  wm(1,j),wm(2,j),wm(3,j),	  &
					  wwm(1,j),wwm(2,j),wwm(3,j),  &
					  wwm(4,j),wwm(5,j),wwm(6,j)
			end do
			close(s_IOunitNumber)
			
			!gradU
			open(UNIT=s_IOunitNumber,FILE=adjustl(trim(dirName)//'/stats_gradU'//char_reg),&
				 STATUS='REPLACE',ACTION='WRITE')
			do j=1,ny	
				write(s_IOunitNumber,'(21'//s_doubleFormat(2:10)//')') &
					  dum(1,j),dum(2,j),dum(3,j),  		&
					  dum(4,j),dum(5,j),dum(6,j),  		&
					  dum(7,j),dum(8,j),dum(9,j),  		&
					  dudum(1,j),dudum(2,j),dudum(3,j),  &
					  dudum(4,j),dudum(5,j),dudum(6,j),  &
					  dudum(7,j),dudum(8,j),dudum(9,j),	&
					  dudum(10,j),dudum(11,j),dudum(12,j)
			end do
			close(s_IOunitNumber)			
		
		end if
		
		
	end subroutine
!========================================================================================!

!========================================================================================!
	subroutine shearStress_liquid(this,t)
		type(statistics), intent(in) :: this
		real(DP), intent(in) :: t
		real(DP) :: y1,y2,u1,u2,dudy,tw
		integer :: ny
		
		if (IS_MASTER) then
			
			ny = this%gmesh_%ny_
				
			!bottom wall
			y1=this%gmesh_%yc_(1)
			y2=this%gmesh_%yc_(2)
			u1=um_(1,1)
			u2=um_(1,2)
			dudy = -(y2*y2*u1-y1*y1*u2)/(y1*y2*(y1-y2))

			tw = this%nu_*dudy

			!top wall
			y1=this%gmesh_%yf_(ny)-this%gmesh_%yc_(ny)
			y2=this%gmesh_%yf_(ny)-this%gmesh_%yc_(ny-1)
			u1=um_(1,ny)
			u2=um_(1,ny-1)
			dudy = -(y2*y2*u1-y1*y1*u2)/(y1*y2*(y1-y2))	
			
			tw = tw + this%nu_*dudy	
			
			write(*,'(A,'//s_doubleFormat(2:10)//',A,'//s_doubleFormat(2:10)//')') &
			        '	Total tau_w: ', tw, ' t= ', t
			
		end if
		

	end subroutine
!========================================================================================!

!========================================================================================!
	subroutine shearStress(this,t,dt)
		type(statistics), intent(in) :: this
		real(DP), intent(in) :: t,dt
		type(vfield), pointer :: u => NULL()
		type(field), pointer :: mu => NULL()
		type(mpiControl), pointer :: mpic => NULL()
		integer :: i,j,k,ip,im,jp,jm,kp,km,q,nx,ny,nz,ierror,pmin,pmax,pcoord
		real(DP), dimension(2) :: sw
		integer, dimension(2) :: plate,qloop
		real(DP) :: mur,mul,duxdy,duydx,A,dx,dz,tau,tau_g,Lx,Lz
	

		u  => this%u_
		mu => this%mu_
		mpic => this%u_%ptrMesh_%ptrMPIC_

		nx = u%ptrMesh_%nx_
		ny = u%ptrMesh_%ny_
		nz = u%ptrMesh_%nz_
		
		Lx = u%ptrMesh_%Lxg_
		Lz = u%ptrMesh_%Lzg_

		tau=0.d0
		plate(1)=1
		plate(2)=ny
		
		!check procs adjacent to wall only
		pmin = 0
		pmax = mpic%nProcsAxis_(2)-1
		pcoord = mpic%procCoord_(2)
		
		if ((pcoord==pmin).AND.(pcoord==pmax)) then
			qloop(1)=1
			qloop(2)=2
		elseif (pcoord==pmin) then
			qloop(1)=1
			qloop(2)=1
		elseif (pcoord==pmax) then
			qloop(1)=2
			qloop(2)=2	
		else
			qloop(1)=1
			qloop(2)=0		
		end if		

		
		do q=qloop(1),qloop(2)
			
			j=plate(q)
			jm=j-1
			jp=j+1
		
			if (q==1) then
				sw(1)=0.d0
				sw(2)=1.d0
			else
				sw(1)=1.d0
				sw(2)=0.d0
			end if
		
			do k = 0,nz-1
				kp = k + 1
				km = k - 1
				do i = 0,nx-1
			
					ip = i + 1
					im = i - 1
				
					dx=u%ptrMesh_%dxc_(ip)
					dz=u%ptrMesh_%dzf_(k)
					A=dx*dz
								
					mur = sw(1)*0.25d0*(mu%f_(ip,j,k)+mu%f_(i,j,k)+mu%f_(ip,jp,k)+mu%f_(i,jp,k))
					mul = sw(2)*0.25d0*(mu%f_(ip,j,k)+mu%f_(i,j,k)+mu%f_(ip,jm,k)+mu%f_(i,jm,k)) 	  
					duxdy =   mur*(u%ux_%f_(i,jp,k)-u%ux_%f_(i,j,k))/(u%ptrMesh_%dyc_(jp))	&
				     	    - mul*(u%ux_%f_(i,j,k)-u%ux_%f_(i,jm,k))/(u%ptrMesh_%dyc_(j))

					tau = tau + A*duxdy
				
				end do
			end do
			
		end do
		
        call Mpi_reduce(tau, tau_g, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, mpic%cartComm_, ierror)

		if (IS_MASTER) then
		
			tau_g=tau_g/(Lx*Lz)
			write(*,'(A,'//s_doubleFormat(2:10)//')') '	(tau_w)_A:  ', tau_g		
		
			if (this%isTimeAverage_) then
				tauw_tav = ( tauw_tav*(t-this%Ts_)+tau_g*dt ) / (t-this%Ts_+dt)
				write(*,'(A,'//s_doubleFormat(2:10)//')') '	(tau_w)_A,t:', tauw_tav
			end if

		end if
	
	end subroutine
!========================================================================================!

!========================================================================================!
    subroutine checkTimeAverage(this,t)
    	type(statistics), intent(inout) :: this
    	real(DP), intent(in) :: t

    	
    	if (t >= this%Ts_) then    	
    		this%isTimeAverage_ = .TRUE.
    		this%Ts_=t
    		this%hasTimeAvStarted_=.TRUE.
		else		
			this%isTimeAverage_ = .FALSE.			
		end if
    	
    end subroutine
!========================================================================================!

!========================================================================================!
    subroutine info(cpuTime)
    	real(DP) :: cpuTime

		if (IS_MASTER) then
			write(*,'(A,'//s_outputFormat(2:9)//')') '	Statistics: CPU time = ', cpuTime
		end if
		
    	
    end subroutine
!========================================================================================!

	
end module statisticsMod


	




