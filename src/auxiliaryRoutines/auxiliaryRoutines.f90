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

module auxiliaryRoutinesMod

	use vfieldMod

	implicit none

	public :: computeContinuityError
	public :: computeCFLmax
	public :: compute_dt_CFL
	public :: compute_Umag_max
	public :: updateShearFlow
	public :: shapeError
	public :: setPressGrad
	public :: scaleVelocity
	public :: computeVorticity	
	public :: setFlowRate
	public :: computeFlowRate
	
contains


!========================================================================================!
    subroutine computeContinuityError(u,dt)
		type(vfield), intent(in) :: u
		real(DP), intent(in) :: dt
		type(grid), pointer :: mesh
		type(mpiControl), pointer :: mpic
		real(DP) :: dx, dy, dz
		real(DP) :: cError, cErrorg
		integer :: i, j, k, nx, ny, nz
		integer :: ierror
        
        mesh => u%ptrMesh_
        mpic => mesh%ptrMPIC_
        
        cError = -1.d0
        
        nx = mesh%nx_
        ny = mesh%ny_
        nz = mesh%nz_       
        
		!$OMP PARALLEL DO DEFAULT(none) &
		!$OMP SHARED(u,mesh) &
		!$OMP SHARED(nx,ny,nz) &
		!$OMP PRIVATE(dx,dy,dz) &
		!$OMP PRIVATE(i,j,k) &
		!$OMP REDUCTION(max:cError)
        do k = 1,nz
        	do j = 1,ny
        		do i = 1,nx
        		
        			dx = (u%ux_%f_(i,j,k)-u%ux_%f_(i-1,j,k))/mesh%dxf_(i)
        			dy = (u%uy_%f_(i,j,k)-u%uy_%f_(i,j-1,k))/mesh%dyf_(j)
        			dz = (u%uz_%f_(i,j,k)-u%uz_%f_(i,j,k-1))/mesh%dzf_(k)
        			
        			cError = max(abs(dx+dy+dz),cError)
        			
        		end do
        	end do
        end do
        !$OMP END PARALLEL DO 
        
        cError = cError*dt            
        
        call Mpi_reduce(cError, cErrorg, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, mpic%cartComm_, ierror)
        
        if (IS_MASTER) then
        	write(*,'(A,'//s_outputFormat(2:9)//')') '	continuity error: ', cErrorg
        end if

        
    end subroutine
!========================================================================================!

!========================================================================================!
    function computeCFLmax(u,dt) result(cflg)
		type(vfield), intent(in) :: u
		real(DP), intent(in) :: dt
		type(grid), pointer :: mesh
		type(mpiControl), pointer :: mpic
		real(DP) :: cx, cy, cz
		real(DP) :: cfl, cflg
		integer :: i, j, k, nx, ny, nz
		integer :: ierror
        
        mesh => u%ptrMesh_
        mpic => mesh%ptrMPIC_
        
        cfl = -1.d0
        
        nx = mesh%nx_
        ny = mesh%ny_
        nz = mesh%nz_       
        
		!$OMP PARALLEL DO DEFAULT(none) &
		!$OMP SHARED(u,mesh) &
		!$OMP SHARED(nx,ny,nz) &
		!$OMP PRIVATE(cx,cy,cz) &
		!$OMP PRIVATE(i,j,k) &
		!$OMP REDUCTION(max:cfl)
        do k = 1,nz
        	do j = 1,ny
        		do i = 1,nx
        		
        			cx = abs(0.5d0*(u%ux_%f_(i,j,k)+u%ux_%f_(i-1,j,k)))/mesh%dxf_(i)
        			cy = abs(0.5d0*(u%uy_%f_(i,j,k)+u%uy_%f_(i,j-1,k)))/mesh%dyf_(j)
        			cz = abs(0.5d0*(u%uz_%f_(i,j,k)+u%uz_%f_(i,j,k-1)))/mesh%dzf_(k)
        			
        			cfl = max(cx+cy+cz,cfl)
        			
        			
        		end do
        	end do
        end do   
        !$OMP END PARALLEL DO           
        
        cfl = cfl*dt
        
        call Mpi_Allreduce(cfl, cflg, 1, MPI_DOUBLE_PRECISION, MPI_MAX, mpic%cartComm_, ierror)

        
    end function
!========================================================================================!

!========================================================================================!
    function compute_dt_CFL(u,cfl) result(dt_g)
		type(vfield), intent(in) :: u
		real(DP), intent(in) :: cfl
		type(grid), pointer :: mesh
		type(mpiControl), pointer :: mpic
		real(DP) :: cx,cy,cz,dt,dt_g,small
		integer :: i,j,k,nx,ny,nz
		integer :: ierror
        
        mesh => u%ptrMesh_
        mpic => mesh%ptrMPIC_
        
        dt = huge(0.d0)
        
        small = tiny(0.d0)
        
        nx = mesh%nx_
        ny = mesh%ny_
        nz = mesh%nz_       
        
		!$OMP PARALLEL DO DEFAULT(none) &
		!$OMP SHARED(u,mesh,cfl,small) &
		!$OMP SHARED(nx,ny,nz) &
		!$OMP PRIVATE(cx,cy,cz) &
		!$OMP PRIVATE(i,j,k) &
		!$OMP REDUCTION(min:dt)
        do k = 1,nz
        	do j = 1,ny
        		do i = 1,nx
        		
        			cx = abs(0.5d0*(u%ux_%f_(i,j,k)+u%ux_%f_(i-1,j,k)))/mesh%dxf_(i)
        			cy = abs(0.5d0*(u%uy_%f_(i,j,k)+u%uy_%f_(i,j-1,k)))/mesh%dyf_(j)
        			cz = abs(0.5d0*(u%uz_%f_(i,j,k)+u%uz_%f_(i,j,k-1)))/mesh%dzf_(k)
        			
        			dt = min(cfl/(cx+cy+cz+small),dt)
        			
        		end do
        	end do
        end do   
        !$OMP END PARALLEL DO
        
        call Mpi_Allreduce(dt, dt_g, 1, MPI_DOUBLE_PRECISION, MPI_MIN, mpic%cartComm_, ierror)

        
    end function
!========================================================================================!

!========================================================================================!
    subroutine compute_Umag_max(u)
		type(vfield), intent(in) :: u
		type(grid), pointer :: mesh
		type(mpiControl), pointer :: mpic
		real(DP) :: ux, uy, uz, um, umax, umax_g
		integer :: i, j, k, nx, ny, nz
		integer :: ierror
        
        mesh => u%ptrMesh_
        mpic => mesh%ptrMPIC_
        
        um = -1.d0
        
        nx = mesh%nx_
        ny = mesh%ny_
        nz = mesh%nz_       
        
		!$OMP PARALLEL DO DEFAULT(none) &
		!$OMP SHARED(u,mesh) &
		!$OMP SHARED(nx,ny,nz) &
		!$OMP PRIVATE(ux,uy,uz,um) &
		!$OMP PRIVATE(i,j,k) &
		!$OMP REDUCTION(max:umax)
        do k = 1,nz
        	do j = 1,ny
        		do i = 1,nx  		
        			ux = 0.5d0*(u%ux_%f_(i,j,k)+u%ux_%f_(i-1,j,k))
        			uy = 0.5d0*(u%uy_%f_(i,j,k)+u%uy_%f_(i,j-1,k))
        			uz = 0.5d0*(u%uz_%f_(i,j,k)+u%uz_%f_(i,j,k-1))
        			um=sqrt(ux*ux+uy*uy+uz*uz)
        			umax = max(umax,um)
        		end do
        	end do
        end do   
        !$OMP END PARALLEL DO
        
        call Mpi_reduce(umax,umax_g,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,mpic%cartComm_,ierror)
        
        if (IS_MASTER) then
        	write(*,'(A,'//s_outputFormat(2:9)//')') '	Umax: ', umax_g
        end if

        
    end subroutine
!========================================================================================!

!========================================================================================!
    subroutine updateShearFlow(u,u0,t,dt)
		type(vfield), intent(inout) :: u
		type(vfield), intent(in) :: u0
		real(DP), intent(in) :: t, dt
		integer :: i, j, k
		real(DP) :: Ts, th
		real(DP), parameter :: pi = 4.d0*DATAN(1.d0)
		
		Ts = 3.d0
		th = t - 0.5d0*dt
		
        
		!ux
		do k=u%ux_%ks_,u%ux_%ke_
			do j=u%ux_%js_,u%ux_%je_
				do i=u%ux_%is_,u%ux_%ie_
					
					u%ux_%f_(i,j,k) = u0%ux_%f_(i,j,k)*cos(pi*th/Ts)
				
				end do
			end do
		end do
		
		!uy
		do k=u%uy_%ks_,u%uy_%ke_
			do j=u%uy_%js_,u%uy_%je_
				do i=u%uz_%is_,u%uz_%ie_
					
					u%uy_%f_(i,j,k) = u0%uy_%f_(i,j,k)*cos(pi*th/Ts)
				
				end do
			end do
		end do
		
		!uz
		do k=u%uz_%ks_,u%uz_%ke_
			do j=u%uz_%js_,u%uz_%je_
				do i=u%uz_%is_,u%uz_%ie_
					
					u%uz_%f_(i,j,k) = u0%uz_%f_(i,j,k)*cos(pi*th/Ts)
				
				end do
			end do
		end do

        
    end subroutine
!========================================================================================!

!========================================================================================!
    subroutine shapeError(c,c0,t)
		type(field), intent(in) :: c, c0
		real(DP), intent(in) :: t
		type(grid), pointer :: mesh
		type(mpiControl), pointer :: mpic
		integer :: n, nx, ny, nz
		real(DP) :: lerr, gerr
		integer :: ierror
		logical :: exist
		
		lerr = 0.d0
		
		mesh => c%ptrMesh_
		mpic => mesh%ptrMPIC_
		
		nx = mesh%nx_
		ny = mesh%ny_
		nz = mesh%nz_
		
		n = mesh%nxg_*mesh%nyg_*mesh%nzg_
		
		lerr = lerr + sum(abs(c%f_(1:nx,1:ny,1:nz)-c0%f_(1:nx,1:ny,1:nz)))/n
		
		call Mpi_reduce(lerr, gerr, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, mpic%cartComm_, ierror)
		
		!output to a text file
		if (IS_MASTER) then
		
			inquire(file="shError", exist=exist)
		
			if (exist) then
				open(UNIT=s_IOunitNumber,FILE='shError',STATUS='OLD',&
					 POSITION="append",ACTION='WRITE')
					write(s_IOunitNumber,'(2'//s_doubleFormat(2:10)//')') gerr,t
				close(s_IOunitNumber)	
			else
				open(UNIT=s_IOunitNumber,FILE='shError',STATUS='NEW',ACTION='WRITE')
					write(s_IOunitNumber,'(2'//s_doubleFormat(2:10)//')') gerr,t
				close(s_IOunitNumber)				
			end if	

		end if
		
        
    end subroutine
!========================================================================================!

!========================================================================================!
     subroutine setPressGrad(f_ctrl,Ret,rho,rhol,mul,g,f)
    	integer, intent(in) :: f_ctrl
    	type(field), intent(in) :: rho
    	real(DP), intent(in) :: rhol,mul,g,Ret
    	real(DP), intent(inout) :: f
    	type(grid), pointer :: mesh
		type(mpiControl), pointer :: mpic
		real(DP) :: x,rho_sum,rho_sum_g,rhoAv,nul,Re,tauw
		integer :: nx,ny,nz,i,j,k,ierror
		
		
        mesh => rho%ptrMesh_
        mpic => mesh%ptrMPIC_
		
		nx = mesh%nx_
		ny = mesh%ny_
		nz = mesh%nz_
		
		rho_sum=0.d0
		!$OMP PARALLEL DO DEFAULT(none) &
		!$OMP SHARED(rho,mesh) &
		!$OMP SHARED(nx,ny,nz) &
		!$OMP PRIVATE(i,j,k,x) &
		!$OMP REDUCTION(+:rho_sum)
		do k=1,nz
			do j=1,ny
				do i=0,nx-1
    				x = 0.5d0*(rho%f_(i+1,j,k)+rho%f_(i,j,k))
    				rho_sum = rho_sum + x*mesh%Vsx_(i,j,k)
				end do
			end do
		end do
		!$OMP END PARALLEL DO
		
    	call Mpi_Allreduce(rho_sum, rho_sum_g, 1, MPI_DOUBLE_PRECISION, MPI_SUM,  &
    				   	   mpic%cartComm_, ierror)		
		
		rhoAv = rho_sum_g/s_Vg
			
		if (f_ctrl==1) then
		
			nul = mul/rhol
			Re = 1.d0/nul
			tauw = (Ret/Re)*(Ret/Re)
		
			f = tauw - rhoAv*g
		
		elseif ((f_ctrl==2).OR.(f_ctrl==3)) then
						   
    		f = -rhoAv*g 
    		
		else
		
			f = 0.d0

		end if
        
    end subroutine
!========================================================================================!

!========================================================================================!
    subroutine scaleVelocity(u)
		type(vfield), intent(inout) :: u
		real(DP) :: Ubar
		
		Ubar = 1.85d1
		
		u%ux_%f_ = u%ux_%f_/Ubar
		u%uy_%f_ = u%uy_%f_/Ubar
		u%uz_%f_ = u%uz_%f_/Ubar

        
    end subroutine
!========================================================================================!

!========================================================================================!
    subroutine computeVorticity(u,w)
		type(vfield), intent(in) :: u
		type(vfield), intent(inout) :: w
		type(grid), pointer :: mesh
		real(DP) :: duzdy,duydz,duzdx,duxdz,duydx,duxdy
		integer :: i,j,k,nx,ny,nz
		
		mesh => u%ptrMesh_
		nx = mesh%nx_
		ny = mesh%ny_
		nz = mesh%nz_

		!$OMP PARALLEL DO DEFAULT(none) &
		!$OMP SHARED(u,w,mesh) &
		!$OMP SHARED(nx,ny,nz) &
		!$OMP PRIVATE(i,j,k) &
		!$OMP PRIVATE(duzdy,duydz,duzdx,duxdz,duydx,duxdy) 	
		do k=1,nz
			do j=1,ny
				do i=1,nx
					
					duzdy=(0.5d0*(u%uz_%f_(i,j+1,k)+u%uz_%f_(i,j+1,k-1))-  &
						   0.5d0*(u%uz_%f_(i,j-1,k)+u%uz_%f_(i,j-1,k-1)))/ &
						   (mesh%dyc_(j+1)+mesh%dyc_(j))
					duydz=(0.5d0*(u%uy_%f_(i,j,k+1)+u%uy_%f_(i,j-1,k+1))-  &
						   0.5d0*(u%uy_%f_(i,j,k-1)+u%uy_%f_(i,j-1,k-1)))/ &
						   (mesh%dzc_(k+1)+mesh%dzc_(k))
					duzdx=(0.5d0*(u%uz_%f_(i+1,j,k)+u%uz_%f_(i+1,j,k-1))-  &
						   0.5d0*(u%uz_%f_(i-1,j,k)+u%uz_%f_(i-1,j,k-1)))/ &
						   (mesh%dxc_(i+1)+mesh%dxc_(i))
					duxdz=(0.5d0*(u%ux_%f_(i,j,k+1)+u%ux_%f_(i-1,j,k+1))-  &
						   0.5d0*(u%ux_%f_(i,j,k-1)+u%ux_%f_(i-1,j,k-1)))/ &
						   (mesh%dzc_(k+1)+mesh%dzc_(k))
					duydx=(0.5d0*(u%uy_%f_(i+1,j,k)+u%uy_%f_(i+1,j-1,k))-  &
						   0.5d0*(u%uy_%f_(i-1,j,k)+u%uy_%f_(i-1,j-1,k)))/ &
						   (mesh%dxc_(i+1)+mesh%dxc_(i))
					duxdy=(0.5d0*(u%ux_%f_(i,j+1,k)+u%ux_%f_(i-1,j+1,k))-  &
						   0.5d0*(u%ux_%f_(i,j-1,k)+u%ux_%f_(i-1,j-1,k)))/ &
						   (mesh%dyc_(j+1)+mesh%dyc_(j))
						   
					w%ux_%f_(i,j,k)=duzdy-duydz
					w%uy_%f_(i,j,k)=duxdz-duzdx
					w%uz_%f_(i,j,k)=duydx-duxdy
					
					
				end do
			end do
		end do
		!$OMP END PARALLEL DO

        
    end subroutine
!========================================================================================!

!========================================================================================!
    subroutine setFlowRate(ctrl,u,rho,Q,dt,alpha)
    	integer, intent(in) :: ctrl
		type(vfield), intent(inout) :: u
		type(field), intent(in) :: rho
		real(DP), intent(in) :: Q,dt,alpha
		type(grid), pointer :: mesh
		type(mpiControl), pointer :: mpic
		real(DP) :: x,irho,uf,Qs,irho_g,r
		integer :: nx,ny,nz,i,j,k,ierror

		if ((ctrl<2).OR.(ctrl>2)) then
			return
		end if

        mesh => u%ptrMesh_
        mpic => mesh%ptrMPIC_
        
        nx=mesh%nx_
        ny=mesh%ny_
        nz=mesh%nz_
    
    	call computeFlowRate(u,Qs)

    	uf = (Q-Qs)/s_Vg
    	
		!$OMP PARALLEL DO DEFAULT(none) &
		!$OMP PRIVATE(i,j,k) &
		!$OMP SHARED(u,uf,nx,ny,nz)
    	do k=1,nz
    		do j=1,ny
    			do i=0,nx		
    				
    				u%ux_%f_(i,j,k) = u%ux_%f_(i,j,k) + uf
    				
    			end do
    		end do
    	end do
    	!$OMP END PARALLEL DO
    	call updateBoundariesV(u)

    	if (IS_MASTER) then
    		!write(*,'(A,'//s_outputFormat(2:9)//')') '	Press grad set to: ', uf/(alpha*dt)
    	end if
    
    end subroutine
!========================================================================================!

!========================================================================================!
    subroutine computeFlowRate(u,Q_g)
		type(vfield), intent(in) :: u
		real(DP), intent(out) :: Q_g
		type(grid), pointer :: mesh
		type(mpiControl), pointer :: mpic
		integer :: i,j,k,nx,ny,nz,ierror
		real(DP) :: Q

        mesh => u%ptrMesh_
        mpic => mesh%ptrMPIC_
        
        nx = mesh%nx_
        ny = mesh%ny_
        nz = mesh%nz_ 
    
    	Q = 0.d0
		!$OMP PARALLEL DO DEFAULT(none) &
		!$OMP PRIVATE(i,j,k) &
		!$OMP SHARED(u,mesh,nx,ny,nz) &
		!$OMP REDUCTION(+:Q)
    	do k=1,nz
    		do j=1,ny
    			do i=0,nx-1			
    				
    				Q = Q + u%ux_%f_(i,j,k)*mesh%Vsx_(i,j,k)
    				
    			end do
    		end do
    	end do
    	!$OMP END PARALLEL DO
    	
    	call Mpi_Allreduce(Q, Q_g, 1, MPI_DOUBLE_PRECISION, MPI_SUM,  &
    					   mpic%cartComm_, ierror)
    
    end subroutine
!========================================================================================!

	
end module auxiliaryRoutinesMod


	



