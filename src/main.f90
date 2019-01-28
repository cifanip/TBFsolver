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

PROGRAM main

	USE momentumEqnMod
	USE poissonEqnMod
	USE initialConditionsMod
	USE vofMOD
	USE statisticsMod
	USE rampUpPropMod
	
	IMPLICIT NONE
	
	integer :: ierror,flow_solver
	type(mpiControl) :: mpiCTRL
	type(time) :: runTime
	type(grid) :: gMesh, mesh
	type(field) :: gp, p, gpsi, psi
	type(field) :: gc, c, gcs, cs, gcurv, curv
	type(field) :: rho, mu
	type(vfield) :: gU, U, gW, w, gST, st
	type(momentumEqn) :: uEqn
	type(poissonEqn) :: pEqn
	type(VOF) :: vofS
	type(statistics) :: stats
	type(rampUpProp) :: rhoRamp,muRamp
	type(parFile) :: file_fsolver
	real(DP) :: t_S, t_E, t_S0, t_E0


	call MPI_INIT(ierror)
	call mpiGVAR()
	call mpiControlCTOR(mpiCTRL)
	call timeCTOR(runTime,u,mpiCTRL)	
	
	t_S0 = MPI_Wtime()

	INCLUDE 'createFields_H.f90'
	
	call init_main_solver()

	call info_run_start()

	!***************************** FLOW SOLVER ******************************!
	if (flow_solver==SINGLE_PHASE_FLOW) then
		call sph_flow_solver()
	end if 

	if (flow_solver==TWO_PHASE_FLOW) then
		call tph_flow_solver()
	end if		
	!************************************************************************!
	
	!write final before exit
	INCLUDE 'writeFields_H.f90'

	t_E0 = MPI_Wtime()
	
	call info_run_end()
	
	
	!clean up
	call deallocateBlocks(vofBlocks)
	call MPI_FINALIZE(ierror)	

	
100 FORMAT(E18.10)
101 FORMAT(I4)


contains

!========================================================================================!
	subroutine info_run_start()
	
		if (IS_MASTER) then
			write(*,*) ''
			write(*,'(A,'//s_intFormat(2:3)//',A)') &
					'INIT TIME-INTEGRATOR on ', mpiCTRL%nProcs_*N_THREADS, ' cores'
		end if	
	
	end subroutine
!========================================================================================!

!========================================================================================!
	subroutine info_run_end()
	
		if (IS_MASTER) then
			write(*,*) ''
			write(*,'(A,'//s_outputFormat(2:9)//')') &
					'EXIT RUN NORMAL. SIMULATION TIME: ', t_E0-t_S0
		end if		
	
	end subroutine
!========================================================================================!

!========================================================================================!
	subroutine info_run_cpu_time()

		if (IS_MASTER) then
			write(*,'(A,'//s_outputFormat(2:9)//')') '	TOTAL CPU TIME: ', t_E-t_S
		end if
	
    end subroutine
!========================================================================================!

!========================================================================================!
	subroutine init_main_solver()
	
		!read flow solver type
		call parFileCTOR(file_fsolver,'flowSolver','specs')
		call readParameter(file_fsolver,flow_solver,'flow_solver')
	
		call vofCTOR(vofS,gmesh,mesh,runTime,flow_solver)
		call momentumEqnCTOR(uEqn,gMesh,mesh,gu,u,runTime)
#ifdef FAST_MODE
		call poissonEqnCTOR(pEqn,mesh,gMesh,psi,runTime,vofS%rhol_,vofS%rhog_)
#endif
#ifdef MG_MODE
		call poissonEqnCTOR(pEqn,mesh,cs,psi,runTime)
#endif

		!build ramps
		call rampUpPropCTOR(rhoRamp,vofS%rhog_,vofS%rhol_,vofS%rhog_)
		call rampUpPropCTOR(muRamp,vofS%mug_,vofS%mul_,vofS%mug_)
		call updateMaterialProps(vofS,c,cs,rho,mu)

		call statisticsCTOR(stats,u,w,p,c,mu,vofS%mul_/vofS%rhol_,gMesh,runTime)
		
		!store time-step restrictions
		call compute_timestep_restrictions(runTime,u,gMesh,vofS%rhol_,vofS%rhog_,&
									       vofS%mul_,vofS%mug_,vofS%sigma_,flow_solver)
	
	end subroutine
!========================================================================================!

!========================================================================================!
	subroutine reset_matprops_sph()
	
    	c%f_   = 0.d0
    	cs%f_  = 0.d0
    	vofS%rhog_ = vofS%rhol_
    	vofS%mug_ = vofS%mul_
    	rho%f_ = vofS%rhol_
    	mu%f_  = vofS%mul_	
	
	end subroutine
!========================================================================================!

!========================================================================================!
    subroutine tph_flow_solver()

	!cycle time loop
	do while (timeloop(runTime))
	
		t_S = MPI_Wtime()
		
		!update ramped props
		call updateProp(rhoRamp,runTime%t_,vofS%rhog_)
		call updateProp(muRamp,runTime%t_,vofS%mug_)

		!update stats
		call updateStats(stats,runTime%t_,runTime%dt_)
		
		do while (timeRkStep(runTime))	

			!advect VOF
			call solveVOF(vofS,c,u)
			
			!update surface tension and mat pros
			call computeSurfaceTension(vofS,st,curv)			
			call updateMaterialProps(vofS,c,cs,rho,mu)
			
			call setPressGrad(uEqn%flowCtrl_,uEqn%Ret_,rho,vofS%rhol_,vofS%mul_,&
							  uEqn%gCH_,uEqn%fs_)

			!solve momentum equation
			call solveMomentumEqn(uEqn,u,p,mu,rho,st,c)

			!solve poisson equation
			call solvePoissonEqn(pEqn,psi,rho,u)	
			
			!divergence free velocity			
#ifdef FAST_MODE
			call makeVelocityDivFree(uEqn,u,psi,rho,pEqn%rho0_,pEqn%nl_)
#endif
#ifdef MG_MODE
			call makeVelocityDivFree(uEqn,u,psi,rho)
#endif

			!set flow rate
			call setFlowRate(uEqn%flowCtrl_,u,rho,uEqn%Q0_,runTime%dt_,&
					         alphaRKS(runTime))

			!print out continuity error
    		call computeContinuityError(u,runTime%dt_)
			
			!correct pressure
#ifdef FAST_MODE
			call updatePressure(pEqn,p,psi)
#endif
#ifdef MG_MODE
			call updatePressure(p,psi)
#endif
			
		end do	

		call computeVorticity(u,w)

		if (timeOutput(runTime)) then
			INCLUDE 'writeFields_H.f90'
		end if
			
		t_E = MPI_Wtime()
		
		call info_run_cpu_time()
		
		
		!*********
		call mpiAbort('APPOSTO')
		!*********
	end do

    	
    end subroutine
!========================================================================================!

!========================================================================================!
    subroutine sph_flow_solver()
    
	call reset_matprops_sph()
    
	!cycle time loop
	do while (timeloop(runTime))
	
		t_S = MPI_Wtime()

		!update stats
		call updateStats(stats,runTime%t_,runTime%dt_)
		
		do while (timeRkStep(runTime))	
			
			call setPressGrad(uEqn%flowCtrl_,uEqn%Ret_,rho,vofS%rhol_,vofS%mul_,&
							  uEqn%gCH_,uEqn%fs_)

			!solve momentum equation
			call solveMomentumEqn(uEqn,u,p,mu,rho,st,c)

			!solve poisson equation
			call solvePoissonEqn(pEqn,psi,rho,u)	
			
			!divergence free velocity			
#ifdef FAST_MODE
			call makeVelocityDivFree(uEqn,u,psi,rho,pEqn%rho0_,pEqn%nl_)
#endif
#ifdef MG_MODE
			call makeVelocityDivFree(uEqn,u,psi,rho)
#endif

			!set flow rate
			call setFlowRate(uEqn%flowCtrl_,u,rho,uEqn%Q0_,runTime%dt_,&
					         alphaRKS(runTime))

			!print out continuity error
    		call computeContinuityError(u,runTime%dt_)
			
			!correct pressure
#ifdef FAST_MODE
			call updatePressure(pEqn,p,psi)
#endif
#ifdef MG_MODE
			call updatePressure(p,psi)
#endif
			
		end do	

		call computeVorticity(u,w)

		if (timeOutput(runTime)) then
			INCLUDE 'writeFields_H.f90'
		end if
			
		t_E = MPI_Wtime()
		
		call info_run_cpu_time()
		
	end do
    	
    end subroutine
!========================================================================================!
    

END PROGRAM main	

