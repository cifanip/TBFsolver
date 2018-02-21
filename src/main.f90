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
	INTEGER :: ierror
	TYPE(mpiControl) :: mpiCTRL
	type(time) :: runTime
	type(grid) :: gMesh, mesh
	type(scalarField) :: gp, p, gpsi, psi
	type(scalarField) :: gc, c, gcs, cs, gcurv, curv
	type(scalarField) :: rho, mu
	type(vectorField) :: gU, U, gW, w, gST, st
	type(momentumEqn) :: uEqn
	type(poissonEqn) :: pEqn
	type(VOF) :: vofS
	type(statistics) :: stats
	type(rampUpProp) :: rhoRamp,muRamp
	real(DP) :: t_S, t_E, t_S0, t_E0


	call MPI_INIT(ierror)
	call mpiGVAR()
	call mpiControlCTOR(mpiCTRL)
	call timeCTOR(runTime,u,mpiCTRL)	
	
	t_S0 = MPI_Wtime()

	INCLUDE 'createFields_H.f90'
	
	call vofCTOR(vofS,gmesh,mesh,runTime)
	call momentumEqnCTOR(uEqn,gMesh,mesh,u,runTime)
!DIR$ IF DEFINED (FAST_MODE)	
	call poissonEqnCTOR(pEqn,mesh,gMesh,psi,runTime,vofS%rhol_,vofS%rhog_)
!DIR$ ELSEIF DEFINED (MG_MODE)
	call poissonEqnCTOR(pEqn,mesh,cs,psi,runTime)
!DIR$ ENDIF
	
	!build ramps
	call rampUpPropCTOR(rhoRamp,vofS%rhog_,vofS%rhol_,vofS%rhog_)
	call rampUpPropCTOR(muRamp,vofS%mug_,vofS%mul_,vofS%mug_)
	call updateMaterialProps(vofS,c,cs,rho,mu)

	call statisticsCTOR(stats,u,w,p,c,vofS%mul_/vofS%rhol_,gMesh)

	
	if (IS_MASTER) then
		write(*,*) ''
		write(*,'(A,'//s_intFormat(2:3)//',A)') &
				'INIT TIME-INTEGRATOR on ', mpiCTRL%nProcs_*N_THREADS, ' cores'
	end if

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
			
			call setPressGrad(uEqn%flowCtrl_,cs,rho,vofS%rhog_,vofS%rhol_,&
			 		          vofS%mul_,uEqn%gCH_,uEqn%fs_)

			!solve momentum equation
			call solveMomentumEqn(uEqn,u,p,mu,rho,st,c)

			!solve poisson equation
			call solvePoissonEqn(pEqn,psi,rho,u)	
			
			!divergence free velocity			
!DIR$ IF DEFINED (FAST_MODE)
			call makeVelocityDivFree(uEqn,u,psi,rho,pEqn%rho0_,pEqn%nl_)
!DIR$ ELSEIF DEFINED (MG_MODE)
			call makeVelocityDivFree(uEqn,u,psi,rho)
!DIR$ ENDIF

			!set flow rate
			call setFlowRate(uEqn%flowCtrl_,u,rho,uEqn%Q0_,runTime%dt_,&
					         alphaRKS(runTime))

			!print out continuity error
    		call computeContinuityError(u,runTime%dt_)
			
			!correct pressure
!DIR$ IF DEFINED (FAST_MODE)
			call updatePressure(pEqn,p,psi)
!DIR$ ELSEIF DEFINED (MG_MODE)
			call updatePressure(p,psi)
!DIR$ ENDIF
			
		end do	

		!call computeVorticity(u,w)

		if (timeOutput(runTime)) then
			INCLUDE 'writeFields_H.f90'
		end if
			
		t_E = MPI_Wtime()
		
		if (IS_MASTER) then
			write(*,'(A,'//s_outputFormat(2:9)//')') '	TOTAL CPU TIME: ', t_E-t_S
		end if
		
	end do
	
	!write final before exit
	INCLUDE 'writeFields_H.f90'

	t_E0 = MPI_Wtime()
	
	
	if (IS_MASTER) then
		write(*,*) ''
		write(*,'(A,'//s_outputFormat(2:9)//')') &
				'EXIT RUN NORMAL. SIMULATION TIME: ', t_E0-t_S0
		!write(*,'(A)'), 'EXIT RUN NORMAL.'
	end if
	
	
	!clean up
	call deallocateBlocks(vofBlocks)
	call MPI_FINALIZE(ierror)	

	
100 FORMAT(E18.10)
101 FORMAT(I4)


contains


END PROGRAM main	

