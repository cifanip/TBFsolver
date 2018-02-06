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

module multiGridMod
	
	use rbgsMod
	use poissMatMod
	
	implicit none


	type, public :: multiGrid
		
		type(dictionary), private :: dict_ 
		
		integer :: nLevels_
		integer, private :: nPreSweep_
		integer, private :: nPostSweep_
		
		real(DP), private :: tol_
		integer :: iter_
		integer, private :: maxIter_
		
		!keep a pointer to grid
		type(grid), pointer, private :: ptrMesh_ => NULL()
		
		type(rbgs) :: smoother_
		type(poissMat), private :: directS_
		
		logical, private :: fullInfo_

	end type
	
	private :: checkLevelsNumber
	private :: initDirectSolver
	private :: mgVcycle
	private :: continueIterating
	
	public :: multiGridCTOR
	public :: solveMG
	public :: residualMG

  	
contains

!========================================================================================!
    elemental function residualMG(this) result(r)
        class(multiGrid), intent(in) :: this
		real(DP) :: r
		
		r = this%smoother_%res_

    end function
!========================================================================================!

!========================================================================================!
	subroutine multiGridCTOR(this,mesh,p,beta,s)
		type(multiGrid) :: this
		type(scalarField), intent(inout) :: p, beta, s
		type(grid), target :: mesh
		
		this%ptrMesh_ => mesh
		
		call dictionaryCTOR(this%dict_,'pcg_solver','specs')
		
		call readParameter(this%dict_,this%nLevels_,'nLevels')
		call readParameter(this%dict_,this%nPreSweep_,'nPreSweep')
		call readParameter(this%dict_,this%nPostSweep_,'nPostSweep')
		call readParameter(this%dict_,this%tol_,'tolMG')
		call readParameter(this%dict_,this%maxIter_,'maxIterMG')
		call readParameter(this%dict_,this%fullInfo_,'fullInfoMG')
		
		call checkLevelsNumber(this)
		
		!init smoother
		call rbgsCTOR(this%smoother_,mesh,p)
						
		!init coarse meshes and fields
		call coarsenGrids(mesh,this%nLevels_)
		call coarsenFields(p,this%nLevels_)
		call coarsenFields(beta,this%nLevels_)
		call coarsenFields(s,this%nLevels_)
		
		!set mg levels
		call setMgLevels(this%ptrMesh_,this%nLevels_)
		
		!init coarse Gauss-Seidel solvers
		call coarsenRbgsSolvers(this%smoother_,mesh,p,this%nLevels_)
		
		!init direct solver for coarsest level
		call initDirectSolver(this,mesh,p)
		
		!homogeneous b.c. for error fields
		call resetBCerrorField(p)
		
					    
	end subroutine
!========================================================================================!

!========================================================================================!
    subroutine checkLevelsNumber(this)
        type(multiGrid), intent(inout) :: this
        integer :: r, i

        !check direction-wise
        !x
        do i=1,this%nLevels_-1
        	r = mod(this%ptrMesh_%nx_,2**i)
        	if (r /= 0) then
        		call mpiABORT('Attempt to initialize a MG grid with a non-integer number of cells along x') 
        	end if
        end do
        
        !y
        do i=1,this%nLevels_-1
        	r = mod(this%ptrMesh_%ny_,2**i)
        	if (r /= 0) then
        		call mpiABORT('Attempt to initialize a MG grid with a non-integer number of cells along y') 
        	end if
        end do
        
        !z
        do i=1,this%nLevels_-1
        	r = mod(this%ptrMesh_%nz_,2**i)
        	if (r /= 0) then
        		call mpiABORT('Attempt to initialize a MG grid with a non-integer number of cells along z') 
        	end if
        end do
        
    end subroutine
!========================================================================================!

!========================================================================================!
    recursive subroutine initDirectSolver(this,mesh,p)
        type(multiGrid), intent(inout) :: this
        type(grid), intent(in) :: mesh
		type(scalarField), intent(in) :: p
		
		if ( mesh%level_ == 1) then
			 call poissMatCTOR(this%directS_,mesh,p)
		else
			call initDirectSolver(this,mesh%ptrGrid_,p%ptrf_)
		end if

    end subroutine
!========================================================================================!


!========================================================================================!
    subroutine solveMG(this,mesh,p,beta,s)
        type(multiGrid), intent(inout) :: this
        type(grid), intent(in) :: mesh
		type(scalarField), intent(inout) :: p
		type(scalarField), intent(inout) :: beta, s
		
		this%iter_ = 0
		
		!compute initial residual
		call computeResiduals(this%smoother_,p,beta,s)
		
    	do while(continueIterating(this))
    	
    		call mgVcycle(this,mesh,p,beta,s,this%smoother_)
    		
    	end do
    
    end subroutine
!========================================================================================!

!========================================================================================!
    recursive subroutine mgVcycle(this,mesh,p,beta,s,smoother)
        type(multiGrid), intent(inout) :: this
        type(grid), intent(in) :: mesh
		type(scalarField), intent(inout) :: p
		type(scalarField), intent(inout) :: beta, s
		type(rbgs), intent(inout) :: smoother
		
		
		if (mesh%level_ == 1) then
			!direct solver for coarsest level
			call solveMat(this%directS_,p,beta,s)			
		else 
		
			!pre-smoothing
			if (mesh%level_ == this%nLevels_) then
				call solveRBGS(smoother,p,beta,s,this%nPreSweep_,isToBeReset=.FALSE.)			
			else
				call solveRBGS(smoother,p,beta,s,this%nPreSweep_,isToBeReset=.TRUE.)
			end if
			
			
			!residual restriction
			call restrictField(s,smoother%resV_)	
			!beta field restriction
			call restrictField(beta)
			

			!recursive call to coarser level
			call mgVcycle(this,mesh%ptrGrid_,		&
							   p%ptrf_,				&
							   beta%ptrf_,			&
							   s%ptrf_,				&
							   smoother%ptrRbgs_)
							   						   

			!coarse grid correction
			call prolongateField(p,p%prol_)
			call unarySum_omp(p%f_,p%prol_,1,p%nx_,1,p%ny_,1,p%nz_)
			call updateBoundaries(p)

			!post smoothing
			call solveRBGS(smoother,p,beta,s,this%nPostSweep_,isToBeReset=.FALSE.)
		
		end if
		

    end subroutine
!========================================================================================!


!========================================================================================!
    function continueIterating(this) RESULT(isVar)
        type(multiGrid), intent(inout) :: this
        logical :: isVar
        
        
	!check max iteration limit
	if (this%iter_ == this%maxIter_) then
		if (IS_MASTER) then 
			write(*,*) '!*************** WARNING *****************'
			write(*,*) 'EXIT multiGrid iterations: max iter reached'
		end if
		isVar = .FALSE.
		return
	end if
	
	!check tolerance
	if (this%smoother_%res_ > this%tol_) then
	
		if (IS_MASTER) then 
			if (this%fullInfo_) then
				write(*,*) 'multiGrid solver: iteration ', this%iter_, ' residual = ', &
				this%smoother_%res_
			end if
		end if
	   isVar = .TRUE.
	   
	else
	
		if (IS_MASTER) then
			if (this%fullInfo_) then
				write(*,*) 'Criterion met at iteration: ', this%iter_, ' residual = ', &
				this%smoother_%res_
			end if
		end if
		
		isVar = .FALSE.
		
	end if
	
	this%iter_ = this%iter_ + 1
		

    end function
!========================================================================================!





end module multiGridMod

