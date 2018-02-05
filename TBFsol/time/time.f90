! ************************************************************************************** !
!    TBFsol - DNS turbulent bubbly flow solver
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

module timeMod

	use auxiliaryRoutinesMod
	
	implicit none
	
	!enumerate time integration scheme
	integer, parameter :: s_AB2 = 0
	integer, parameter :: s_RK3 = 1
	
	type, public :: time
	
		real(DP) :: t_
		real(DP), private :: Tf_
		real(DP) :: dt_,dtout_,tout_
		integer, private :: writeInterval_
		integer, private :: writeIter_
		
		!adaptive time-step based on CFL
		type(vectorField), pointer, private :: ptrU_ => NULL()
		logical, private :: adaptiveTimeStep_
		real(DP), private :: cflLim_
		real(DP), private :: cflMax_
	
		!timeControl dictionary
		type(dictionary), private :: dict_
		
		!counter time iterations
		integer, private :: iter_
		
		!input-output folder
		integer :: inputFold_,outputFold_
		
		!time integr. scheme
		integer :: scheme_
		
		!rk3 coefficients
		real(DP), dimension(3) :: alpha_
		real(DP), dimension(3) :: gamma_
		real(DP), dimension(3) :: xi_
		integer, private :: rkn_
		
		integer, private :: rkIter_
		
		!redistribution VOF blocks time interval
		real(DP), private :: tVOFB_,dtVOFB_


		contains	
	end type
	
	private :: update 
	private :: info
	private :: initRK3coef
	private :: initTimeLevel
	
	public :: timeCTOR
	public :: timeLoop
	public :: timeOutput
	public :: timeRkStep
	public :: writeTimeFolder
	public :: alphaRKS
	public :: gammaRKS
	public :: xiRKS
	public :: vofBlocksRed

	
contains

!========================================================================================!
    elemental function alphaRKS(this) result(r)
        class(time), intent(in) :: this
        real(DP) :: r
        
        r = this%alpha_(this%rkIter_)
		
    end function
!========================================================================================!

!========================================================================================!
    elemental function gammaRKS(this) result(r)
        class(time), intent(in) :: this
        real(DP) :: r
        
        r = this%gamma_(this%rkIter_)
		
    end function
!========================================================================================!

!========================================================================================!
    elemental function xiRKS(this) result(r)
        class(time), intent(in) :: this
        real(DP) :: r
        
        r = this%xi_(this%rkIter_)
		
    end function
!========================================================================================!

!========================================================================================!
    function vofBlocksRed(this) result(r)
        type(time), intent(inout) :: this
        logical :: r
        
        if (this%t_>=this%tVOFB_) then
        	r = .TRUE.
        	this%tVOFB_ = this%tVOFB_ + this%dtVOFB_
        else
        	r = .FALSE.
        end if
		
    end function
!========================================================================================!

!========================================================================================!
	subroutine timeCTOR(this,u,mpic)
		type(time), intent(out) :: this
		type(vectorField), intent(in), target :: u
		type(mpiControl), intent(in) :: mpic
		
		this%ptrU_ => u
		
		call dictionaryCTOR(this%dict_,'timeControl','specs')

		call readParameter(this%dict_,this%Tf_,'Tf')
		call readParameter(this%dict_,this%dt_,'dt')
		call readParameter(this%dict_,this%inputFold_,'input_folder')
		call readParameter(this%dict_,this%writeInterval_,'writeInterval')
		call readParameter(this%dict_,this%dtout_,'dtout')
		call readParameter(this%dict_,this%adaptiveTimeStep_,'adaptiveTimeStep')
		call readParameter(this%dict_,this%cflLim_,'cflMax')
		call readParameter(this%dict_,this%dtVOFB_,'vofBlocksRedInterval')
		
		this%tout_ = this%dtout_
		this%iter_ = 0
		this%writeIter_ = 0
		this%outputFold_ = this%inputFold_
		
		this%tVOFB_ = this%dtVOFB_

		!init time level
		call initTimeLevel(this,mpic)	
		!init RK parameters
		call initRK3coef(this)
		this%rkIter_ = 0


	end subroutine
!========================================================================================!

!========================================================================================!
    subroutine update(this)
        type(time), intent(inout) :: this
        
        this%iter_ = this%iter_ + 1
        this%writeIter_ = this%writeIter_ + 1
        this%t_ = this%t_ + this%dt_
        
        !compute cfl max
        this%cflMax_ = computeCFLmax(this%ptrU_,this%dt_)
        
        if ((this%iter_>1) .AND. (this%adaptiveTimeStep_)) then
        	this%dt_ = this%dt_*this%cflLim_/this%cflMax_
        	!this%cflMax_ = this%cflLim_
        end if
		
    end subroutine
!========================================================================================!

!========================================================================================!
    function timeLoop(this) result(isRun)
        type(time), intent(inout) :: this
        logical :: isRun
        real(DP) :: small
        
        small = 1.d-9
        
        call update(this)
        
		if (this%t_ >= this%Tf_ + small) then
			isRun = .FALSE.
		else
			call info(this)
			isRun = .TRUE.
		end if
		
    end function
!========================================================================================!

!========================================================================================!
    function timeRkStep(this) result(run)
        type(time), intent(inout) :: this
        logical :: run
        
        this%rkIter_ = this%rkIter_ + 1
        
        if (this%rkIter_ > this%rkn_) then
        	run = .FALSE.
        	this%rkIter_ = 0
        else
        	run = .TRUE.
        	if (IS_MASTER) then
        		if (this%scheme_==s_RK3) then
        			write(*,'(A,I2)') '	rk step: ', this%rkIter_
        		else
        			write(*,'(A)') ' AB2 step'
        		end if
        	end if
        end if
		
    end function
!========================================================================================!

!========================================================================================!
    function timeOutput(this) result(isOutput)
        type(time), intent(inout) :: this
        logical :: isOutput
        
        if ((this%writeInterval_ == this%writeIter_).OR. &
        	(this%t_ >= this%tout_)) then
        	isOutput = .TRUE.
        	this%writeIter_ = 0
        	this%tout_ = this%tout_ + this%dtout_
        else
        	isOutput = .FALSE.
        end if

    end function
!========================================================================================!

!========================================================================================!
    subroutine writeTimeFolder(this)
        type(time), intent(inout) :: this
        character(len=10) :: dirName
        integer :: CSTAT
        
        this%outputFold_ = this%outputFold_ + 1
        
        write(dirName,s_intFormat) this%outputFold_
        
        if (IS_MASTER) then 	
        	
        	!mkdir
        	call execute_command_line('mkdir -p ./' // adjustl(trim(dirName)),CMDSTAT=CSTAT )

        	if ((CSTAT > 0) .OR. (CSTAT < 0)) then
        		call mpiABORT('mkdir on new time level failed ')
        	end if
        	
        	call execute_command_line(adjustl('touch ./'//trim(dirName)//'/time'),CMDSTAT=CSTAT )

        	
        	!write time
        	open(UNIT=s_IOunitNumber,FILE=trim(adjustl(dirName))//'/time',STATUS='REPLACE',ACTION='WRITE')
				write(s_IOunitNumber,s_doubleFormat) this%t_
			close(s_IOunitNumber)

        
        end if
        

    end subroutine
!========================================================================================!

!========================================================================================!
    subroutine info(this)
    	type(time), intent(inout) :: this

		if (IS_MASTER) then
			write(*,*) ''
		
			write(*,'(3(A,'//s_outputFormat(2:9)//'))') 'SOLVING FOR TIME: t = ', this%t_, '  dt = ', this%dt_, &
			  	'  CFL max = ', this%cflMax_
		end if

    end subroutine
!========================================================================================!

!========================================================================================!
    subroutine initRK3coef(this)
    	type(time), intent(inout) :: this
    	type(dictionary) :: dict
    	
    	
    	call dictionaryCTOR(dict,'schemes','specs')
    	call readParameter(dict,this%scheme_,'time_scheme')
    
    	if (this%scheme_==1) then
    	
    		this%rkn_ = 3
    	
        	this%gamma_(1) = 8.d0/15.d0
       	 	this%gamma_(2) = 5.d0/12.d0
       	 	this%gamma_(3) = 3.d0/4.d0
        
        	this%xi_(1) = 0.d0
        	this%xi_(2) = -17.d0/60.d0
        	this%xi_(3) = -5.d0/12.d0
        
        else if (this%scheme_==0) then
        
    		this%rkn_ = 1
    	
        	this%gamma_(1) = 1.5d0
        	this%gamma_(2) = 0.d0
        	this%gamma_(3) = 0.d0
        
        	this%xi_(1) = -0.5d0
        	this%xi_(2) = 0.d0
        	this%xi_(3) = 0.d0
        
        else
        	call mpiABORT('invalid time integration scheme ')
        end if
        
        this%alpha_(1) = this%gamma_(1)+this%xi_(1)
        this%alpha_(2) = this%gamma_(2)+this%xi_(2)
        this%alpha_(3) = this%gamma_(3)+this%xi_(3)
		
    end subroutine
!========================================================================================!

!========================================================================================!
    subroutine initTimeLevel(this,mpic)
    	type(time), intent(inout) :: this
    	type(mpiControl), intent(in) :: mpic
    	character(len=10) :: dirName	
    	integer :: ierror
    	
    	if (IS_MASTER) then
    	
    		write(dirName,s_intFormat) this%inputFold_
    	
    		if (this%inputFold_ == 0) then
    			this%t_ = 0.d0
    		else
        		open(UNIT=s_IOunitNumber,FILE=adjustl(trim(dirName)//'/time'),STATUS='old',ACTION='read')
					read(s_IOunitNumber,s_doubleFormat) this%t_
				close(s_IOunitNumber)    		
    		end if
    	end if
    	
    	call MPI_BCAST(this%t_, 1, MPI_DOUBLE_PRECISION, 0, mpic%cartComm_, ierror)
    	

    end subroutine
!========================================================================================!




	
end module timeMod


	



