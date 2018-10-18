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

module timeMod

	use solverTypesMod
	use auxiliaryRoutinesMod
	use initialConditionsMod, only: pi
	
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
		type(vfield), pointer, private :: ptrU_ => NULL()
		logical, private :: adaptiveTimeStep_
		real(DP), private :: cflLim_
		real(DP), private :: cflMax_
		
		!time-step restr. (viscosity, surface tension)
		logical, private :: setTimeStep_
		real(DP), private :: dtLim_
	
		!timeControl parFile
		type(parFile), private :: pfile_
		
		!counter time iterations
		integer :: iter_
		
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
		
		!restart boxes
		logical :: restart_boxes_


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
	public :: compute_timestep_restrictions

	
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
		type(vfield), intent(in), target :: u
		type(mpiControl), intent(in) :: mpic
		
		this%ptrU_ => u
		
		call parFileCTOR(this%pfile_,'timeControl','specs')

		call readParameter(this%pfile_,this%Tf_,'Tf')
		call readParameter(this%pfile_,this%dt_,'dt')
		call readParameter(this%pfile_,this%inputFold_,'input_folder')
		call readParameter(this%pfile_,this%writeInterval_,'writeInterval')
		call readParameter(this%pfile_,this%dtout_,'dtout')
		call readParameter(this%pfile_,this%adaptiveTimeStep_,'adaptiveTimeStep')
		call readParameter(this%pfile_,this%setTimeStep_,'setTimeStep')
		call readParameter(this%pfile_,this%cflLim_,'cflMax')
		call readParameter(this%pfile_,this%dtVOFB_,'vofBlocksRedInterval')
		call readParameter(this%pfile_,this%restart_boxes_,'restart_boxes')
		
		this%iter_ = 0
		this%writeIter_ = 0
		this%outputFold_ = this%inputFold_

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
        real(DP) :: dt_cfl
        
        this%iter_ = this%iter_ + 1
        this%writeIter_ = this%writeIter_ + 1
        
        !compute cfl max
        this%cflMax_ = computeCFLmax(this%ptrU_,this%dt_)
        
        if ((this%adaptiveTimeStep_).AND.(this%iter_>1)) then
        	dt_cfl=this%dt_*this%cflLim_/(this%cflMax_+tiny(0.d0))
        	this%dt_ = min(dt_cfl,this%dtLim_)
        end if
	
	this%t_ = this%t_ + this%dt_
		
    end subroutine
!========================================================================================!

!========================================================================================!
    subroutine compute_timestep_restrictions(this,u,gmesh,rhol,rhog,mul,mug,sigma,solver)
    	type(time), intent(inout) :: this
    	type(vfield), intent(in) :: u
    	type(grid), intent(in) :: gmesh
    	real(DP), intent(in) :: rhol,rhog,mul,mug,sigma
    	real(DP) :: dxm,dym,dzm,d,dt_cfl,dt_nul,dt_nug,dt_sigma,dt_lim
    	integer, intent(in) :: solver
    	integer :: nx,ny,nz
    	
    	nx=gmesh%nx_
    	ny=gmesh%ny_
    	nz=gmesh%nz_
    	
    	dxm=minval(gmesh%dxf_(0:nx+1))
    	dym=minval(gmesh%dyf_(0:ny+1))
    	dzm=minval(gmesh%dzf_(0:nz+1))
    	d=minval((/dxm,dym,dzm/))
    	
    	dt_nul=rhol*(1.d0/6.d0)*d*d/mul
    	dt_nug=rhog*(1.d0/6.d0)*d*d/mug
    	if (solver==TWO_PHASE_FLOW) then
    		dt_sigma=sqrt(d*d*d*(rhol+rhog)/(4.d0*pi*sigma))
    	else
    		dt_sigma=huge(0.d0)
    	end if

    	dt_cfl=compute_dt_CFL(u,this%cflLim_)
    	
    	dt_nul=dt_nul/2.d0
    	dt_nug=dt_nug/2.d0
    	dt_sigma=dt_sigma/2.5d0
    	
    	dt_lim=minval((/dt_nul,dt_nug,dt_sigma/))
		
		if (this%setTimeStep_) then
			this%dtLim_=dt_lim
			this%dt_=min(this%dtLim_,dt_cfl)
		else
			this%dtLim_=this%dt_
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
			this%t_=this%t_-this%dt_
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
        			write(*,'(A,I2)') ' 	RK STEP: ', this%rkIter_
        		else
        			write(*,'(A)') ' 	AB2 STEP'
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
    subroutine writeTimeFolder(this,nb)
        type(time), intent(inout) :: this
        integer, intent(in) :: nb
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
        	
        	call execute_command_line(adjustl('touch ./'//trim(dirName)//'/info_restart'),CMDSTAT=CSTAT )

        	
        	!write time and total number of bubbles
        	open(UNIT=s_IOunitNumber,FILE=trim(adjustl(dirName))//'/info_restart',&
        		 STATUS='REPLACE',ACTION='WRITE')
				write(s_IOunitNumber,s_doubleFormat) this%t_
				write(s_IOunitNumber,s_intFormat) nb
			close(s_IOunitNumber)

        
        end if
        

    end subroutine
!========================================================================================!

!========================================================================================!
    subroutine info(this)
    	type(time), intent(inout) :: this

		if (IS_MASTER) then
			write(*,*) ''
			
			write(*,*) 'SOLVING FOR TIME:'
			
			write(*,'(A,'//s_doubleFormat(2:10)//')') '	t  =', this%t_
			write(*,'(A,'//s_outputFormat(2:9)//')') '	dt =  ', this%dt_
			write(*,'(A,'//s_outputFormat(2:9)//')') '	CFL max = ', this%cflMax_

		end if

    end subroutine
!========================================================================================!

!========================================================================================!
    subroutine initRK3coef(this)
    	type(time), intent(inout) :: this
    	type(parFile) :: dict
    	
    	
    	call parFileCTOR(dict,'schemes','specs')
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
    			this%tout_ = this%dtout_
    			this%tVOFB_ = this%dtVOFB_
    		else
        		open(UNIT=s_IOunitNumber,FILE=adjustl(trim(dirName)//'/info_restart'),&
        			 STATUS='old',ACTION='read')
					read(s_IOunitNumber,s_doubleFormat) this%t_
				close(s_IOunitNumber)  
				this%tout_ = this%dtout_+this%t_  
				this%tVOFB_ = this%dtVOFB_+this%t_	
    		end if
    	end if
    	
    	call MPI_BCAST(this%t_, 1, MPI_DOUBLE_PRECISION, 0, mpic%cartComm_, ierror)
    	call MPI_BCAST(this%tout_, 1, MPI_DOUBLE_PRECISION, 0, mpic%cartComm_, ierror)
    	call MPI_BCAST(this%tVOFB_, 1, MPI_DOUBLE_PRECISION, 0, mpic%cartComm_, ierror)
    	

    end subroutine
!========================================================================================!




	
end module timeMod


	



