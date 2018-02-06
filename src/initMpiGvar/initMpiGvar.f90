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

module initMpiGvarMod

	implicit none
	
	logical, protected :: IS_MASTER, IS_PAR 
	integer, protected :: N_THREADS
	
	INCLUDE 'mpif.h'
	INCLUDE 'omp_lib.h'
	
	public :: mpiGVAR
	
contains

!========================================================================================!
	subroutine mpiGVAR() 
		integer :: rank, nProcs, ierror

        call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror)
        call MPI_COMM_SIZE(MPI_COMM_WORLD, nProcs, ierror)
        
        !set isMaster
        if (rank == 0) then
        	IS_MASTER = .TRUE.
        else
        	IS_MASTER = .FALSE.
        end if
        
        !set isParallel
        if (nProcs == 1) then
        	IS_PAR = .FALSE.
        else
        	IS_PAR = .TRUE.
        end if
        
        !set number of threads
		!$OMP PARALLEL DEFAULT(none) &
		!$OMP SHARED(N_THREADS)
		 N_THREADS = omp_get_num_threads() 
		!$OMP END PARALLEL	
		
	end subroutine
!========================================================================================!

	
end module initMpiGvarMod


	



