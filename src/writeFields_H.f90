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

!reconstruct and write fields
call writeTimeFolder(runTime,s_nb)

call reconstructAndWriteField(c,gc,runTime%outputFold_)
call reconstructAndWriteField(cs,gcs,runTime%outputFold_)
call reconstructAndWriteField(curv,gcurv,runTime%outputFold_)
call reconstructAndWriteField(p,gp,runTime%outputFold_)
call reconstructAndWriteField(psi,gpsi,runTime%outputFold_)
call reconstructAndWriteFieldV(u,gu,runTime%outputFold_)
!call reconstructAndWriteFieldV(st,gst,runTime%outputFold_)
call reconstructAndWriteFieldV(w,gw,runTime%outputFold_)

!write old time fluxes
call printOldTimeFlux(uEqn,runTime%scheme_,runTime%outputFold_)
	
!write stats
call writeStats(stats,runTime%outputFold_)

!write boxes: make sure every MPI proc can read and write in the main folder.
call printVOFblocks(runTime%outputFold_)
