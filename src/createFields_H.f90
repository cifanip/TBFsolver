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

!init mesh and global fields
if (IS_MASTER) then
	write(*,'(A)') 'BUILD MESH AND FIELDS'
end if

!build mesh
if (IS_MASTER) then
	call gridCTOR(gMesh,mpiCTRL)
end if
call decomposeGrid(gMesh,mesh,mpiCTRL)
call broadCastGlobalGrid(gMesh,mesh) 

!init volume fraction
if (IS_MASTER) then
	call scalarFieldCTOR(gc,'c',gMesh,'cl',1,initOpt=1,nFolder=runTime%inputFold_)
end if
call scalarFieldCTOR(c,'c',mesh,'cl',1,initOpt=-1)
call decomposeField(gc,c)

!init smooth volume fraction
if (IS_MASTER) then
	call scalarFieldCTOR(gcs,'cs',gMesh,'cl',1,initOpt=0)
end if
call scalarFieldCTOR(cs,'cs',mesh,'cl',1,initOpt=-1)
call decomposeField(gcs,cs)

!init curvature
if (IS_MASTER) then
	call scalarFieldCTOR(gcurv,'k',gMesh,'cl',1,initOpt=0)
end if
call scalarFieldCTOR(curv,'k',mesh,'cl',1,initOpt=-1)
call decomposeField(gcurv,curv)

!init pressure
if (IS_MASTER) then
	call scalarFieldCTOR(gp,'p',gMesh,'cl',1,initOpt=1,nFolder=runTime%inputFold_)
end if
call scalarFieldCTOR(p,'p',mesh,'cl',1,initOpt=-1)
call decomposeField(gp,p)

!init phi (delta pressure)
if (IS_MASTER) then
	call scalarFieldCTOR(gpsi,'psi',gMesh,'cl',1,initOpt=1,nFolder=runTime%inputFold_)
end if
call scalarFieldCTOR(psi,'psi',mesh,'cl',1,initOpt=-1)
call decomposeField(gpsi,psi)

!init velocity
if (IS_MASTER) then
	call vectorFieldCTOR(gu,'u',gMesh,'sx','sy','sz',2,initOpt=1,nFolder=runTime%inputFold_)
	!********************************!
	!call initChFlowVelocity(gu,gMesh)
	!********************************!
end if
call vectorFieldCTOR(u,'u',mesh,'sx','sy','sz',2,initOpt=-1)
call decomposeFieldV(gu,u)

!init vorticity
if (IS_MASTER) then
	call vectorFieldCTOR(gw,'w',gMesh,'cl','cl','cl',1,initOpt=0)
end if
call vectorFieldCTOR(w,'w',mesh,'cl','cl','cl',1,initOpt=-1)
call decomposeFieldV(gw,w)

!init surface tension force
if (IS_MASTER) then
	call vectorFieldCTOR(gst,'st',gMesh,'sx','sy','sz',1,initOpt=0)
end if
call vectorFieldCTOR(st,'st',mesh,'sx','sy','sz',1,initOpt=-1)
call decomposeFieldV(gst,st)

!init material properties fields
call scalarFieldCTOR(rho,'rho',mesh,'cl',1,initOpt=-1)
call scalarFieldCTOR(mu,'mu',mesh,'cl',1,initOpt=-1)
call copyBoundary(rho,cs)
call copyBoundary(mu,cs)

if (IS_MASTER) then
	write(*,'(A)') 'END BUILD MESH AND FIELDS'
end if

