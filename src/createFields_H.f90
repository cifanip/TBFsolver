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

!*********************** primary fields
!init volume fraction
call read_file_field(gc,c,gMesh,mesh,'c',runTime%inputFold_,halo_size=1)

!init velocity
call read_file_vfield(gu,u,gMesh,mesh,'u',runTime%inputFold_,halo_size=2)

!init phi
call read_file_field(gpsi,psi,gMesh,mesh,'psi',runTime%inputFold_,halo_size=1)

!*********************** derived fields
!init smooth volume fraction
if (IS_MASTER) then
	call fieldCTOR(gcs,'cs',gMesh,'cl',halo_size=1,initOpt=0)
end if
call fieldCTOR(cs,'cs',mesh,'cl',halo_size=1,initOpt=-1)
call decomposeField(gcs,cs)

!init curvature
if (IS_MASTER) then
	call fieldCTOR(gcurv,'k',gMesh,'cl',halo_size=1,initOpt=0)
end if
call fieldCTOR(curv,'k',mesh,'cl',halo_size=1,initOpt=-1)
call decomposeField(gcurv,curv)

!init pressure
if (IS_MASTER) then
	call fieldCTOR(gp,'p',gMesh,'cl',halo_size=1,initOpt=3,nFolder=runTime%inputFold_)
	call copyBoundary(gp,gpsi)
end if
call fieldCTOR(p,'p',mesh,'cl',halo_size=1,initOpt=-1)
call decomposeField(gp,p)

!init vorticity
if (IS_MASTER) then
	call vfieldCTOR(gw,'w',gMesh,'cl','cl','cl',halo_size=1,initOpt=0)
end if
call vfieldCTOR(w,'w',mesh,'cl','cl','cl',halo_size=1,initOpt=-1)
call decomposeFieldV(gw,w)

!init surface tension force
if (IS_MASTER) then
	call vfieldCTOR(gst,'st',gMesh,'sx','sy','sz',halo_size=1,initOpt=0)
end if
call vfieldCTOR(st,'st',mesh,'sx','sy','sz',halo_size=1,initOpt=-1)
call decomposeFieldV(gst,st)

!init material properties fields
call fieldCTOR(rho,'rho',mesh,'cl',halo_size=1,initOpt=-1)
call fieldCTOR(mu,'mu',mesh,'cl',halo_size=1,initOpt=-1)
call copyBoundary(rho,cs)
call copyBoundary(mu,cs)

if (IS_MASTER) then
	write(*,'(A)') 'END BUILD MESH AND FIELDS'
end if
