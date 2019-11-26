function test60

clear
clc

nex=5;
ney=4;
nez=3;

%--- Set Up --------------------------------

nx=nex+1;
ny=ney+1;
nz=nez+1;

n=nx*ny*nz;

nxny=nx*ny;
nexney=nex*ney;

ne=nex*ney*nez;

%---Unique for sphere---
nx_=nx-2;
ny_=ny-2;
nz_=nz-2;
n_nodesurf=2*(nx_*ny_+nx_*nz_+ny_*nz_)+4*(nx_+ny_+nz_)+8;

n_nodessurf_x=ny*nz;
n_nodessurf_y=nx_*nz;
n_nodessurf_z=nx_*ny_;

n_nSurf=zeros(1,6);
n_nSurf(1)=n_nodessurf_x;
n_nSurf(2)=2*n_nodessurf_x;
n_nSurf(3)=n_nSurf(2)+n_nodessurf_y;
n_nSurf(4)=n_nSurf(3)+n_nodessurf_y;
n_nSurf(5)=n_nSurf(4)+n_nodessurf_z;
n_nSurf(6)=n_nSurf(5)+n_nodessurf_z;

if n_nSurf(6)~=n_nodesurf
    error('dafafafdagasd')
end

n_eSurf=zeros(1,6);

n_eSurf_x=ney*nez;
n_eSurf_y=nex*nez;
n_eSurf_z=nex*ney;

n_eSurf(1)=n_eSurf_x;
n_eSurf(2)=n_eSurf_x*2;
n_eSurf(3)=n_eSurf(2)+n_eSurf_y;
n_eSurf(4)=n_eSurf(3)+n_eSurf_y;
n_eSurf(5)=n_eSurf(4)+n_eSurf_z;
n_eSurf(6)=n_eSurf(5)+n_eSurf_z;

n_eSurf;


neight=8*n;

neightTotal=neight+8*n_nSurf(6);

nTotal=n+n_nSurf(6);

neightTotal=8*nTotal;

%-----------------------------

e=75

nodes=get_nodes_of_element_sphere(e,ne,nex,ney,nexney,nx,ny,nz,n,nxny,n_eSurf,n_nSurf,nx_);

nodes'









end