function [nx,ny,nz,n,neight]=setUp_3d(nex,ney,nez)

nx=nex+1;
ny=ney+1;
nz=nez+1;

n=nx*ny*nz;
neight=n*8;

end