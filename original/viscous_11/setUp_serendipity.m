function [nx,ny,nz,n,nfour]=setUp_serendipity(nex,ney,nez)

nx=nex+1;
ny=ney+1;
nz=nez+1;

n=nx*ny*nz;
nfour=n*4;

end