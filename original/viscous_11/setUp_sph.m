function [nxny,nx_,ny_,n_nSurf,n_eSurf,neTotal,nTotal,nnz_,neightTotal]=setUp_sph(...
    n,nx,ny,nz,ne,nex,ney,nez)
nxny=nx*ny;
nx_=nx-2;
ny_=ny-2;
[n_nSurf,n_eSurf]=get_n_nSurf_n_eSurf(nex,ney,nez,ny,nz,nx_,ny_);
neTotal=ne+n_eSurf(6);
nTotal=n+n_nSurf(6);
nnz_=nnz_sj_sph(nx,ny,nz,n_nSurf);
neightTotal=nTotal*8;
end