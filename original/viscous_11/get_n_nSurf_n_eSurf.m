function [n_nSurf,n_eSurf]=get_n_nSurf_n_eSurf(nex,ney,nez,ny,nz,nx_,ny_)
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
n_eSurf_x=ney*nez;
n_eSurf_y=nex*nez;
n_eSurf_z=nex*ney;
n_eSurf(1)=n_eSurf_x;
n_eSurf(2)=n_eSurf_x*2;
n_eSurf(3)=n_eSurf(2)+n_eSurf_y;
n_eSurf(4)=n_eSurf(3)+n_eSurf_y;
n_eSurf(5)=n_eSurf(4)+n_eSurf_z;
n_eSurf(6)=n_eSurf(5)+n_eSurf_z;
end