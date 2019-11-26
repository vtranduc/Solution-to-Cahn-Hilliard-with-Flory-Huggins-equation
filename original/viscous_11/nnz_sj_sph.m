function nnz_=nnz_sj_sph(nx,ny,nz,n_nSurf)
nnz_=nnz_sj_3d(nx,ny,nz)+1728*(n_nSurf(6)-8)+10752;
end