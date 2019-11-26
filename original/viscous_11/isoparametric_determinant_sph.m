function sol=isoparametric_determinant_sph(ne,neTotal,nex,ney,n,nx,ny,nz,nxny,nexney,X,Y,Z,n_eSurf,n_nSurf,nx_,parallel_computing)
gps=[(-sqrt(3/5)+1)/2 0.5 (sqrt(3/5)+1)/2];
sol=zeros(neTotal,3,3,3);
if parallel_computing==1
    parfor e=1:1:neTotal
        [eX,eY,eZ]=get_eXYZ_sph(e,ne,nex,ney,nexney,nx,ny,nz,n,nxny,n_eSurf,n_nSurf,nx_,X,Y,Z);
        sol(e,:,:,:)=eDeterminant_3d(eX,eY,eZ,gps);
    end
elseif parallel_computing==0
    for e=1:1:neTotal
        [eX,eY,eZ]=get_eXYZ_sph(e,ne,nex,ney,nexney,nx,ny,nz,n,nxny,n_eSurf,n_nSurf,nx_,X,Y,Z);
        sol(e,:,:,:)=eDeterminant_3d(eX,eY,eZ,gps);
    end
end

end

function determinants=eDeterminant_3d(eX,eY,eZ,gps)
J=zeros(3,3);
determinants=zeros(3,3,3);
for iz=1:1:3
    for iy=1:1:3
        for ix=1:1:3

            J(1,1)=map_to_global_1_compt_3d(gps(ix),gps(iy),gps(iz),eX,[1 0 0]);
            J(1,2)=map_to_global_1_compt_3d(gps(ix),gps(iy),gps(iz),eX,[0 1 0]);
            J(1,3)=map_to_global_1_compt_3d(gps(ix),gps(iy),gps(iz),eX,[0 0 1]);
            J(2,1)=map_to_global_1_compt_3d(gps(ix),gps(iy),gps(iz),eY,[1 0 0]);
            J(2,2)=map_to_global_1_compt_3d(gps(ix),gps(iy),gps(iz),eY,[0 1 0]);
            J(2,3)=map_to_global_1_compt_3d(gps(ix),gps(iy),gps(iz),eY,[0 0 1]);
            J(3,1)=map_to_global_1_compt_3d(gps(ix),gps(iy),gps(iz),eZ,[1 0 0]);
            J(3,2)=map_to_global_1_compt_3d(gps(ix),gps(iy),gps(iz),eZ,[0 1 0]);
            J(3,3)=map_to_global_1_compt_3d(gps(ix),gps(iy),gps(iz),eZ,[0 0 1]);
            
            determinants(ix,iy,iz)=det(J);
            
        end
    end
end

end