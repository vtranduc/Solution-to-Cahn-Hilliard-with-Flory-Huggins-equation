function sol=isoparametric_determinant_surface_sph(...
    ne,nex,ney,n,nx,ny,nz,nxny,nexney,X,Y,Z,n_eSurf,n_nSurf,...
    nx_,parallel_computing)
gps=[(-sqrt(3/5)+1)/2 0.5 (sqrt(3/5)+1)/2];
sol=zeros(n_eSurf(6),3,3);

if parallel_computing==1
    parfor e=1:1:n_eSurf(1)
        [~,eY,eZ]=get_eXYZ_sph(e+ne,ne,nex,ney,nexney,nx,ny,nz,n,nxny,n_eSurf,n_nSurf,nx_,X,Y,Z);
        sol(e,:,:)=eDeterminant_1_sph(eY,eZ,gps);
    end
    parfor e=1+n_eSurf(1):1:n_eSurf(2)
        [~,eY,eZ]=get_eXYZ_sph(e+ne,ne,nex,ney,nexney,nx,ny,nz,n,nxny,n_eSurf,n_nSurf,nx_,X,Y,Z);
        sol(e,:,:)=eDeterminant_2_sph(eY,eZ,gps);
    end
    parfor e=1+n_eSurf(2):1:n_eSurf(3)
        [eX,~,eZ]=get_eXYZ_sph(e+ne,ne,nex,ney,nexney,nx,ny,nz,n,nxny,n_eSurf,n_nSurf,nx_,X,Y,Z);
        sol(e,:,:)=eDeterminant_3_sph(eX,eZ,gps);
    end
    parfor e=1+n_eSurf(3):1:n_eSurf(4)
        [eX,~,eZ]=get_eXYZ_sph(e+ne,ne,nex,ney,nexney,nx,ny,nz,n,nxny,n_eSurf,n_nSurf,nx_,X,Y,Z);
        sol(e,:,:)=eDeterminant_4_sph(eX,eZ,gps);
    end
    parfor e=1+n_eSurf(4):1:n_eSurf(5)
        [eX,eY,~]=get_eXYZ_sph(e+ne,ne,nex,ney,nexney,nx,ny,nz,n,nxny,n_eSurf,n_nSurf,nx_,X,Y,Z);
        sol(e,:,:)=eDeterminant_5_sph(eX,eY,gps);
    end
    parfor e=1+n_eSurf(5):1:n_eSurf(6)
        [eX,eY,~]=get_eXYZ_sph(e+ne,ne,nex,ney,nexney,nx,ny,nz,n,nxny,n_eSurf,n_nSurf,nx_,X,Y,Z);
        sol(e,:,:)=eDeterminant_6_sph(eX,eY,gps);
    end
elseif parallel_computing==0
    for e=1:1:n_eSurf(1)
        [~,eY,eZ]=get_eXYZ_sph(e+ne,ne,nex,ney,nexney,nx,ny,nz,n,nxny,n_eSurf,n_nSurf,nx_,X,Y,Z);
        sol(e,:,:)=eDeterminant_1_sph(eY,eZ,gps);
    end
    for e=1+n_eSurf(1):1:n_eSurf(2)
        [~,eY,eZ]=get_eXYZ_sph(e+ne,ne,nex,ney,nexney,nx,ny,nz,n,nxny,n_eSurf,n_nSurf,nx_,X,Y,Z);
        sol(e,:,:)=eDeterminant_2_sph(eY,eZ,gps);
    end
    for e=1+n_eSurf(2):1:n_eSurf(3)
        [eX,~,eZ]=get_eXYZ_sph(e+ne,ne,nex,ney,nexney,nx,ny,nz,n,nxny,n_eSurf,n_nSurf,nx_,X,Y,Z);
        sol(e,:,:)=eDeterminant_3_sph(eX,eZ,gps);
    end
    for e=1+n_eSurf(3):1:n_eSurf(4)
        [eX,~,eZ]=get_eXYZ_sph(e+ne,ne,nex,ney,nexney,nx,ny,nz,n,nxny,n_eSurf,n_nSurf,nx_,X,Y,Z);
        sol(e,:,:)=eDeterminant_4_sph(eX,eZ,gps);
    end
    for e=1+n_eSurf(4):1:n_eSurf(5)
        [eX,eY,~]=get_eXYZ_sph(e+ne,ne,nex,ney,nexney,nx,ny,nz,n,nxny,n_eSurf,n_nSurf,nx_,X,Y,Z);
        sol(e,:,:)=eDeterminant_5_sph(eX,eY,gps);
    end
    for e=1+n_eSurf(5):1:n_eSurf(6)
        [eX,eY,~]=get_eXYZ_sph(e+ne,ne,nex,ney,nexney,nx,ny,nz,n,nxny,n_eSurf,n_nSurf,nx_,X,Y,Z);
        sol(e,:,:)=eDeterminant_6_sph(eX,eY,gps);
    end
end

end

function determinants=eDeterminant_1_sph(eY,eZ,gps)
J=zeros(2,2);
determinants=zeros(3,3);
for iz=1:1:3
    for iy=1:1:3
        J(1,1)=map_to_global_1_compt_3d(0,gps(iy),gps(iz),eY,[0 1 0]);
        J(1,2)=map_to_global_1_compt_3d(0,gps(iy),gps(iz),eY,[0 0 1]);
        J(2,1)=map_to_global_1_compt_3d(0,gps(iy),gps(iz),eZ,[0 1 0]);
        J(2,2)=map_to_global_1_compt_3d(0,gps(iy),gps(iz),eZ,[0 0 1]);
        determinants(iy,iz)=det(J);
    end
end
end

function determinants=eDeterminant_2_sph(eY,eZ,gps)
J=zeros(2,2);
determinants=zeros(3,3);
for iz=1:1:3
    for iy=1:1:3
        J(1,1)=map_to_global_1_compt_3d(1,gps(iy),gps(iz),eY,[0 1 0]);
        J(1,2)=map_to_global_1_compt_3d(1,gps(iy),gps(iz),eY,[0 0 1]);
        J(2,1)=map_to_global_1_compt_3d(1,gps(iy),gps(iz),eZ,[0 1 0]);
        J(2,2)=map_to_global_1_compt_3d(1,gps(iy),gps(iz),eZ,[0 0 1]);
        determinants(iy,iz)=det(J);
    end
end
end

function determinants=eDeterminant_3_sph(eX,eZ,gps)
J=zeros(2,2);
determinants=zeros(3,3);
for iz=1:1:3
    for ix=1:1:3
        J(1,1)=map_to_global_1_compt_3d(gps(ix),0,gps(iz),eX,[1 0 0]);
        J(1,2)=map_to_global_1_compt_3d(gps(ix),0,gps(iz),eX,[0 0 1]);
        J(2,1)=map_to_global_1_compt_3d(gps(ix),0,gps(iz),eZ,[1 0 0]);
        J(2,2)=map_to_global_1_compt_3d(gps(ix),0,gps(iz),eZ,[0 0 1]);
        determinants(ix,iz)=det(J);
    end
end
end

function determinants=eDeterminant_4_sph(eX,eZ,gps)
J=zeros(2,2);
determinants=zeros(3,3);
for iz=1:1:3
    for ix=1:1:3
        J(1,1)=map_to_global_1_compt_3d(gps(ix),1,gps(iz),eX,[1 0 0]);
        J(1,2)=map_to_global_1_compt_3d(gps(ix),1,gps(iz),eX,[0 0 1]);
        J(2,1)=map_to_global_1_compt_3d(gps(ix),1,gps(iz),eZ,[1 0 0]);
        J(2,2)=map_to_global_1_compt_3d(gps(ix),1,gps(iz),eZ,[0 0 1]);
        determinants(ix,iz)=det(J);
    end
end
end

function determinants=eDeterminant_5_sph(eX,eY,gps)
J=zeros(2,2);
determinants=zeros(3,3);
for iy=1:1:3
    for ix=1:1:3
        J(1,1)=map_to_global_1_compt_3d(gps(ix),gps(iy),0,eX,[1 0 0]);
        J(1,2)=map_to_global_1_compt_3d(gps(ix),gps(iy),0,eX,[0 1 0]);
        J(2,1)=map_to_global_1_compt_3d(gps(ix),gps(iy),0,eY,[1 0 0]);
        J(2,2)=map_to_global_1_compt_3d(gps(ix),gps(iy),0,eY,[0 1 0]);
        determinants(ix,iy)=det(J);
    end
end
end

function determinants=eDeterminant_6_sph(eX,eY,gps)
J=zeros(2,2);
determinants=zeros(3,3);
for iy=1:1:3
    for ix=1:1:3
        J(1,1)=map_to_global_1_compt_3d(gps(ix),gps(iy),1,eX,[1 0 0]);
        J(1,2)=map_to_global_1_compt_3d(gps(ix),gps(iy),1,eX,[0 1 0]);
        J(2,1)=map_to_global_1_compt_3d(gps(ix),gps(iy),1,eY,[1 0 0]);
        J(2,2)=map_to_global_1_compt_3d(gps(ix),gps(iy),1,eY,[0 1 0]);
        determinants(ix,iy)=det(J);
    end
end
end