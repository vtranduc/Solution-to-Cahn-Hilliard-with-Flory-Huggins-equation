function weights=isoparametrically_adjust_weights_surface_sph(...
    ne,nex,ney,n,nx,ny,nz,nxny,nexney,X,Y,Z,parallel_computing,...
    n_eSurf,n_nSurf,nx_)

weights=zeros(n_eSurf(6),8,8,7,3,3);

gps=[(-sqrt(3/5)+1)/2 0.5 (sqrt(3/5)+1)/2];
w=[5/18 4/9 5/18];
orientation_list=...
    [0 0 0
    1 0 0
    0 1 0
    1 1 0
    0 0 1
    1 0 1
    0 1 1
    1 1 1];

%Work on each surface

if parallel_computing==1
    parfor e=ne+1:1:ne+n_eSurf(1)
        [eX,eY,eZ]=get_eXYZ_sph(e,ne,nex,ney,nexney,nx,ny,nz,n,nxny,n_eSurf,n_nSurf,nx_,X,Y,Z);
        weights(e-ne,:,:,:,:,:)=isoparametric_surface_minus_x_sph(gps,w,eX,eY,eZ,orientation_list);
    end
    parfor e=ne+1+n_eSurf(1):1:ne+n_eSurf(2)
        [eX,eY,eZ]=get_eXYZ_sph(e,ne,nex,ney,nexney,nx,ny,nz,n,nxny,n_eSurf,n_nSurf,nx_,X,Y,Z);
        weights(e-ne,:,:,:,:,:)=isoparametric_surface_plus_x_sph(gps,w,eX,eY,eZ,orientation_list);
    end
    parfor e=ne+1+n_eSurf(2):1:ne+n_eSurf(3)
        [eX,eY,eZ]=get_eXYZ_sph(e,ne,nex,ney,nexney,nx,ny,nz,n,nxny,n_eSurf,n_nSurf,nx_,X,Y,Z);
        weights(e-ne,:,:,:,:,:)=isoparametric_surface_minus_y_sph(gps,w,eX,eY,eZ,orientation_list);
    end
    parfor e=ne+1+n_eSurf(3):1:ne+n_eSurf(4)
        [eX,eY,eZ]=get_eXYZ_sph(e,ne,nex,ney,nexney,nx,ny,nz,n,nxny,n_eSurf,n_nSurf,nx_,X,Y,Z);
        weights(e-ne,:,:,:,:,:)=isoparametric_surface_plus_y_sph(gps,w,eX,eY,eZ,orientation_list);
    end
    parfor e=ne+1+n_eSurf(4):1:ne+n_eSurf(5)
        [eX,eY,eZ]=get_eXYZ_sph(e,ne,nex,ney,nexney,nx,ny,nz,n,nxny,n_eSurf,n_nSurf,nx_,X,Y,Z);
        weights(e-ne,:,:,:,:,:)=isoparametric_surface_minus_z_sph(gps,w,eX,eY,eZ,orientation_list);
    end
    parfor e=ne+1+n_eSurf(5):1:ne+n_eSurf(6)
        [eX,eY,eZ]=get_eXYZ_sph(e,ne,nex,ney,nexney,nx,ny,nz,n,nxny,n_eSurf,n_nSurf,nx_,X,Y,Z);
        weights(e-ne,:,:,:,:,:)=isoparametric_surface_plus_z_sph(gps,w,eX,eY,eZ,orientation_list);
    end
elseif parallel_computing==0
    for e=ne+1:1:ne+n_eSurf(1)
        [eX,eY,eZ]=get_eXYZ_sph(e,ne,nex,ney,nexney,nx,ny,nz,n,nxny,n_eSurf,n_nSurf,nx_,X,Y,Z);
        weights(e-ne,:,:,:,:,:)=isoparametric_surface_minus_x_sph(gps,w,eX,eY,eZ,orientation_list);
    end
    for e=ne+1+n_eSurf(1):1:ne+n_eSurf(2)
        [eX,eY,eZ]=get_eXYZ_sph(e,ne,nex,ney,nexney,nx,ny,nz,n,nxny,n_eSurf,n_nSurf,nx_,X,Y,Z);
        weights(e-ne,:,:,:,:,:)=isoparametric_surface_plus_x_sph(gps,w,eX,eY,eZ,orientation_list);
    end
    for e=ne+1+n_eSurf(2):1:ne+n_eSurf(3)
        [eX,eY,eZ]=get_eXYZ_sph(e,ne,nex,ney,nexney,nx,ny,nz,n,nxny,n_eSurf,n_nSurf,nx_,X,Y,Z);
        weights(e-ne,:,:,:,:,:)=isoparametric_surface_minus_y_sph(gps,w,eX,eY,eZ,orientation_list);
    end
    for e=ne+1+n_eSurf(3):1:ne+n_eSurf(4)
        [eX,eY,eZ]=get_eXYZ_sph(e,ne,nex,ney,nexney,nx,ny,nz,n,nxny,n_eSurf,n_nSurf,nx_,X,Y,Z);
        weights(e-ne,:,:,:,:,:)=isoparametric_surface_plus_y_sph(gps,w,eX,eY,eZ,orientation_list);
    end
    for e=ne+1+n_eSurf(4):1:ne+n_eSurf(5)
        [eX,eY,eZ]=get_eXYZ_sph(e,ne,nex,ney,nexney,nx,ny,nz,n,nxny,n_eSurf,n_nSurf,nx_,X,Y,Z);
        weights(e-ne,:,:,:,:,:)=isoparametric_surface_minus_z_sph(gps,w,eX,eY,eZ,orientation_list);
    end
    for e=ne+1+n_eSurf(5):1:ne+n_eSurf(6)
        [eX,eY,eZ]=get_eXYZ_sph(e,ne,nex,ney,nexney,nx,ny,nz,n,nxny,n_eSurf,n_nSurf,nx_,X,Y,Z);
        weights(e-ne,:,:,:,:,:)=isoparametric_surface_plus_z_sph(gps,w,eX,eY,eZ,orientation_list);
    end
end

end

function sol=isoparametric_surface_minus_x_sph(gps,w,eX,eY,eZ,orientation_list)
sol=zeros(8,8,7,3,3);
for iy=1:1:3
    for iz=1:1:3
        transformation=compute_weights_specific_3d(0,gps(iy),gps(iz),gps,w,eX,eY,eZ,orientation_list);
        for iorientation=1:1:8
            for itype=1:1:8
                for iorder=1:1:7
                    sol(iorientation,itype,iorder,iy,iz)=transformation(iorientation,itype,iorder);
                end
            end
        end
    end
end
end

function sol=isoparametric_surface_plus_x_sph(gps,w,eX,eY,eZ,orientation_list)
sol=zeros(8,8,7,3,3);
for iy=1:1:3
    for iz=1:1:3
        transformation=compute_weights_specific_3d(1,gps(iy),gps(iz),gps,w,eX,eY,eZ,orientation_list);
        for iorientation=1:1:8
            for itype=1:1:8
                for iorder=1:1:7
                    sol(iorientation,itype,iorder,iy,iz)=transformation(iorientation,itype,iorder);
                end
            end
        end
    end
end
end

function sol=isoparametric_surface_minus_y_sph(gps,w,eX,eY,eZ,orientation_list)
sol=zeros(8,8,7,3,3);
for ix=1:1:3
    for iz=1:1:3
        transformation=compute_weights_specific_3d(gps(ix),0,gps(iz),gps,w,eX,eY,eZ,orientation_list);
        for iorientation=1:1:8
            for itype=1:1:8
                for iorder=1:1:7
                    sol(iorientation,itype,iorder,ix,iz)=transformation(iorientation,itype,iorder);
                end
            end
        end
    end
end
end

function sol=isoparametric_surface_plus_y_sph(gps,w,eX,eY,eZ,orientation_list)
sol=zeros(8,8,7,3,3);
for ix=1:1:3
    for iz=1:1:3
        transformation=compute_weights_specific_3d(gps(ix),1,gps(iz),gps,w,eX,eY,eZ,orientation_list);
        for iorientation=1:1:8
            for itype=1:1:8
                for iorder=1:1:7
                    sol(iorientation,itype,iorder,ix,iz)=transformation(iorientation,itype,iorder);
                end
            end
        end
    end
end
end

function sol=isoparametric_surface_minus_z_sph(gps,w,eX,eY,eZ,orientation_list)
sol=zeros(8,8,7,3,3);
for ix=1:1:3
    for iy=1:1:3
        transformation=compute_weights_specific_3d(gps(ix),gps(iy),0,gps,w,eX,eY,eZ,orientation_list);
        for iorientation=1:1:8
            for itype=1:1:8
                for iorder=1:1:7
                    sol(iorientation,itype,iorder,ix,iy)=transformation(iorientation,itype,iorder);
                end
            end
        end
    end
end
end

function sol=isoparametric_surface_plus_z_sph(gps,w,eX,eY,eZ,orientation_list)
sol=zeros(8,8,7,3,3);
for ix=1:1:3
    for iy=1:1:3
        transformation=compute_weights_specific_3d(gps(ix),gps(iy),1,gps,w,eX,eY,eZ,orientation_list);
        for iorientation=1:1:8
            for itype=1:1:8
                for iorder=1:1:7
                    sol(iorientation,itype,iorder,ix,iy)=transformation(iorientation,itype,iorder);
                end
            end
        end
    end
end
end