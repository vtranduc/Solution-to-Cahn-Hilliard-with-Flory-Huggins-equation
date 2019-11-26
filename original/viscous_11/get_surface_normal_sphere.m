function normals=get_surface_normal_sphere(n_eSurf,ne,nex,ney,nexney,nx,ny,nz,n,nxny,n_nSurf,nx_,X,Y,Z)
normals=zeros(n_eSurf(6),3,3,3);
gps=[(-sqrt(3/5)+1)/2 0.5 (sqrt(3/5)+1)/2];
index=0;
for e=ne+1:1:ne+n_eSurf(1)
    index=index+1;
    for dim1=1:1:3
        for dim2=1:1:3
            [eX,eY,eZ]=get_eXYZ_sph(...
                e,ne,nex,ney,nexney,nx,ny,nz,n,nxny,n_eSurf,n_nSurf,nx_,X,Y,Z);
            normals(index,dim1,dim2,:)=...
                convert_to_unit_vector(map_to_global_3_compts_3d(...
                0,gps(dim1),gps(dim2),eX,eY,eZ,[0 0 0]));
        end
    end 
end
for e=ne+n_eSurf(1)+1:1:ne+n_eSurf(2)
    index=index+1;
    for dim1=1:1:3
        for dim2=1:1:3
            [eX,eY,eZ]=get_eXYZ_sph(...
                e,ne,nex,ney,nexney,nx,ny,nz,n,nxny,n_eSurf,n_nSurf,nx_,X,Y,Z);
            normals(index,dim1,dim2,:)=...
                convert_to_unit_vector(map_to_global_3_compts_3d(...
                1,gps(dim1),gps(dim2),eX,eY,eZ,[0 0 0]));
        end
    end 
end
for e=ne+n_eSurf(2)+1:1:ne+n_eSurf(3)
    index=index+1;
    for dim1=1:1:3
        for dim2=1:1:3
            [eX,eY,eZ]=get_eXYZ_sph(...
                e,ne,nex,ney,nexney,nx,ny,nz,n,nxny,n_eSurf,n_nSurf,nx_,X,Y,Z);
            normals(index,dim1,dim2,:)=...
                convert_to_unit_vector(map_to_global_3_compts_3d(...
                gps(dim1),0,gps(dim2),eX,eY,eZ,[0 0 0]));
        end
    end 
end
for e=ne+n_eSurf(3)+1:1:ne+n_eSurf(4)
    index=index+1;
    for dim1=1:1:3
        for dim2=1:1:3
            [eX,eY,eZ]=get_eXYZ_sph(...
                e,ne,nex,ney,nexney,nx,ny,nz,n,nxny,n_eSurf,n_nSurf,nx_,X,Y,Z);
            normals(index,dim1,dim2,:)=...
                convert_to_unit_vector(map_to_global_3_compts_3d(...
                gps(dim1),1,gps(dim2),eX,eY,eZ,[0 0 0]));
        end
    end 
end
for e=ne+n_eSurf(4)+1:1:ne+n_eSurf(5)
    index=index+1;
    for dim1=1:1:3
        for dim2=1:1:3
            [eX,eY,eZ]=get_eXYZ_sph(...
                e,ne,nex,ney,nexney,nx,ny,nz,n,nxny,n_eSurf,n_nSurf,nx_,X,Y,Z);
            normals(index,dim1,dim2,:)=...
                convert_to_unit_vector(map_to_global_3_compts_3d(...
                gps(dim1),gps(dim2),0,eX,eY,eZ,[0 0 0]));
        end
    end 
end
for e=ne+n_eSurf(5)+1:1:ne+n_eSurf(6)
    index=index+1;
    for dim1=1:1:3
        for dim2=1:1:3
            [eX,eY,eZ]=get_eXYZ_sph(...
                e,ne,nex,ney,nexney,nx,ny,nz,n,nxny,n_eSurf,n_nSurf,nx_,X,Y,Z);
            normals(index,dim1,dim2,:)=...
                convert_to_unit_vector(map_to_global_3_compts_3d(...
                gps(dim1),gps(dim2),1,eX,eY,eZ,[0 0 0]));
        end
    end 
end
end

% function convert