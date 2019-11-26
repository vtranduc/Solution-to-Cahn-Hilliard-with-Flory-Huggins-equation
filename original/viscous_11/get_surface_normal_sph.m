function normals=get_surface_normal_sph(ne,nex,ney,nexney,nx,ny,nz,n,nxny,n_eSurf,n_nSurf,nx_,X,Y,Z,parallel_computing)
normals=zeros(n_eSurf(6),3);
if parallel_computing==1
    parfor e=1:1:n_eSurf(1)
        normals(e,:)=get_surface_elemental_normal_1_sph(e+ne,ne,nex,ney,nexney,nx,ny,nz,n,nxny,n_eSurf,n_nSurf,nx_,X,Y,Z);
    end
    parfor e=n_eSurf(1)+1:1:n_eSurf(2)
        normals(e,:)=get_surface_elemental_normal_2_sph(e+ne,ne,nex,ney,nexney,nx,ny,nz,n,nxny,n_eSurf,n_nSurf,nx_,X,Y,Z);
    end
    parfor e=n_eSurf(2)+1:1:n_eSurf(3)
        normals(e,:)=get_surface_elemental_normal_3_sph(e+ne,ne,nex,ney,nexney,nx,ny,nz,n,nxny,n_eSurf,n_nSurf,nx_,X,Y,Z);
    end
    parfor e=n_eSurf(3)+1:1:n_eSurf(4)
        normals(e,:)=get_surface_elemental_normal_4_sph(e+ne,ne,nex,ney,nexney,nx,ny,nz,n,nxny,n_eSurf,n_nSurf,nx_,X,Y,Z);
    end
    parfor e=n_eSurf(4)+1:1:n_eSurf(5)
        normals(e,:)=get_surface_elemental_normal_5_sph(e+ne,ne,nex,ney,nexney,nx,ny,nz,n,nxny,n_eSurf,n_nSurf,nx_,X,Y,Z);
    end
    parfor e=n_eSurf(5)+1:1:n_eSurf(6)
        normals(e,:)=get_surface_elemental_normal_6_sph(e+ne,ne,nex,ney,nexney,nx,ny,nz,n,nxny,n_eSurf,n_nSurf,nx_,X,Y,Z);
    end
elseif parallel_computing==0
    for e=1:1:n_eSurf(1)
        normals(e,:)=get_surface_elemental_normal_1_sph(e+ne,ne,nex,ney,nexney,nx,ny,nz,n,nxny,n_eSurf,n_nSurf,nx_,X,Y,Z);
    end
    for e=n_eSurf(1)+1:1:n_eSurf(2)
        normals(e,:)=get_surface_elemental_normal_2_sph(e+ne,ne,nex,ney,nexney,nx,ny,nz,n,nxny,n_eSurf,n_nSurf,nx_,X,Y,Z);
    end
    for e=n_eSurf(2)+1:1:n_eSurf(3)
        normals(e,:)=get_surface_elemental_normal_3_sph(e+ne,ne,nex,ney,nexney,nx,ny,nz,n,nxny,n_eSurf,n_nSurf,nx_,X,Y,Z);
    end
    for e=n_eSurf(3)+1:1:n_eSurf(4)
        normals(e,:)=get_surface_elemental_normal_4_sph(e+ne,ne,nex,ney,nexney,nx,ny,nz,n,nxny,n_eSurf,n_nSurf,nx_,X,Y,Z);
    end
    for e=n_eSurf(4)+1:1:n_eSurf(5)
        normals(e,:)=get_surface_elemental_normal_5_sph(e+ne,ne,nex,ney,nexney,nx,ny,nz,n,nxny,n_eSurf,n_nSurf,nx_,X,Y,Z);
    end
    for e=n_eSurf(5)+1:1:n_eSurf(6)
        normals(e,:)=get_surface_elemental_normal_6_sph(e+ne,ne,nex,ney,nexney,nx,ny,nz,n,nxny,n_eSurf,n_nSurf,nx_,X,Y,Z);
    end
end
end

function sol=get_surface_elemental_normal_1_sph(e,ne,nex,ney,nexney,nx,ny,nz,n,nxny,n_eSurf,n_nSurf,nx_,X,Y,Z)
nodes=get_nodes_of_element_sph(e,ne,nex,ney,nexney,nx,ny,nz,n,nxny,n_eSurf,n_nSurf,nx_);
line1=[X(nodes(5))-X(nodes(1)),Y(nodes(5))-Y(nodes(1)),Z(nodes(5))-Z(nodes(1))];
line2=[X(nodes(3))-X(nodes(1)),Y(nodes(3))-Y(nodes(1)),Z(nodes(3))-Z(nodes(1))];
normal=cross(line1,line2);
sol=convert_to_unit_vector(normal);
end

function sol=get_surface_elemental_normal_2_sph(e,ne,nex,ney,nexney,nx,ny,nz,n,nxny,n_eSurf,n_nSurf,nx_,X,Y,Z)
nodes=get_nodes_of_element_sph(e,ne,nex,ney,nexney,nx,ny,nz,n,nxny,n_eSurf,n_nSurf,nx_);
line1=[X(nodes(4))-X(nodes(2)),Y(nodes(4))-Y(nodes(2)),Z(nodes(4))-Z(nodes(2))];
line2=[X(nodes(6))-X(nodes(2)),Y(nodes(6))-Y(nodes(2)),Z(nodes(6))-Z(nodes(2))];
normal=cross(line1,line2);
sol=convert_to_unit_vector(normal);
end

function sol=get_surface_elemental_normal_3_sph(e,ne,nex,ney,nexney,nx,ny,nz,n,nxny,n_eSurf,n_nSurf,nx_,X,Y,Z)
nodes=get_nodes_of_element_sph(e,ne,nex,ney,nexney,nx,ny,nz,n,nxny,n_eSurf,n_nSurf,nx_);
line1=[X(nodes(6))-X(nodes(2)),Y(nodes(6))-Y(nodes(2)),Z(nodes(6))-Z(nodes(2))];
line2=[X(nodes(1))-X(nodes(2)),Y(nodes(1))-Y(nodes(2)),Z(nodes(1))-Z(nodes(2))];
normal=cross(line1,line2);
sol=convert_to_unit_vector(normal);
end

function sol=get_surface_elemental_normal_4_sph(e,ne,nex,ney,nexney,nx,ny,nz,n,nxny,n_eSurf,n_nSurf,nx_,X,Y,Z)
nodes=get_nodes_of_element_sph(e,ne,nex,ney,nexney,nx,ny,nz,n,nxny,n_eSurf,n_nSurf,nx_);
line1=[X(nodes(3))-X(nodes(4)),Y(nodes(3))-Y(nodes(4)),Z(nodes(3))-Z(nodes(4))];
line2=[X(nodes(8))-X(nodes(4)),Y(nodes(8))-Y(nodes(4)),Z(nodes(8))-Z(nodes(4))];
normal=cross(line1,line2);
sol=convert_to_unit_vector(normal);
end

function sol=get_surface_elemental_normal_5_sph(e,ne,nex,ney,nexney,nx,ny,nz,n,nxny,n_eSurf,n_nSurf,nx_,X,Y,Z)
nodes=get_nodes_of_element_sph(e,ne,nex,ney,nexney,nx,ny,nz,n,nxny,n_eSurf,n_nSurf,nx_);
line1=[X(nodes(1))-X(nodes(2)),Y(nodes(1))-Y(nodes(2)),Z(nodes(1))-Z(nodes(2))];
line2=[X(nodes(4))-X(nodes(2)),Y(nodes(4))-Y(nodes(2)),Z(nodes(4))-Z(nodes(2))];
normal=cross(line1,line2);
sol=convert_to_unit_vector(normal);
end

function sol=get_surface_elemental_normal_6_sph(e,ne,nex,ney,nexney,nx,ny,nz,n,nxny,n_eSurf,n_nSurf,nx_,X,Y,Z)
nodes=get_nodes_of_element_sph(e,ne,nex,ney,nexney,nx,ny,nz,n,nxny,n_eSurf,n_nSurf,nx_);
line1=[X(nodes(8))-X(nodes(6)),Y(nodes(8))-Y(nodes(6)),Z(nodes(8))-Z(nodes(6))];
line2=[X(nodes(5))-X(nodes(6)),Y(nodes(5))-Y(nodes(6)),Z(nodes(5))-Z(nodes(6))];
normal=cross(line1,line2);
sol=convert_to_unit_vector(normal);
end