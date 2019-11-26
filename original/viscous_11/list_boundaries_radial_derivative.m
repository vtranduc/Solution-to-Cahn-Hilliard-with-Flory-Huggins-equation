function [bc_gbfs_r,n_bc_nodes]=list_boundaries_radial_derivative(nx,ny,nz)

n_bc_nodes=2*((nx-2)*(ny-2)+(nx-2)*(nz-2)+(ny-2)*(nz-2))+4*(nx+ny+nz)-16;
bc_gbfs_r=zeros(1,n_bc_nodes);

index=0;
for yth=1:1:ny
    for xth=1:1:nx
        index=index+1;
        bc_gbfs_r(index)=(index-1)*4+2;
    end
end

for zth=2:1:nz-1
    pre_node=get_node_3d(nx,ny,zth-1,nx,ny)-2;
    for xth=1:1:nx
        pre_node=pre_node+1;
        index=index+1;
        bc_gbfs_r(index)=pre_node*4+2;
    end
    for yth=2:1:ny-1
        index=index+1;
        bc_gbfs_r(index)=(get_node_3d(1,yth,zth,nx,ny)-1)*4+2;
        index=index+1;
        bc_gbfs_r(index)=(get_node_3d(nx,yth,zth,nx,ny)-1)*4+2;
    end
    pre_node=get_node_3d(1,ny,zth,nx,ny)-2;
    for xth=1:1:nx
        pre_node=pre_node+1;
        index=index+1;
        bc_gbfs_r(index)=pre_node*4+2;
    end
end

pre_node=get_node_3d(1,1,nz,nx,ny)-2;

for yth=1:1:ny
    for xth=1:1:nx
        pre_node=pre_node+1;
        index=index+1;
        bc_gbfs_r(index)=pre_node*4+2;
    end
end

end