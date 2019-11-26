function sol=test72(bc_gbfs_r,n_bc_nodes,X,Y,Z)
sol=zeros(n_bc_nodes,3,3);
for inode=1:1:n_bc_nodes
    [node,~]=analyze_gbs_3d(bc_gbfs_r(inode));
    spherical_coord=cartesian_to_spherical([X(node),Y(node),Z(node)]);
    sol(inode,:,:)=rotation_matrix(spherical_coord);
end

end

function sol=rotation_matrix(spherical_coord)
sol=zeros(3,3);
sol(1,1)=cos(spherical_coord(2))*sin(spherical_coord(3));
sol(1,2)=sin(spherical_coord(2))*sin(spherical_coord(3));
sol(1,3)=cos(spherical_coord(3));
sol(2,1)=-spherical_coord(1)*sin(spherical_coord(2))*sin(spherical_coord(3));
sol(2,2)=spherical_coord(1)*cos(spherical_coord(2))*sin(spherical_coord(3));
sol(3,1)=spherical_coord(1)*cos(spherical_coord(2))*cos(spherical_coord(3));
sol(3,2)=spherical_coord(1)*sin(spherical_coord(2))*cos(spherical_coord(3));
sol(3,3)=-spherical_coord(1)*sin(spherical_coord(3));
end
