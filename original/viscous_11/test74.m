function sol=test74(bc_gbfs_r,n_bc_nodes,X,Y,Z)
sol=zeros(n_bc_nodes,3);
for inode=1:1:n_bc_nodes
    [node,~]=analyze_gbs_3d(bc_gbfs_r(inode));
    spherical_coord=cartesian_to_spherical([X(node),Y(node),Z(node)]);
    sol(inode,:)=get_normal(spherical_coord);
end

end

function sol=get_normal(spherical_coord)
sol=zeros(1,3);
sol(1)=cos(spherical_coord(2))*sin(spherical_coord(3));
sol(2)=sin(spherical_coord(2))*sin(spherical_coord(3));
sol(3)=cos(spherical_coord(3));
end