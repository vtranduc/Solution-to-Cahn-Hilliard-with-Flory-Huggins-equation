function normals=test73(n_bc_nodes,bc_gbfs_r,X,Y,Z)

% bc_gbfs_r,n_bc_nodes

normals=zeros(n_bc_nodes,3);

for i=1:1:n_bc_nodes
    [node,~]=analyze_gbs_serendipity(bc_gbfs_r(i));
    normals(i,:)=[X(node) Y(node) Z(node)];
end

end