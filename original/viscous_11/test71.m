function cr=test71(c,bc_gbfs_r,n_bc_nodes,rotation_mx)

cr=c;

for inode=1:1:n_bc_nodes
    c_r_theta_phi=c(bc_gbfs_r(inode):1:bc_gbfs_r(inode)+2);
    for compt=1:1:3
        cr(bc_gbfs_r(inode)+compt-1)=...
            rotation_mx(inode,compt,1)*c_r_theta_phi(1)...
            +rotation_mx(inode,compt,2)*c_r_theta_phi(2)...
            +rotation_mx(inode,compt,3)*c_r_theta_phi(3);
    end
end

end