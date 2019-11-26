function sol=get_gbfs_of_element_sph(e,ne,nex,ney,nexney,nx,ny,nz,n,nxny,n_eSurf,n_nSurf,nx_)

sol=zeros(1,64);

nodes=get_nodes_of_element_sph(e,ne,nex,ney,nexney,nx,ny,nz,n,nxny,n_eSurf,n_nSurf,nx_);

lbf=0;

for local_node=1:1:8
    gbs_nodal=get_gbs_3d(nodes(local_node));
    for type=1:1:8
        lbf=lbf+1;
        sol(lbf)=gbs_nodal(type);
    end
end

end