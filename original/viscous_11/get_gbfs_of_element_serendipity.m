function gbfs=get_gbfs_of_element_serendipity(e,nex,nexney,nxny,nx)

gbfs=zeros(1,32);

nodes=get_nodes_of_element_3d(e,nex,nexney,nxny,nx);

pre_gbf=(nodes(1)-1)*4+1;
gbfs(1:1:8)=pre_gbf:1:pre_gbf+7;

pre_gbf=(nodes(3)-1)*4+1;
gbfs(9:1:16)=pre_gbf:1:pre_gbf+7;

pre_gbf=(nodes(5)-1)*4+1;
gbfs(17:1:24)=pre_gbf:1:pre_gbf+7;

pre_gbf=(nodes(7)-1)*4+1;
gbfs(25:1:32)=pre_gbf:1:pre_gbf+7;

end