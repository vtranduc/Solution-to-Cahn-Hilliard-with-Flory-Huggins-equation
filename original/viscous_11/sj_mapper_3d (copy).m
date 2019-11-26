function [gbfs1,gbfs2]=sj_mapper_3d(nx,ny,nz)
nnz_=nnz_sj_3d(nx,ny,nz);
gbfs1=zeros(1,nnz_);
gbfs2=zeros(1,nnz_);
sjth=0;



for zth=1:1:nz
    for yth=1:1:ny
        for xth=1:1:nx
            node1=get_node_3d(xth,yth,zth,nx,ny);
            gbfs1_nodal=get_gbs_3d(node1);
            node_adjacent=get_node_adjacent_3d(node1,nx,ny,nz);
            for type1th=1:1:8
                for iAdjacency=1:1:27
                    if node_adjacent(iAdjacency)~=0
                        node2=node_adjacent(iAdjacency);
                        gbfs2_nodal=get_gbs_3d(node2);
                        for type2th=1:1:8
                            sjth=sjth+1;
                            gbfs1(sjth)=gbfs1_nodal(type1th);
                            gbfs2(sjth)=gbfs2_nodal(type2th);
                        end
                    end
                end
            end
        end
    end
end
end

% function edge_degree=get_edge_degree_3d(xth,yth,zth,nx,ny,nz)
% edge_degree=0;
% if xth==1 || xth==nx
%     edge_degree=edge_degree+1;
% end
% if yth==1 || yth==ny
%     edge_degree=edge_degree+1;
% end
% if zth==1 || zth==nz
%     edge_degree=edge_degree+1;
% end
% end
% 
% function nInteractions=get_nInteractions_3d(edge_degree)
% if edge_degree==0
%     nInteractions=216;
% elseif edge_degree==1
%     nInteractions=144;
% elseif edge_degree==2
%     nInteractions=96;
% elseif edge_degree==3
%     nInteractions=64;
% end
% end