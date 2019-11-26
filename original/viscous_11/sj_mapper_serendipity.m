function [gbfs1,gbfs2]=sj_mapper_serendipity(nx,ny,nz,nnz_)
% function [gbfs1,gbfs2,iZeros,iOnes]=sj_mapper_serendipity(nx,ny,nz)
%--------------------
% nnz_=nnz_sj_serendipity(nx,ny,nz);
%--------------------
gbfs1=zeros(1,nnz_);
gbfs2=zeros(1,nnz_);
sjth=0;

% n_edge_deg_1=2*((nx-2)*(ny-2)+(nx-2)*(nz-2)+(ny-2)*(nz-2));
% n_edge_deg_2=4*(nx+ny+nz)-24;
% n_edge_deg_3=8;

% iZeros=zeros(1,572*n_edge_deg_1+570*n_edge_deg_2+441*n_edge_deg_3);
% iOnes=zeros(1,4*n_edge_deg_1+6*n_edge_deg_2+7*n_edge_deg_3);
% 
% iZero=0;
% iOne=0;

for zth=1:1:nz
    for yth=1:1:ny
        for xth=1:1:nx
            
            %===============================
            edge_status=get_edge_status_3d(xth,yth,zth,nx,ny,nz);
            %===============================
            
            if edge_status==0
            
                node1=get_node_3d(xth,yth,zth,nx,ny);
                
                %--------------------------
%                 gbfs1_nodal=get_gbs_3d(node1);
                gbfs1_nodal=get_gbs_serendipity(node1);
                %---------------------------
                node_adjacent=get_node_adjacent_3d(node1,nx,ny,nz);
                
                
                %-----------
                for type1th=1:1:4
                    %--------------------
                    for iAdjacency=1:1:27
                        %--
                        node2=node_adjacent(iAdjacency);
                        %--------------------------
%                         gbfs2_nodal=get_gbs_3d(node2);
                        gbfs2_nodal=get_gbs_serendipity(node2);
                        %--------------------------
                        for type2th=1:1:4
                            sjth=sjth+1;
                            gbfs1(sjth)=gbfs1_nodal(type1th);
                            gbfs2(sjth)=gbfs2_nodal(type2th);
                        end
                    end
                end
            else
                %=================================================
                node1=get_node_3d(xth,yth,zth,nx,ny);
                %--------------------------
%                 gbfs1_nodal=get_gbs_3d(node1);
                gbfs1_nodal=get_gbs_serendipity(node1);
                %------------------------------
                node_adjacent=get_node_adjacent_3d(node1,nx,ny,nz);
                
                for type1th=1:1:4
                    for iAdjacency=1:1:27
                        if node_adjacent(iAdjacency)~=0
                            node2=node_adjacent(iAdjacency);
                            gbfs2_nodal=get_gbs_serendipity(node2);
                            for type2th=1:1:4
                                sjth=sjth+1;
                                gbfs1(sjth)=gbfs1_nodal(type1th);
                                gbfs2(sjth)=gbfs2_nodal(type2th);
                            end
                            
                        end
                    end
                end
                
%                 for type1th=1:1:4
%                     if edge_status(type1th)==0
%                     
%                         for iAdjacency=1:1:27
%                             if node_adjacent(iAdjacency)~=0
%                                 node2=node_adjacent(iAdjacency);
%                                 gbfs2_nodal=get_gbs_3d(node2);
%                                 for type2th=1:1:8
%                                     sjth=sjth+1;
%                                     gbfs1(sjth)=gbfs1_nodal(type1th);
%                                     gbfs2(sjth)=gbfs2_nodal(type2th);
%                                 end
%                             end
%                         end
%                         
%                     elseif edge_status(type1th)==1
%                         
%                         for iAdjacency=1:1:27
%                             if node_adjacent(iAdjacency)~=0
%                                 node2=node_adjacent(iAdjacency);
%                                 gbfs2_nodal=get_gbs_3d(node2);
%                                 for type2th=1:1:8
%                                     sjth=sjth+1;
%                                     gbfs1(sjth)=gbfs1_nodal(type1th);
%                                     gbfs2(sjth)=gbfs2_nodal(type2th);
%                                     
% %                                     if gbfs1(sjth)~=gbfs2(sjth)
% %                                         iZero=iZero+1;
% %                                         iZeros(iZero)=sjth;
% %                                     elseif gbfs1(sjth)==gbfs2(sjth)
% %                                         iOne=iOne+1;
% %                                         iOnes(iOne)=sjth;
% %                                     end
%                                 end
%                             end
%                         end
%                         
%                     end
%                 end
                %========================================================
            end
        end
    end
end
end



function status=get_edge_status_3d(xth,yth,zth,nx,ny,nz)
if xth==1 || xth==nx || yth==1 || yth==ny || zth==1 || zth==nz
    status=1;
else
    status=0;
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