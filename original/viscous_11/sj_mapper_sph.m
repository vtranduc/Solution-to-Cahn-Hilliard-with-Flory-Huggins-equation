function [gbfs1,gbfs2]=sj_mapper_sph(n,nx,ny,nz,nnz_,nxny,n_nSurf,nx_,nTotal)

gbfs1=zeros(1,nnz_);
gbfs2=zeros(1,nnz_);

sjth=0;
gbf1=0;

for node=1:1:nTotal
    
    [nodes_adjacent,extruding]=get_node_adjacent_sph(node,n,nx,ny,nz,nxny,n_nSurf,nx_,nTotal);

    for itype1=1:1:8 
        gbf1=gbf1+1;
        
        %-----------------------------------------
        
        if ~isnan(nodes_adjacent)
            for iAdjacency=1:1:27
                if nodes_adjacent(iAdjacency)~=0
                    gbfs2_nodal=get_gbs_3d(nodes_adjacent(iAdjacency));
                    for itype2=1:1:8
                        sjth=sjth+1;
                        gbfs1(sjth)=gbf1;
                        gbfs2(sjth)=gbfs2_nodal(itype2);
                    end
                end
            end
        end
            
            %-----------------------------------------
        if ~isnan(extruding)
            for iAdjacency=1:1:length(extruding)
                gbfs2_nodal=get_gbs_3d(extruding(iAdjacency));
                for itype2=1:1:8
                    sjth=sjth+1;
                    gbfs1(sjth)=gbf1;
                    gbfs2(sjth)=gbfs2_nodal(itype2);
                end
            end
            %-----------------------------------------
        end
    end
end
end