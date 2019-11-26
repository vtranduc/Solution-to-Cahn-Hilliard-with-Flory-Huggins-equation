function solution=compute_sf(nfour,weights,nx,ny,n,...
    adjusted_chi,c,co,adjusted_diffT,n1,n2,dt,dx,dy)


%===== Inputs =================

%==============================
solution=zeros(1,nfour);
% display('Enter sf')

% for gbf=43:1:43
for gbf=1:1:nfour %PARALLEL COMPUTING MUST BE FEASIBLE!
    [nodeID,class]=classify_gbf(gbf);
    [adjacency,inx,iny]=elemental_adjacency(nodeID,ny,n);
%     x_coord_=x_coord;y_coord_=y_coord;
    sf_val=0;
    ctemp=c;
    cotemp=co;
    c_=zeros(1,16);
    co_=zeros(1,16);
    if inx>1
        [chi,diffT]=chi_diffT_extractor(adjusted_chi,adjusted_diffT,inx-1);
        if adjacency(1)~=0
    %         dx=x_coord_(inx)-x_coord_(inx-1);
    %         dy=y_coord_(iny)-y_coord_(iny-1);
            gbfs=elemental_gbf(inx-1,iny-1,ny);
            for i=1:1:16
                c_(i)=ctemp(gbfs(i));
                co_(i)=cotemp(gbfs(i));
            end
            ilocal=12+class;
            sf_val=sf_val+...
                sf_oneElement(chi,c_,co_,weights,ilocal,diffT,n1,n2,dx,dy,dt);
        end
        if adjacency(2)~=0
    %         dx=x_coord_(inx)-x_coord_(inx-1);
    %         dy=y_coord_(iny+1)-y_coord_(iny);
            gbfs=elemental_gbf(inx-1,iny,ny);
            for i=1:1:16
                c_(i)=ctemp(gbfs(i));
                co_(i)=cotemp(gbfs(i));
            end
            ilocal=8+class;
            sf_val=sf_val+...
                sf_oneElement(chi,c_,co_,weights,ilocal,diffT,n1,n2,dx,dy,dt);
        end
    elseif inx<nx
        [chi,diffT]=chi_diffT_extractor(adjusted_chi,adjusted_diffT,inx);
        if adjacency(3)~=0
            gbfs=elemental_gbf(inx,iny-1,ny);
            for i=1:1:16
                c_(i)=ctemp(gbfs(i));
                co_(i)=cotemp(gbfs(i));
            end
            ilocal=4+class;
            sf_val=sf_val+...
                sf_oneElement(chi,c_,co_,weights,ilocal,diffT,n1,n2,dx,dy,dt);
        end
        if adjacency(4)~=0
            gbfs=elemental_gbf(inx,iny,ny);
            for i=1:1:16
                c_(i)=ctemp(gbfs(i));
                co_(i)=cotemp(gbfs(i));
            end
            ilocal=class;
            sf_val=sf_val+...
                sf_oneElement(chi,c_,co_,weights,ilocal,diffT,n1,n2,dx,dy,dt);
        end
    end
    solution(gbf)=sf_val;
end

end