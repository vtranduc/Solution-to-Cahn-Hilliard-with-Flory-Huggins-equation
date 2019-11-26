function solution=compute_sf(nfour,weights,nx,ny,n,c,co,n1,n2,dt,dx,dy,...
    adjusted_chi,adjusted_diffT)

solution=zeros(1,nfour);

parfor gbf=1:1:nfour %PARALLEL COMPUTING MUST BE FEASIBLE!
    [nodeID,class]=classify_gbf(gbf);
    [adjacency,inx,iny]=elemental_adjacency(nodeID,ny,n);
    sf_val=0;
    ctemp=c;
    cotemp=co;
    c_=zeros(1,16);
    co_=zeros(1,16);
    if inx>1
        [chi,diffT]=chi_diffT_extractor(adjusted_chi,adjusted_diffT,inx-1);
        if adjacency(1)~=0
            gbfs=elemental_gbf(inx-1,iny-1,ny);
            for i=1:1:16
                c_(i)=ctemp(gbfs(i));
                co_(i)=cotemp(gbfs(i));
            end
            ilocal=12+class;
            sf_val=sf_val+...
                sf_oneElement(c_,co_,weights,ilocal,n1,n2,dx,dy,dt,chi,diffT);
        end
        if adjacency(2)~=0
            gbfs=elemental_gbf(inx-1,iny,ny);
            for i=1:1:16
                c_(i)=ctemp(gbfs(i));
                co_(i)=cotemp(gbfs(i));
            end
            ilocal=8+class;
            sf_val=sf_val+...
                sf_oneElement(c_,co_,weights,ilocal,n1,n2,dx,dy,dt,chi,diffT);
        end
    end
    if inx<nx
        [chi,diffT]=chi_diffT_extractor(adjusted_chi,adjusted_diffT,inx);
        if adjacency(3)~=0
            gbfs=elemental_gbf(inx,iny-1,ny);
            for i=1:1:16
                c_(i)=ctemp(gbfs(i));
                co_(i)=cotemp(gbfs(i));
            end
            ilocal=4+class;
            sf_val=sf_val+...
                sf_oneElement(c_,co_,weights,ilocal,n1,n2,dx,dy,dt,chi,diffT);
        end
        if adjacency(4)~=0
            gbfs=elemental_gbf(inx,iny,ny);
            for i=1:1:16
                c_(i)=ctemp(gbfs(i));
                co_(i)=cotemp(gbfs(i));
            end
            ilocal=class;
            sf_val=sf_val+...
                sf_oneElement(c_,co_,weights,ilocal,n1,n2,dx,dy,dt,chi,diffT);
        end
    end
    solution(gbf)=sf_val;
end

end