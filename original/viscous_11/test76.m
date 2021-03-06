function solution=test76(...
    nfour,weights,nx,ny,n,co,n1,n2,dt,dxdy,conc_,adjusted_chi,adjusted_diffT,...
    terms,wTerms)

solution=zeros(1,nfour);

parfor gbf=1:1:nfour %PARALLEL COMPUTING MUST BE FEASIBLE!
    [nodeID,class]=classify_gbf(gbf);
    [adjacency,inx,iny]=elemental_adjacency(nodeID,ny,n);
%     chi__=chi;chi_=chi__(inx);
%     T__=T;T_=T__(inx);
    sf_val=0;
    cotemp=co;
    co_=zeros(1,16);
    if inx>1
        [chi,diffT]=chi_diffT_extractor(adjusted_chi,adjusted_diffT,inx-1);
        if adjacency(1)~=0
            gbfs=elemental_gbf(inx-1,iny-1,ny);
            for i=1:1:16
                co_(i)=cotemp(gbfs(i));
            end
            ilocal=12+class;
            sf_val=sf_val+...
                sf_oneElement(co_,weights,ilocal,n1,n2,dxdy,dt,conc_,...
                adjacency(1),chi,diffT,...
                terms,wTerms);
        end
        if adjacency(2)~=0
            gbfs=elemental_gbf(inx-1,iny,ny);
            for i=1:1:16
                co_(i)=cotemp(gbfs(i));
            end
            ilocal=8+class;
            sf_val=sf_val+...
                sf_oneElement(co_,weights,ilocal,n1,n2,dxdy,dt,conc_,...
                adjacency(2),chi,diffT,...
                terms,wTerms);
        end
    end
    if inx<nx
        [chi,diffT]=chi_diffT_extractor(adjusted_chi,adjusted_diffT,inx);
        if adjacency(3)~=0
            gbfs=elemental_gbf(inx,iny-1,ny);
            for i=1:1:16
                co_(i)=cotemp(gbfs(i));
            end
            ilocal=4+class;
            sf_val=sf_val+...
                sf_oneElement(co_,weights,ilocal,n1,n2,dxdy,dt,conc_,...
                adjacency(3),chi,diffT,...
                terms,wTerms);
        end
        if adjacency(4)~=0
            gbfs=elemental_gbf(inx,iny,ny);
            for i=1:1:16
                co_(i)=cotemp(gbfs(i));
            end
            ilocal=class;
            sf_val=sf_val+...
                sf_oneElement(co_,weights,ilocal,n1,n2,dxdy,dt,conc_,...
                adjacency(4),chi,diffT,...
                terms,wTerms);
        end
    end
    solution(gbf)=sf_val;
end

end

function solution=sf_oneElement(co_,weights,ilocal,~,~,dxdy,dt,conc_,...
            element,~,diffT,...
            terms,wTerms)

w=[5/18 4/9 5/18];

con=conc_(:,:,1,element);
cono=conc(co_,weights,1);
% conx=conc_(:,:,2,element);
% cony=conc_(:,:,3,element);
% conxx=conc_(:,:,4,element);
% conyy=conc_(:,:,5,element);

cont=(con-cono)/dt;

%No computation is required here, so it does not have to be input
phi=weights(:,:,ilocal,1);
% phix=weights(:,:,ilocal,2)
% phiy=weights(:,:,ilocal,3)
% phixx=weights(:,:,ilocal,4);
% phiyy=weights(:,:,ilocal,5);

solution=0;
for ix=1:1:3
    for iy=1:1:3
        solution=solution+w(ix)*w(iy)*(...
            phi(ix,iy)*cont(ix,iy)-diffT(ix)*phi(ix,iy)...
            *(terms(element,ix,iy,1)+terms(element,ix,iy,2))...
            +terms(element,ix,iy,3)*wTerms(ix,iy,ilocal)...
            );
    end
end
solution=dxdy*solution;
end