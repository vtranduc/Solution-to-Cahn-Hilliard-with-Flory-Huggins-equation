function solution=test92(...
    ne,nfour,weights,nx,ny,n,co,dt,dxdy,conc_,coef_T,...
    terms,wTerms,...
    grad_T,diff,...
    terms_LS)


solution=zeros(1,nfour);
cont=get_cont(ne,nx,ny,conc_,co,weights,dt);

if grad_T==0
    
    parfor gbf=1:1:nfour
        [nodeID,class]=classify_gbf(gbf);
        [adjacency,inx,~]=elemental_adjacency(nodeID,ny,n);
        sf_val=0;
        if inx>1

            if adjacency(1)~=0
                sf_val=...
                    sf_oneElement(weights,12+class,dxdy,cont,...
                    adjacency(1),coef_T,...
                    terms,wTerms);
            end
            if adjacency(2)~=0
                sf_val=sf_val+...
                    sf_oneElement(weights,8+class,dxdy,cont,...
                    adjacency(2),coef_T,...
                    terms,wTerms);
            end
        end
        if inx<nx

            if adjacency(3)~=0
                sf_val=sf_val+...
                    sf_oneElement(weights,4+class,dxdy,cont,...
                    adjacency(3),coef_T,...
                    terms,wTerms);
            end
            if adjacency(4)~=0
                sf_val=sf_val+...
                    sf_oneElement(weights,class,dxdy,cont,...
                    adjacency(4),coef_T,...
                    terms,wTerms);
            end
        end
        solution(gbf)=sf_val;
    end
elseif grad_T==1
    terms_=terms(:,:,:,1:1:2);
    for gbf=1:1:nfour
        [nodeID,class]=classify_gbf(gbf);
        [adjacency,inx,~]=elemental_adjacency(nodeID,ny,n);
        sf_val=0;
        if inx>1

            if adjacency(1)~=0
                sf_val=...
                    sf_oneElement_grad(weights,12+class,dxdy,cont,...
                    adjacency(1),...
                    terms_,wTerms,diff,...
                    terms_LS);
            end
            if adjacency(2)~=0
                sf_val=sf_val+...
                    sf_oneElement_grad(weights,8+class,dxdy,cont,...
                    adjacency(2),...
                    terms_,wTerms,diff,...
                    terms_LS);
            end
        end
        if inx<nx
            
            if adjacency(3)~=0
                sf_val=sf_val+...
                    sf_oneElement_grad(weights,4+class,dxdy,cont,...
                    adjacency(3),...
                    terms_,wTerms,diff,...
                    terms_LS);
            end
            if adjacency(4)~=0
                sf_val=sf_val+...
                    sf_oneElement_grad(weights,class,dxdy,cont,...
                    adjacency(4),...
                    terms_,wTerms,diff,...
                    terms_LS);
            end
        end
        solution(gbf)=sf_val;
    end
end

end

function solution=sf_oneElement_grad(weights,ilocal,dxdy,cont,...
            element,terms,wTerms,diff,...
            terms_LS)
        
w=[5/18 4/9 5/18];
solution=0;
for ix=1:1:3
    for iy=1:1:3
        solution=solution+w(ix)*w(iy)*(...
            weights(ix,iy,ilocal,1)*(cont(element,ix,iy)...
            -diff...
            *terms(element,ix,iy,1))...
            +terms(element,ix,iy,2)*wTerms(ix,iy,ilocal)...
            ...
            -weights(ix,iy,ilocal,1)*terms_LS(element,ix,iy,1)...
            ...
            );
    end
end
solution=dxdy*solution;
end

function solution=sf_oneElement(weights,ilocal,dxdy,cont,...
            element,diffT,terms,wTerms)
w=[5/18 4/9 5/18];
solution=0;
for ix=1:1:3
    for iy=1:1:3
        solution=solution+w(ix)*w(iy)*(...
            weights(ix,iy,ilocal,1)*cont(element,ix,iy)...
            -diffT*weights(ix,iy,ilocal,1)...
            *(terms(element,ix,iy,1)+terms(element,ix,iy,2))...
            +terms(element,ix,iy,3)*wTerms(ix,iy,ilocal)...
            );
    end
end
solution=dxdy*solution;
end

function solution=get_cont(ne,nx,ny,conc_,co,weights,dt)
solution=zeros(ne,3,3);
ney=ny-1;
co_=zeros(1,16);
e=0;
for ix=1:1:nx-1
    for iy=1:1:ney
        gbfs=elemental_gbf(ix,iy,ny);
        for i=1:1:16
            co_(i)=co(gbfs(i));
        end
        e=e+1;
        solution(e,:,:)=conc_(:,:,1,e)-conc(co_,weights,1);
    end
end
solution=solution/dt;
end