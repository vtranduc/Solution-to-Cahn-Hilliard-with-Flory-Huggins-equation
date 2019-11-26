function solution=sf_2d(...
    ne,nfour,weights,nx,ny,n,co,dt,dxdy,conc_,...
    terms,wTerms,...
    grad_T,wetting)


solution=zeros(1,nfour);
cont=get_cont(ne,nx,ny,conc_,co,weights,dt);
w=[5/18 4/9 5/18];


if wetting==0
    
    if grad_T==0

        parfor gbf=1:1:nfour
            [nodeID,class]=classify_gbf(gbf);
            [adjacency,inx,~]=elemental_adjacency(nodeID,ny,n);
            sf_val=0;
            if inx>1

                if adjacency(1)~=0
                    sf_val=...
                        sf_oneElement(w,weights,12+class,dxdy,cont,...
                        adjacency(1),...
                        terms,wTerms);
                end
                if adjacency(2)~=0
                    sf_val=sf_val+...
                        sf_oneElement(w,weights,8+class,dxdy,cont,...
                        adjacency(2),...
                        terms,wTerms);
                end
            end
            if inx<nx
                if adjacency(3)~=0
                    sf_val=sf_val+...
                        sf_oneElement(w,weights,4+class,dxdy,cont,...
                        adjacency(3),...
                        terms,wTerms);
                end
                if adjacency(4)~=0
                    sf_val=sf_val+...
                        sf_oneElement(w,weights,class,dxdy,cont,...
                        adjacency(4),...
                        terms,wTerms);
                end
            end
            solution(gbf)=sf_val;
        end

    elseif grad_T==1
        parfor gbf=1:1:nfour
            [nodeID,class]=classify_gbf(gbf);
            [adjacency,inx,~]=elemental_adjacency(nodeID,ny,n);
            sf_val=0;
            if inx>1

                if adjacency(1)~=0
                    sf_val=...
                        sf_oneElement_grad(w,weights,12+class,dxdy,cont,...
                        adjacency(1),...
                        terms,wTerms);
                end
                if adjacency(2)~=0
                    sf_val=sf_val+...
                        sf_oneElement_grad(w,weights,8+class,dxdy,cont,...
                        adjacency(2),...
                        terms,wTerms);
                end
            end
            if inx<nx

                if adjacency(3)~=0
                    sf_val=sf_val+...
                        sf_oneElement_grad(w,weights,4+class,dxdy,cont,...
                        adjacency(3),...
                        terms,wTerms);
                end
                if adjacency(4)~=0
                    sf_val=sf_val+...
                        sf_oneElement_grad(w,weights,class,dxdy,cont,...
                        adjacency(4),...
                        terms,wTerms);
                end
            end
            solution(gbf)=sf_val;
        end
    end
    
elseif wetting==1

    parfor gbf=1:1:nfour
        [nodeID,class]=classify_gbf(gbf);
        [adjacency,inx,~]=elemental_adjacency(nodeID,ny,n);
        sf_val=0;
        if inx>1

            if adjacency(1)~=0
                sf_val=...
                    sf_oneElement_wetting(w,weights,12+class,dxdy,cont,...
                    adjacency(1),...
                    terms,wTerms);
            end
            if adjacency(2)~=0
                sf_val=sf_val+...
                    sf_oneElement_wetting(w,weights,8+class,dxdy,cont,...
                    adjacency(2),...
                    terms,wTerms);
            end
        end
        if inx<nx
            if adjacency(3)~=0
                sf_val=sf_val+...
                    sf_oneElement_wetting(w,weights,4+class,dxdy,cont,...
                    adjacency(3),...
                    terms,wTerms);
            end
            if adjacency(4)~=0
                sf_val=sf_val+...
                    sf_oneElement_wetting(w,weights,class,dxdy,cont,...
                    adjacency(4),...
                    terms,wTerms);
            end
        end
        solution(gbf)=sf_val;
    end
    
end

end

function solution=sf_oneElement_wetting(w,weights,ilocal,dxdy,cont,...
            element,terms,wTerms)

solution=0;
for ix=1:1:3
    for iy=1:1:3
        
        solution=solution+w(ix)*w(iy)*(...
            ...
            weights(ix,iy,ilocal,1)*cont(element,ix,iy)...
            ...
            +terms(element,ix,iy,1)*weights(ix,iy,ilocal,2)...
            +terms(element,ix,iy,2)*weights(ix,iy,ilocal,3)...
            ...
            +terms(element,ix,iy,3)*wTerms(ix,iy,ilocal)...
            ...
            ...
            );
        
        
%         solution=solution+w(ix)*w(iy)*(...
%             ...
%             weights(ix,iy,ilocal,1)*(...
%             cont(element,ix,iy)-terms(element,ix,iy,1))...
%             ...
%             +(terms(element,ix,iy,2)*weights(ix,iy,ilocal,2)...
%             +terms(element,ix,iy,3)*weights(ix,iy,ilocal,3))...
%             ...
%             +terms(element,ix,iy,4)*wTerms(ix,iy,ilocal)...
%             ...
%             );
        
        
    end
end
solution=dxdy*solution;
end

function solution=sf_oneElement_grad(w,weights,ilocal,dxdy,cont,...
            element,terms,wTerms)

solution=0;
for ix=1:1:3
    for iy=1:1:3
        
        solution=solution+w(ix)*w(iy)*(...
            ...
            weights(ix,iy,ilocal,1)*(...
            cont(element,ix,iy)-terms(element,ix,iy,1))...
            ...
            +(terms(element,ix,iy,2)*weights(ix,iy,ilocal,2)...
            +terms(element,ix,iy,3)*weights(ix,iy,ilocal,3))...
            ...
            +terms(element,ix,iy,4)*wTerms(ix,iy,ilocal)...
            ...
            );
        
        
    end
end
solution=dxdy*solution;
end

function solution=sf_oneElement(w,weights,ilocal,dxdy,cont,...
            element,terms,wTerms)

solution=0;
for ix=1:1:3
    for iy=1:1:3
        
        solution=solution+w(ix)*w(iy)*(...
            ...
            weights(ix,iy,ilocal,1)*(...
            cont(element,ix,iy)-terms(element,ix,iy,1))...
            ...
            +terms(element,ix,iy,2)*wTerms(ix,iy,ilocal)...
            ...
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