function solution=sf_2d(...
    ne,nfour,weights,nx,ny,n,co,dt,dxdy,conc_,...
    terms,wTerms,viscosity,sfsj_assist)


solution=zeros(1,nfour);
cont=get_cont(ne,nx,ny,conc_,co,weights,dt);

[contx,conty]=get_cont_xy(ne,nx,ny,conc_,co,weights,dt);

contx=viscosity*contx;
conty=viscosity*conty;

w=[5/18 4/9 5/18];
    
parfor gbf=1:1:nfour
    [nodeID,class]=classify_gbf(gbf);
    [adjacency,inx,~]=elemental_adjacency(nodeID,ny,n);
    sf_val=0;
    if inx>1

        if adjacency(1)~=0
            sf_val=...
                sf_oneElement(w,weights,12+class,dxdy,cont,...
                adjacency(1),...
                terms,wTerms,contx,conty,sfsj_assist);
        end
        if adjacency(2)~=0
            sf_val=sf_val+...
                sf_oneElement(w,weights,8+class,dxdy,cont,...
                adjacency(2),...
                terms,wTerms,contx,conty,sfsj_assist);
        end
    end
    if inx<nx
        if adjacency(3)~=0
            sf_val=sf_val+...
                sf_oneElement(w,weights,4+class,dxdy,cont,...
                adjacency(3),...
                terms,wTerms,contx,conty,sfsj_assist);
        end
        if adjacency(4)~=0
            sf_val=sf_val+...
                sf_oneElement(w,weights,class,dxdy,cont,...
                adjacency(4),...
                terms,wTerms,contx,conty,sfsj_assist);
        end
    end
    solution(gbf)=sf_val;
end

end


function solution=sf_oneElement(w,weights,ilocal,dxdy,cont,...
            e,terms,wTerms,contx,conty,sfsj_assist)

solution=0;
for ix=1:1:3
    for iy=1:1:3
        
        solution=solution+w(ix)*w(iy)*(...
            ...
            cont(e,ix,iy)*sfsj_assist(e,ix,iy,ilocal,1)...
            ...
            -terms(e,ix,iy,1)*sfsj_assist(e,ix,iy,ilocal,2)...
            ...
            +terms(e,ix,iy,2)*sfsj_assist(e,ix,iy,ilocal,3)...
            ...
            );
        
%         solution=solution+w(ix)*w(iy)*(...
%             ...
%             weights(ix,iy,ilocal,1)*cont(e,ix,iy)...
%             ...
%             ...
%             +weights(ix,iy,ilocal,2)*terms(e,ix,iy,1)...
%             +weights(ix,iy,ilocal,3)*terms(e,ix,iy,2)...
%             ...
%             ...
%             +wTerms(ix,iy,ilocal)*terms(e,ix,iy,3)...
%             ...
%             ...
%             ...%------------------------------------------------
%             +weights(ix,iy,ilocal,2)*contx(e,ix,iy)...
%             +weights(ix,iy,ilocal,3)*conty(e,ix,iy)...
%             ...%------------------------------------------------
%             ...
%             ...
%             );
        
    end
end
solution=dxdy*solution;
end

function [contx,conty]=get_cont_xy(ne,nx,ny,conc_,co,weights,dt)

contx=zeros(ne,3,3);
conty=zeros(ne,3,3);
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
        contx(e,:,:)=conc_(:,:,2,e)-conc(co_,weights,2);
        conty(e,:,:)=conc_(:,:,3,e)-conc(co_,weights,3);
    end
end
contx=contx/dt;
conty=conty/dt;
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