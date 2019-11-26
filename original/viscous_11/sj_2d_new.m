function solution=sj_2d_new(nnz_,irow,icol,ny,n,weights,...
    dt,dxdy,terms,wwTerms,...
    nworkers,...
    ...
    viscosity)

w=[5/18 4/9 5/18];
solution=zeros(1,nnz_);

[allocation,extra,load_less,load_more]=wise_task_splitter_v2(length(irow),nworkers);

temp=zeros(nworkers,load_more);

parfor worker=1:1:nworkers
    if worker<=extra
        temp(worker,:)=sj_2d_new_assist(allocation(worker,:),...
            load_more,load_more,irow,icol,ny,n,weights,...
            dt,dxdy,terms,wwTerms,w,viscosity);
    elseif worker>extra
        temp(worker,:)=sj_2d_new_assist(allocation(worker,:),...
            load_less,load_more,irow,icol,ny,n,weights,...
            dt,dxdy,terms,wwTerms,w,viscosity);
    end
end

index=0;
for worker=1:1:nworkers
    if worker<=extra
        load=load_more;
    elseif worker>extra
        load=load_less;
    end
    for i=1:1:load
        index=index+1;
        solution(index)=temp(worker,i);
    end
end

end


%==========================================================================

function sol=sj_2d_new_assist(keys,load,load_more,irow,icol,ny,n,weights,...
    dt,dxdy,terms,wwTerms,w,viscosity)

sol=zeros(1,load_more);

for ikey=1:1:load
    [nodeID1,class1]=classify_gbf(irow(keys(ikey)));
    [nodeID2,class2]=classify_gbf(icol(keys(ikey)));
    [relationship,adjacency]=identify_relationship(nodeID1,nodeID2,ny,n);
    relation_mx=local_relationship(relationship,adjacency,class1,class2);
    sj_val=0;
    for i=1:1:4
        element=relation_mx(1,i);
        if element==0
            continue
        end
        sj_val=sj_val+sj_oneElement(element,relation_mx(2,i),relation_mx(3,i),weights,...
            dt,w,...
            terms,wwTerms,viscosity);
    end
    sol(ikey)=sj_val*dxdy;
end

end

function solution=sj_oneElement(e,ilocal,jlocal,weights,...
    dt,w,...
    terms,wwTerms,viscosity)

solution=0;

for ix=1:1:3
    for iy=1:1:3
        
        solution=solution+w(ix)*w(iy)*(...
            ...
            ...
            (weights(ix,iy,ilocal,1)*weights(ix,iy,jlocal,1)...
            ...
            ...
            +viscosity*wwTerms(ix,iy,ilocal,jlocal,1))/dt...
            ...
            ...
            +weights(ix,iy,jlocal,1)...
            *(weights(ix,iy,ilocal,2)*terms(e,ix,iy,1)...
            +weights(ix,iy,ilocal,3)*terms(e,ix,iy,2))...
            ...
            ...
            +terms(e,ix,iy,3)*wwTerms(ix,iy,ilocal,jlocal,1)...
            ...
            ...
            +wwTerms(ix,iy,ilocal,jlocal,2)...
            ...
            ...
            ...%------------------------------------------------------
            ...%------------------------------------------------------
            ...
            ...
            );
    end
end

end

%==========================================================================
%==========================================================================
%==========================================================================
%==========================================================================
%==========================================================================
%==========================================================================
%==========================================================================
%==========================================================================
%==========================================================================
%==========================================================================

function solution=local_relationship(relationship,adjacency,class1,class2)
solution=zeros(3,4);
if relationship==1
    if adjacency(1)~=0
        solution(1,1)=adjacency(1);
        solution(2,1)=12+class1;
        solution(3,1)=class2;
    end
elseif relationship==2
    if adjacency(1)~=0
        solution(1,1)=adjacency(1);
        solution(2,1)=12+class1;
        solution(3,1)=4+class2;
    end
    if adjacency(2)~=0
        solution(1,2)=adjacency(2);
        solution(2,2)=8+class1;
        solution(3,2)=class2;
    end
elseif relationship==3
    if adjacency(2)~=0
        solution(1,2)=adjacency(2);
        solution(2,2)=8+class1;
        solution(3,2)=4+class2;
    end
elseif relationship==4
    if adjacency(1)~=0
        solution(1,1)=adjacency(1);
        solution(2,1)=12+class1;
        solution(3,1)=8+class2;
    end
    if adjacency(3)~=0
        solution(1,3)=adjacency(3);
        solution(2,3)=4+class1;
        solution(3,3)=class2;
    end
elseif relationship==5
    solution(1,:)=adjacency;
    if adjacency(1)~=0
        solution(2,1)=12+class1;
        solution(3,1)=12+class2;
    end
    if adjacency(2)~=0
        solution(2,2)=8+class1;
        solution(3,2)=8+class2;
    end
    if adjacency(3)~=0
        solution(2,3)=4+class1;
        solution(3,3)=4+class2;
    end
    if adjacency(4)~=0
        solution(2,4)=class1;
        solution(3,4)=class2;
    end
elseif relationship==6
    if adjacency(2)~=0
        solution(1,2)=adjacency(2);
        solution(2,2)=8+class1;
        solution(3,2)=12+class2;
    end
    if adjacency(4)~=0
        solution(1,4)=adjacency(4);
        solution(2,4)=class1;
        solution(3,4)=4+class2;
    end
elseif relationship==7
    if adjacency(3)~=0
        solution(1,3)=adjacency(3);
        solution(2,3)=4+class1;
        solution(3,3)=8+class2;
    end
elseif relationship==8
    if adjacency(3)~=0
        solution(1,3)=adjacency(3);
        solution(2,3)=4+class1;
        solution(3,3)=12+class2;
    end
    if adjacency(4)~=0
        solution(1,4)=adjacency(4);
        solution(2,4)=class1;
        solution(3,4)=8+class2;
    end
elseif relationship==9
    if adjacency(4)~=0
        solution(1,4)=adjacency(4);
        solution(2,4)=class1;
        solution(3,4)=12+class2;
    end
end

end

function [relationship,adjacency]=identify_relationship(nodeID1,nodeID2,ny,n)
[adjacency,inx1,iny1]=elemental_adjacency(nodeID1,ny,n);
[inx2,iny2]=inxiny(nodeID2,ny);
if inx2<inx1
    if iny2<iny1
        relationship=1;
    elseif iny2==iny1
        relationship=2;
    elseif iny2>iny1
        relationship=3;
    end
elseif inx2==inx1
    if iny2<iny1
        relationship=4;
    elseif iny2==iny1
        relationship=5;
    elseif iny2>iny1
        relationship=6;
    end
elseif inx2>inx1
    if iny2<iny1
        relationship=7;
    elseif iny2==iny1
        relationship=8;
    elseif iny2>iny1
        relationship=9;
    end
end
end

function [inx,iny]=inxiny(node,ny)
iny=mod(node,ny);
if iny==0
    iny=ny;
end
inx=(node-iny)/ny+1;
end