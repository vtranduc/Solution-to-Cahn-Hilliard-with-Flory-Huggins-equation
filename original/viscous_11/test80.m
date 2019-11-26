function solution=test80(nnz_,irow,icol,ny,n,weights,...
    ney,dt,adjusted_diffT,dxdy,conc_,terms,wTerms)

w=[5/18 4/9 5/18];
solution=zeros(1,nnz_);
parfor index=1:1:nnz_ %parallel computation must be feasible
    [nodeID1,class1]=classify_gbf(irow(index));
    [nodeID2,class2]=classify_gbf(icol(index));
    [relationship,adjacency]=identify_relationship(nodeID1,nodeID2,ny,n);
    relation_mx=local_relationship(relationship,adjacency,class1,class2);
    sj_val=0;
    for i=1:1:4
        element=relation_mx(1,i);
        if element==0
            continue
        end
        [inx,~]=inxiny_elemental(element,ney);
        sj_val=sj_val+sj_oneElement(element,relation_mx(2,i),relation_mx(3,i),weights,...
            dt,w,diffT_extractor(adjusted_diffT,inx),conc_,...
            terms,wTerms);
    end
    solution(index)=sj_val*dxdy;
end
end

function solution=sj_oneElement(element,ilocal,jlocal,weights,...
    dt,w,diffT,conc_,...
    terms,wTerms)

solution=0;
for ix=1:1:3
    for iy=1:1:3
        solution=solution+w(ix)*w(iy)*(...
            weights(ix,iy,ilocal,1)*weights(ix,iy,jlocal,1)/dt...
            -diffT(ix)*weights(ix,iy,ilocal,1)...
            *(weights(ix,iy,jlocal,1)*terms(element,ix,iy,4)...
            +terms(element,ix,iy,5)...
            *(2*(conc_(ix,iy,2,element)*weights(ix,iy,jlocal,2)...
            +conc_(ix,iy,3,element)*weights(ix,iy,jlocal,3))...
            +weights(ix,iy,jlocal,1)*terms(element,ix,iy,3))...
            +terms(element,ix,iy,6)*wTerms(ix,iy,jlocal))...
            +wTerms(ix,iy,ilocal)*wTerms(ix,iy,jlocal)...
            );
    end
end
end

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