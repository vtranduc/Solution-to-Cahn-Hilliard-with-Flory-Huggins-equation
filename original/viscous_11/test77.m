function solution=test77(nnz_,irow,icol,ny,n,weights,...
    ney,dt,adjusted_chi,n1,n2,adjusted_diffT,dxdy,conc_,...
    terms,wTerms)
w=[5/18 4/9 5/18];
solution=zeros(1,nnz_);
% gp=[(-sqrt(3/5)+1)/2 0.5 (sqrt(3/5)+1)/2];
% adjusted_gp=gp*dx;
parfor index=1:1:nnz_ %parallel computation must be feasible
% for index=100:1000
    gbf1=irow(index);
    gbf2=icol(index);
    [nodeID1,class1]=classify_gbf(gbf1);
    [nodeID2,class2]=classify_gbf(gbf2);
    [relationship,adjacency]=identify_relationship(nodeID1,nodeID2,ny,n);
    
%     [gbf1,gbf2,relationship,adjacency];
    
    relation_mx=local_relationship(relationship,adjacency,class1,class2);
    
    sj_val=0;
    for i=1:1:4
        element=relation_mx(1,i);
        if element==0
            continue
        end
        ilocal=relation_mx(2,i);
        jlocal=relation_mx(3,i);
        [inx,~]=inxiny_elemental(element,ney);
        
        [chi,diffT]=chi_diffT_extractor(adjusted_chi,adjusted_diffT,inx);
        sj_val=sj_val+sj_oneElement(element,ilocal,jlocal,weights,chi,...
            w,n1,n2,dt,diffT,conc_,...
            terms,wTerms);
    end
    
    solution(index)=sj_val*dxdy;
    
end

end

function [chi,diffT]=chi_diffT_extractor(adjusted_chi,adjusted_diffT,inx)
chi=adjusted_chi(inx,:);
diffT=adjusted_diffT(inx,:);
end

function solution=sj_oneElement(element,ilocal,jlocal,weights,~,...
    w,~,~,dt,diffT,conc_,...
    terms,wTerms)
% con=conc_(:,:,1,element);
% conx=conc_(:,:,2,element);
% cony=conc_(:,:,3,element);
% conxx=conc_(:,:,4,element);
% conyy=conc_(:,:,5,element);

iphi=weights(:,:,ilocal,1);
% iphixx=weights(:,:,ilocal,4);
% iphiyy=weights(:,:,ilocal,5);

jphi=weights(:,:,jlocal,1);
jphix=weights(:,:,jlocal,2);
jphiy=weights(:,:,jlocal,3);
% jphixx=weights(:,:,jlocal,4);
% jphiyy=weights(:,:,jlocal,5);

solution=0;
for ix=1:1:3
    for iy=1:1:3
        solution=solution+w(ix)*w(iy)*(...
            iphi(ix,iy)*jphi(ix,iy)/dt...
            -diffT(ix)*iphi(ix,iy)...
            *(jphi(ix,iy)*terms(element,ix,iy,4)...
            +terms(element,ix,iy,5)...
            *(2*(conc_(ix,iy,2,element)*jphix(ix,iy)+conc_(ix,iy,3,element)*jphiy(ix,iy))...
            +jphi(ix,iy)*terms(element,ix,iy,3))...
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