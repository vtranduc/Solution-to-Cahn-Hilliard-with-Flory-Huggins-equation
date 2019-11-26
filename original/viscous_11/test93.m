function solution=test93(nnz_,irow,icol,ny,n,weights,...
    dt,coef_T,dxdy,conc_,terms,wTerms,...
    communication_reliever,nworkers,nCoreTasks,index_array,...
    grad_T,term_persistent,two_slope,diff,ney,...
    terms_LS)
w=[5/18 4/9 5/18];
solution=zeros(1,nnz_);
if grad_T==0
    if nworkers>1
        if communication_reliever==1
            sol_array=zeros(nworkers,nCoreTasks);
            parfor core=1:1:nworkers
                sol_array(core,:)=task_splitter(index_array(core,:),...
                    nCoreTasks,irow,icol,ny,n,weights,...
                    dt,coef_T,dxdy,conc_,terms,wTerms,w);
            end
            for core=1:1:nworkers
                col_index=0;
                for index=index_array(core,1):1:index_array(core,2)
                    col_index=col_index+1;
                    solution(index)=sol_array(core,col_index);
                end
            end
        elseif communication_reliever==0
            parfor index=1:1:nnz_
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
                    sj_val=sj_val+sj_oneElement(element,relation_mx(2,i),relation_mx(3,i),weights,...
                        dt,w,coef_T,conc_,...
                        terms,wTerms);
                end
                solution(index)=sj_val*dxdy;
            end
        end
    elseif nworkers==1
        for index=1:1:nnz_
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
                sj_val=sj_val+sj_oneElement(element,relation_mx(2,i),relation_mx(3,i),weights,...
                    dt,w,coef_T,conc_,...
                    terms,wTerms);
            end
            solution(index)=sj_val*dxdy;
        end
    end
elseif grad_T==1
    if nworkers>1
        if communication_reliever==1
            sol_array=zeros(nworkers,nCoreTasks);
            parfor core=1:1:nworkers
                sol_array(core,:)=task_splitter_grad(index_array(core,:),...
                    nCoreTasks,irow,icol,ny,n,weights,...
                    dt,coef_T,dxdy,conc_,terms,wTerms,w,...
                    diff,term_persistent,two_slope,ney,...
                    terms_LS);
            end
            for core=1:1:nworkers
                col_index=0;
                for index=index_array(core,1):1:index_array(core,2)
                    col_index=col_index+1;
                    solution(index)=sol_array(core,col_index);
                end
            end    
        elseif communication_reliever==0
            parfor index=1:1:nnz_
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
                    sj_val=sj_val+sj_oneElement_grad(element,relation_mx(2,i),relation_mx(3,i),weights,...
                        dt,w,diff,conc_,...
                        terms,wTerms,coef_T,inx,term_persistent,two_slope);
                end
                solution(index)=sj_val*dxdy;
            end
        end
    elseif nworkers==1
        for index=1:1:nnz_
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
                sj_val=sj_val+sj_oneElement_grad(element,relation_mx(2,i),relation_mx(3,i),weights,...
                    dt,w,diff,conc_,...
                    terms,wTerms,coef_T,inx,term_persistent,two_slope);
            end
            solution(index)=sj_val*dxdy;
        end
    end
end
end

function sol=task_splitter(index_ends,nCoreTasks,irow,icol,ny,n,weights,...
    dt,coef_T,dxdy,conc_,terms,wTerms,w)
col_index=0;
sol=zeros(1,nCoreTasks);
for index=index_ends(1):1:index_ends(2)
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
        sj_val=sj_val+sj_oneElement(element,relation_mx(2,i),relation_mx(3,i),weights,...
            dt,w,coef_T,conc_,...
            terms,wTerms);
    end
    col_index=col_index+1;
    sol(col_index)=sj_val*dxdy;
end
end

function sol=task_splitter_grad(index_ends,nCoreTasks,irow,icol,ny,n,weights,...
    dt,coef_T,dxdy,conc_,terms,wTerms,w,...
    diff,term_persistent,two_slope,ney,...
    terms_LS)
col_index=0;
sol=zeros(1,nCoreTasks);
for index=index_ends(1):1:index_ends(2)
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
        sj_val=sj_val+sj_oneElement_grad(element,relation_mx(2,i),relation_mx(3,i),weights,...
            dt,w,diff,conc_,...
            terms,wTerms,coef_T,inx,term_persistent,two_slope,...
            terms_LS);
    end
    col_index=col_index+1;
    sol(col_index)=sj_val*dxdy;
end
end

function solution=sj_oneElement_grad(element,ilocal,jlocal,weights,...
    dt,w,diff,conc_,...
    terms,wTerms,coef_T,inx,term_persistent,two_slope,...
    terms_LS)

solution=0;
for ix=1:1:3
    for iy=1:1:3
        
        solution=solution+w(ix)*w(iy)*(...%second angled
            ...
            weights(ix,iy,ilocal,1)*weights(ix,iy,jlocal,1)/dt...
            ...
            -diff*weights(ix,iy,ilocal,1)...
            *(...%first angled
            2*coef_T(inx,ix)...
            *(...%curly
            weights(ix,iy,jlocal,1)*terms(element,ix,iy,3)-term_persistent(inx,ix)*weights(ix,iy,jlocal,2)...
            )...%curly
            +terms(element,ix,iy,4)*(2*coef_T(inx,ix)*(conc_(ix,iy,2,element)*weights(ix,iy,jlocal,2)+conc_(ix,iy,3,element)*weights(ix,iy,jlocal,3))...
            +weights(ix,iy,jlocal,1)*terms(element,ix,iy,5))...
            +terms(element,ix,iy,6)...
            *(coef_T(inx,ix)*wTerms(ix,iy,jlocal)+two_slope*weights(ix,iy,jlocal,2))...
            )...%first angled
            +wTerms(ix,iy,ilocal)*wTerms(ix,iy,jlocal)...
            ...
            -weights(ix,iy,ilocal,1)*(terms_LS(element,ix,iy,2)*weights(ix,iy,jlocal,1)...
            +terms_LS(element,ix,iy,3)*weights(ix,iy,jlocal,2))...
            ...
            );
    end
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
            -diffT*weights(ix,iy,ilocal,1)...
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