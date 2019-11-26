function solution=sj_2d(nnz_,irow,icol,ny,n,weights,...
    dt,coef_T,dxdy,terms,wTerms,...
    communication_reliever,nworkers,nCoreTasks,index_array,...
    grad_T,...
    thermophoresis,wetting,sj_assist)

w=[5/18 4/9 5/18];
solution=zeros(1,nnz_);


if wetting==0
    if grad_T==0
        if nworkers>1
            if communication_reliever==1
                sol_array=zeros(nworkers,nCoreTasks);
                parfor core=1:1:nworkers
                    sol_array(core,:)=task_splitter(index_array(core,:),...
                        nCoreTasks,irow,icol,ny,n,weights,...
                        dt,coef_T,dxdy,terms,wTerms,w);
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
                            dt,w,coef_T,...
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
                        dt,w,coef_T,...
                        terms,wTerms);
                end
                solution(index)=sj_val*dxdy;
            end
        end
    elseif grad_T==1
        if thermophoresis==0
            if nworkers>1
                if communication_reliever==1
                    sol_array=zeros(nworkers,nCoreTasks);
                    parfor core=1:1:nworkers
                        sol_array(core,:)=task_splitter_grad(index_array(core,:),...
                            nCoreTasks,irow,icol,ny,n,weights,...
                            dt,dxdy,terms,wTerms,w);
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

                            sj_val=sj_val+sj_oneElement_grad(element,relation_mx(2,i),relation_mx(3,i),weights,...
                                dt,w,...
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
        %                 [inx,~]=inxiny_elemental(element,ney);
                        sj_val=sj_val+sj_oneElement_grad(element,relation_mx(2,i),relation_mx(3,i),weights,...
                            dt,w,...
                            terms,wTerms);
                    end
                    solution(index)=sj_val*dxdy;
                end
            end
        elseif thermophoresis~=0
            if nworkers>1
                if communication_reliever==1
                    sol_array=zeros(nworkers,nCoreTasks);
                    parfor core=1:1:nworkers
                        sol_array(core,:)=task_splitter_grad_thermophoresis(...
                            index_array(core,:),...
                            nCoreTasks,irow,icol,ny,n,weights,...
                            dt,dxdy,terms,wTerms,w);
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

                            sj_val=sj_val+sj_oneElement_grad_thermophoresis(...
                                element,relation_mx(2,i),relation_mx(3,i),weights,...
                                dt,w,...
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
        %                 [inx,~]=inxiny_elemental(element,ney);
                        sj_val=sj_val+sj_oneElement_grad_thermophoresis(...
                            element,relation_mx(2,i),relation_mx(3,i),weights,...
                            dt,w,...
                            terms,wTerms);
                    end
                    solution(index)=sj_val*dxdy;
                end
            end
        end
    end
elseif wetting==1
    
    %=============CONSTRUCTION========================
    %=============CONSTRUCTION========================
    %=============CONSTRUCTION========================
    %=============CONSTRUCTION========================
    %=============CONSTRUCTION========================
    %=============CONSTRUCTION========================
    %=============CONSTRUCTION========================
    %=============CONSTRUCTION========================
    %=============CONSTRUCTION========================
    
    if nworkers>1
        if communication_reliever==1
            sol_array=zeros(nworkers,nCoreTasks);
            parfor core=1:1:nworkers
                sol_array(core,:)=task_splitter_wetting(index_array(core,:),...
                    nCoreTasks,irow,icol,ny,n,weights,...
                    dt,dxdy,terms,wTerms,w,sj_assist);
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
                    sj_val=sj_val+sj_oneElement_wetting(...
                        element,relation_mx(2,i),relation_mx(3,i),weights,...
                        dt,w,...
                        terms,wTerms,sj_assist);
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
                sj_val=sj_val+sj_oneElement_wetting(element,relation_mx(2,i),relation_mx(3,i),weights,...
                    dt,w,...
                    terms,wTerms,sj_assist);
            end
            solution(index)=sj_val*dxdy;
        end
    end
    
    %=============CONSTRUCTION========================
    %=============CONSTRUCTION========================
    %=============CONSTRUCTION========================
    %=============CONSTRUCTION========================
    %=============CONSTRUCTION========================
    %=============CONSTRUCTION========================
    %=============CONSTRUCTION========================
    %=============CONSTRUCTION========================
    %=============CONSTRUCTION========================
        
    
end

end

function sol=task_splitter_wetting(index_ends,nCoreTasks,irow,icol,ny,n,weights,...
    dt,dxdy,terms,wTerms,w,sj_assist)
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
        sj_val=sj_val+sj_oneElement_wetting(...
            element,relation_mx(2,i),relation_mx(3,i),weights,...
            dt,w,...
            terms,wTerms,sj_assist);
    end
    col_index=col_index+1;
    sol(col_index)=sj_val*dxdy;
end
end

function solution=sj_oneElement_wetting(e,ilocal,jlocal,weights,...
    dt,w,...
    terms,wTerms,sj_assist)

solution=0;
for ix=1:1:3
    for iy=1:1:3
        
        solution=solution+w(ix)*w(iy)*(...
            ...
            ...
            weights(ix,iy,ilocal,1)*weights(ix,iy,jlocal,1)/dt...
            ...
            ...
            +weights(ix,iy,jlocal,1)*(...
            terms(e,ix,iy,1)*weights(ix,iy,ilocal,2)...
            +terms(e,ix,iy,2)*weights(ix,iy,ilocal,3))...
            ...
            ...
            +terms(e,ix,iy,3)*sj_assist(ix,iy,ilocal,jlocal)...
            ...
            ...
            +wTerms(ix,iy,ilocal)*wTerms(ix,iy,jlocal)...
            ...
            ...
            );
    end
end
end

%=============CONSTRUCTION========================
    %=============CONSTRUCTION========================
    %=============CONSTRUCTION========================
    %=============CONSTRUCTION========================
    %=============CONSTRUCTION========================
    %=============CONSTRUCTION========================
    %=============CONSTRUCTION========================
    %=============CONSTRUCTION========================
    %=============CONSTRUCTION========================
        

function sol=task_splitter(index_ends,nCoreTasks,irow,icol,ny,n,weights,...
    dt,coef_T,dxdy,terms,wTerms,w)
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
            dt,w,coef_T,...
            terms,wTerms);
    end
    col_index=col_index+1;
    sol(col_index)=sj_val*dxdy;
end
end

function sol=task_splitter_grad(index_ends,nCoreTasks,irow,icol,ny,n,weights,...
    dt,dxdy,terms,wTerms,w)
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
        
%         [inx,~]=inxiny_elemental(element,ney);
        sj_val=sj_val+sj_oneElement_grad(element,relation_mx(2,i),relation_mx(3,i),weights,...
            dt,w,...
            terms,wTerms);
    end
    col_index=col_index+1;
    
    sol(col_index)=sj_val*dxdy;
    
end
end

function sol=task_splitter_grad_thermophoresis(index_ends,nCoreTasks,irow,icol,ny,n,weights,...
    dt,dxdy,terms,wTerms,w)
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
        sj_val=sj_val+sj_oneElement_grad_thermophoresis(...
            element,relation_mx(2,i),relation_mx(3,i),weights,...
            dt,w,...
            terms,wTerms);
    end
    col_index=col_index+1;
    sol(col_index)=sj_val*dxdy;
end
end

function solution=sj_oneElement_grad_thermophoresis(...
    e,ilocal,jlocal,weights,...
    dt,w,...
    terms,wTerms)

solution=0;

for ix=1:1:3
    for iy=1:1:3
        
        solution=solution+w(ix)*w(iy)*(weights(ix,iy,ilocal,1)*(...%second angled
            weights(ix,iy,jlocal,1)/dt...
            ...
            -terms(e,ix,iy,5)*weights(ix,iy,jlocal,2)...
            -terms(e,ix,iy,6)*weights(ix,iy,jlocal,3)...
            ...
            ...
            -terms(e,ix,iy,1)*wTerms(ix,iy,jlocal)...
            ...
            -terms(e,ix,iy,2)*(...
            weights(ix,iy,jlocal,1)*terms(e,ix,iy,3)...
            +terms(e,ix,iy,7)*weights(ix,iy,jlocal,2)...
            +terms(e,ix,iy,8)*weights(ix,iy,jlocal,3))...
            -weights(ix,iy,jlocal,1)*terms(e,ix,iy,4))...
            ...
            +weights(ix,iy,jlocal,1)*(...
            terms(e,ix,iy,9)*weights(ix,iy,ilocal,2)...
            +terms(e,ix,iy,10)*weights(ix,iy,ilocal,3))...
            ...
            +wTerms(ix,iy,ilocal)*wTerms(ix,iy,jlocal));
    end
end

end

function solution=sj_oneElement_grad(e,ilocal,jlocal,weights,...
    dt,w,...
    terms,wTerms)

solution=0;

for ix=1:1:3
    for iy=1:1:3
        
        solution=solution+w(ix)*w(iy)*(weights(ix,iy,ilocal,1)*(...%second angled
            weights(ix,iy,jlocal,1)/dt...
            ...
            -terms(e,ix,iy,5)*weights(ix,iy,jlocal,2)...
            -terms(e,ix,iy,6)*weights(ix,iy,jlocal,3)...
            ...
            ...
            -terms(e,ix,iy,1)*wTerms(ix,iy,jlocal)...
            ...
            -terms(e,ix,iy,2)*(...
            weights(ix,iy,jlocal,1)*terms(e,ix,iy,3)...
            +terms(e,ix,iy,7)*weights(ix,iy,jlocal,2)...
            +terms(e,ix,iy,8)*weights(ix,iy,jlocal,3))...
            -weights(ix,iy,jlocal,1)*terms(e,ix,iy,4))...
            ...
            +weights(ix,iy,jlocal,1)*(...
            terms(e,ix,iy,5)*weights(ix,iy,ilocal,2)...
            +terms(e,ix,iy,6)*weights(ix,iy,ilocal,3))...
            ...
            +wTerms(ix,iy,ilocal)*wTerms(ix,iy,jlocal));
        
%         solution=solution+w(ix)*w(iy)*(weights(ix,iy,ilocal,1)*(...%second angled
%             weights(ix,iy,jlocal,1)/dt-diff*(...*first angled
%             ...
%             terms(e,ix,iy,1)*sj_assist(e,jlocal,ix,iy,1)... %SEEMS SUSPICIOUS============
%             ...
%             +terms(e,ix,iy,2)*(...%curved
%             weights(ix,iy,jlocal,1)*terms(e,ix,iy,3)...
%             +terms(e,ix,iy,7)*weights(ix,iy,jlocal,2)...
%             +terms(e,ix,iy,8)*weights(ix,iy,jlocal,3))...
%             +weights(ix,iy,jlocal,1)*terms(e,ix,iy,4)... %MISSING TWO======================
%             +sj_assist(e,jlocal,ix,iy,2)))... %BAD=================================
%             ...
%             +weights(ix,iy,jlocal,1)*(...
%             terms(e,ix,iy,5)*weights(ix,iy,ilocal,2)...
%             +terms(e,ix,iy,6)*weights(ix,iy,ilocal,3))...
%             ...
%             +wTerms(ix,iy,ilocal)*wTerms(ix,iy,jlocal));
        
        
%         if e==56 && ix==2 && iy==3 && ilocal==2 && jlocal==3
%             warning('I have')
%             
%             weights(ix,iy,jlocal,1)/dt
%             terms(e,ix,iy,1)*sj_assist(e,jlocal,ix,iy,1)
%             sj_assist(e,jlocal,ix,iy,1)
%             
%             
%             (...%second angled
%             weights(ix,iy,jlocal,1)/dt-diff*(...*first angled
%             ...
%             terms(e,ix,iy,1)*sj_assist(e,jlocal,ix,iy,1)...
%             ...
%             +terms(e,ix,iy,2)*(...%curved
%             weights(ix,iy,jlocal,1)*terms(e,ix,iy,3)...
%             +terms(e,ix,iy,7)*weights(ix,iy,jlocal,2)...
%             +terms(e,ix,iy,8)*weights(ix,iy,jlocal,3))...
%             +weights(ix,iy,jlocal,1)*terms(e,ix,iy,4)...
%             +sj_assist(e,jlocal,ix,iy,2)));
%         end
        
%         dot_Tj_2=2*(coef_T(e,ix,iy,2)*weights(ix,iy,jlocal,2)...
%             +coef_T(e,ix,iy,3)*weights(ix,iy,jlocal,3));
%         
%         solution=solution+w(ix)*w(iy)*(weights(ix,iy,ilocal,1)*(...%second angled
%             weights(ix,iy,jlocal,1)/dt-diff*(...*first angled
%             ...
%             terms(e,ix,iy,1)*(...%squared
%             weights(ix,iy,jlocal,1)*coef_T(e,ix,iy,4)+dot_Tj_2...
%             +coef_T(e,ix,iy,1)*wTerms(ix,iy,jlocal))...
%             ...
%             +terms(e,ix,iy,2)*(weights(ix,iy,jlocal,1)*terms(e,ix,iy,3)...
%             +2*coef_T(e,ix,iy,1)*(...
%             terms(e,ix,iy,5)*weights(ix,iy,jlocal,2)...
%             +terms(e,ix,iy,6)*weights(ix,iy,jlocal,3)))...
%             ...
%             +weights(ix,iy,jlocal,1)*terms(e,ix,iy,4)...
%             ...
%             -term_persistent(e,ix,iy)*(weights(ix,iy,jlocal,1)*coef_T(e,ix,iy,4)+dot_Tj_2)...
%             ))...
%             ...
%             +wTerms(ix,iy,ilocal)*wTerms(ix,iy,jlocal)...
%             ...
%             ...
%             );
        
        
%         test=test+w(ix)*w(iy)*(weights(ix,iy,ilocal,1)*(...%second angled
%             weights(ix,iy,jlocal,1)/dt-diff*(...*first angled
%             ...
%             terms(e,ix,iy,1)*(...%squared
%             weights(ix,iy,jlocal,1)*coef_T(e,ix,iy,4)+dot_Tj_2...
%             +coef_T(e,ix,iy,1)*wTerms(ix,iy,jlocal))...
%             ...
%             +terms(e,ix,iy,2)*(weights(ix,iy,jlocal,1)*terms(e,ix,iy,3)...
%             +2*coef_T(e,ix,iy,1)*(...
%             conc_(ix,iy,2,e)*weights(ix,iy,jlocal,2)...
%             +conc_(ix,iy,3,e)*weights(ix,iy,jlocal,3)))...
%             ...
%             +weights(ix,iy,jlocal,1)*terms(e,ix,iy,4)...
%             ...
%             -term_persistent(e,ix,iy)*(weights(ix,iy,jlocal,1)*coef_T(e,ix,iy,4)+dot_Tj_2)...
%             ))...
%             ...
%             +wTerms(ix,iy,ilocal)*wTerms(ix,iy,jlocal)...
%             ...
%             ...
%             );

        
%         solution=solution+w(ix)*w(iy)*(...%second angled
%             ...
%             weights(ix,iy,ilocal,1)*weights(ix,iy,jlocal,1)/dt...
%             ...
%             -diff*weights(ix,iy,ilocal,1)...
%             *(...%first angled
%             2*coef_T(inx,ix)...
%             *(...%curly
%             weights(ix,iy,jlocal,1)*terms(element,ix,iy,3)-term_persistent(inx,ix)*weights(ix,iy,jlocal,2)...
%             )...%curly
%             +terms(element,ix,iy,4)*(2*coef_T(inx,ix)*(conc_(ix,iy,2,element)*weights(ix,iy,jlocal,2)+conc_(ix,iy,3,element)*weights(ix,iy,jlocal,3))...
%             +weights(ix,iy,jlocal,1)*terms(element,ix,iy,5))...
%             +terms(element,ix,iy,6)...
%             *(coef_T(inx,ix)*wTerms(ix,iy,jlocal)+two_slope*weights(ix,iy,jlocal,2))...
%             )...%first angled
%             +wTerms(ix,iy,ilocal)*wTerms(ix,iy,jlocal)...
%             ...
%             -weights(ix,iy,ilocal,1)*(terms_LS(element,ix,iy,2)*weights(ix,iy,jlocal,1)...
%             +terms_LS(element,ix,iy,3)*weights(ix,iy,jlocal,2))...
%             ...
%             );
    end
end

end



function solution=sj_oneElement(e,ilocal,jlocal,weights,...
    dt,w,diffT,...
    terms,wTerms)

solution=0;
for ix=1:1:3
    for iy=1:1:3

        
        solution=solution+w(ix)*w(iy)*(...
            ...
            (weights(ix,iy,ilocal,1)*(...%angled
            ...
            weights(ix,iy,jlocal,1)/dt-diffT*(...%curl
            ...
            terms(e,ix,iy,1)*wTerms(ix,iy,jlocal)...
            ...
            +terms(e,ix,iy,2)*(...%squared
            weights(ix,iy,jlocal,1)*terms(e,ix,iy,3)...
            +2*(terms(e,ix,iy,5)*weights(ix,iy,jlocal,2)...
            +terms(e,ix,iy,6)*weights(ix,iy,jlocal,3)))...
            ...
            +weights(ix,iy,jlocal,1)*terms(e,ix,iy,4))))...
            ...
            +wTerms(ix,iy,ilocal)*wTerms(ix,iy,jlocal)...
            ...
            );
        
        
%         solution=solution+w(ix)*w(iy)*(...
%             weights(ix,iy,ilocal,1)*weights(ix,iy,jlocal,1)/dt...
%             -diffT*weights(ix,iy,ilocal,1)...
%             *(weights(ix,iy,jlocal,1)*terms(element,ix,iy,4)...
%             +terms(element,ix,iy,5)...
%             *(2*(conc_(ix,iy,2,element)*weights(ix,iy,jlocal,2)...
%             +conc_(ix,iy,3,element)*weights(ix,iy,jlocal,3))...
%             +weights(ix,iy,jlocal,1)*terms(element,ix,iy,3))...
%             +terms(element,ix,iy,6)*wTerms(ix,iy,jlocal))...
%             +wTerms(ix,iy,ilocal)*wTerms(ix,iy,jlocal)...
%             );
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