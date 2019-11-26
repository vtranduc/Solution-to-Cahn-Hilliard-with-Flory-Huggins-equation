function test127

clear
clc
close all

gp=[(-sqrt(3/5)+1)/2 0.5 (sqrt(3/5)+1)/2];


coords=[0 0
    1 0
    0 1
    1 1];
mx_coef=generate_matrix_serendipity(coords);
mx_constraints=mx_constraints_assembler([0 0 0 0]);
mx_poly_coef=generate_poly_coef(mx_coef,mx_constraints);
% mx_P=zeros(3,3,12,6);
% for ialpha=1:1:3
%     for ibeta=1:1:3
%         mx_P(ialpha,ibeta,:,:)=P(gp(ialpha),gp(ibeta),mx_poly_coef);
%     end
% end

%=========================================IMPORTANT=======

% gp=zeros(1,3);

% mx_poly_coef=getting_odd_poly_coef(gp);

%=========================================================


% mx_poly_coef
% coef=Hermitian_poly_coef_serendipity(alpha,beta,der_compt)

nx=20;
ny=21;

x=linspace(0,1,nx);
y=linspace(0,1,ny);

z=zeros(ny,nx);

% ilocal=3;
der_compt=1;

aruto=zeros(12,ny,nx);
figure(1)
for ilocal=1:1:12
    for ix=1:1:nx
        for iy=1:1:ny
            variables=Hermitian_poly_coef_serendipity(x(ix),y(iy),der_compt);
            z(iy,ix)=dot(variables,mx_poly_coef(ilocal,:));
            aruto(ilocal,iy,ix)=z(iy,ix);
        end
    end
    subplot(4,3,ilocal)
    surf(z)
    xlabel('x');ylabel('y')
    
end


% return

%-----------------------------

alpha=0;
beta=0;

% der_compt=1;
nodal_weights=ones(1,12);

% nodal_weights(5)=0;
% nodal_weights(11)=0;

nodal_weights(1)=0.5;
nodal_weights(2)=1;
nodal_weights(3)=1;

nodal_weights(4)=0.5;
nodal_weights(5)=0;
nodal_weights(6)=-1;

nodal_weights(7)=0.5;
nodal_weights(8)=1;
nodal_weights(9)=1;

nodal_weights(10)=0.5;
nodal_weights(11)=0;
nodal_weights(12)=-1;

% nodal_weights(12)=51;

% sol=get_conc_test(alpha,beta,mx_poly_coef,nodal_weights,der_compt);

for ix=1:1:nx
    for iy=1:1:ny
        z(iy,ix)=get_conc_test(x(ix),y(iy),mx_poly_coef,nodal_weights,der_compt);
    end
end

figure(2)
subplot(1,2,1)
surf(x,y,z)
xlabel('x');ylabel('y')
title('serendipity')
% figure(1)

%---------------------NO SERENDIPITY

% der_compt=2;
nodal_weights=ones(1,16);

% nodal_weights(2)=0;
% nodal_weights(4)=0;
% nodal_weights(6)=0;
% nodal_weights(8)=0;
% 
nodal_weights(10)=0;
nodal_weights(14)=0;

nodal_weights(12)=0;
nodal_weights(16)=0;


der_compt=der_compt+1;
for ix=1:1:nx
    for iy=1:1:ny
        z(iy,ix)=generateWeights_test(x(ix),y(iy),der_compt,nodal_weights);
    end
end

% figure(3)
subplot(1,2,2)
surf(x,y,z)
title('normal')


figure(3)
ifig=0;
for iorientation=1:1:4
    for itype=1:1:4
        ifig=ifig+1;
        subplot(4,4,ifig)
        
        for ix=1:1:nx
            for iy=1:1:ny
                z(iy,ix)=generateWeights_test_2(x(ix),y(iy),der_compt,iorientation,itype);
            end
        end
        
        surf(x,y,z)
        xlabel('x');ylabel('y')
        
    end
end


return
%----------------------------------------------------------


coef=Hermitian_poly_coef_serendipity(alpha,beta,der_compt);


mx_var=zeros(24,12);


alphas=[0 1];
betas=[0 1];
der_compts=[0 1 2];
ilocal=0;
for ibeta=1:1:2
    for ialpha=1:1:2
        for ider_compt=1:1:3
            ilocal=ilocal+1;
            mx_var(ilocal,:)=Hermitian_poly_coef_serendipity(...
                alphas(ialpha),betas(ibeta),der_compts(ider_compt));
        end
    end
end
% mx_var;

for ialpha=1:1:2
    for igp=1:1:3
        ilocal=ilocal+1;
        mx_var(ilocal,:)=Hermitian_poly_coef_serendipity(...
            alphas(ialpha),gp(igp),1);
    end
end

for ibeta=1:1:2
    for igp=1:1:3
        ilocal=ilocal+1;
        mx_var(ilocal,:)=Hermitian_poly_coef_serendipity(...
            gp(igp),betas(ibeta),2);
    end
end


constraints=zeros(24,1);

constraints(1)=1;

size(mx_var)
size(constraints)

mx_var\constraints


hey=getting_odd_poly_coef(gp);

% mx_var

% size(mx_poly_coef)

hey(1,:)

return

z1=zeros(ny,nx);
z2=zeros(ny,nx);


z1(:,:)=aruto(1,:,:);
z2(:,:)=aruto(4,:,:);

% surf(z1+z2)

% surf()

nodal_weights=ones(1,12);

alpha=1;
beta=1;

laplacian=compute_laplacian(alpha,beta,mx_poly_coef,nodal_weights);

laplacian











warning('ending here')
return

nex=255;
ney=255;

dt=10^-10;

rx_to_ry_ratio=1;

nworkers=0;

entropy=1;

%===WETTIN===
wetting_boundary=0.000001; %Does not change with concentration
%============

ci=0.5;
T=0.3;
diff=6000;
n1=1;
n2=1;
fluc=10^-6;
dt_up=1.1;
dt_down=0.5;

%================= EXTRA FEATURE
err_tol=10;
%=================
%-------

tol_nr=10^-6;

export_data=1;
sig_dig=17;
jk_tol=3;

% get_twist([1 -1])*180/pi
% return

% === Set Ups =============================================================

ne=nex*ney;
nx=nex+1;
ny=ney+1;
n=nx*ny;

nthree=n*3;

gp=[(-sqrt(3/5)+1)/2 0.5 (sqrt(3/5)+1)/2];

two_chi_n1=0.5-entropy*(1-1/T);
two_chi_n1=2*two_chi_n1/n1;



diffT=diff*T;

% === Generating things ===================================================

xy=map_square_to_circle(nx,ny,n,rx_to_ry_ratio);

dto=dt;
time=0;

nworkers=cpu_initializer(nworkers);

% === 

nnz_=get_nnz_serendipity(nx,ny);

splitted_tasks=task_splitter(nworkers,nnz_);

[irow,icol]=setUpSparse_serendipity(nnz_,nx,ny,nthree);

disp('trying to determine weights')
[weights,determinants]=compute_weights_serendipity(...
    ne,xy,nx,nex,gp,rx_to_ry_ratio,nworkers);
disp('done with it')
% error('fdafdas')

% return

wTerms=generate_wTerms_serendipity(weights,ne);

if export_data==1
    output_path=folder_setUp_curved(nex,ney,rx_to_ry_ratio,diff,...
        ci,T,n1,n2);
end

%===============================INITIAL CONCENTRATION!!!
co=generate_co_serendipity(ci,fluc,nthree);

for gbf=2:3:nthree-1
    [node,type]=get_node_type_serendipity(gbf);
    [xth,yth]=get_x_y_th_serendipity(node,nx);
    if xth==1 || xth==nx || yth==1 || yth==nx
        if type==2
            co(gbf)=wetting_boundary;
        end
    end
end

% co=dlmread('Run_80');
%===========================

w=[5/18 4/9 5/18];

co=ave_c_adjuster(co,ci,weights,ne,nx,nex,w,determinants,nthree);

%===TEST=========

% [~,~,e1,e2]=get_equilibrium_pts(...
%     n1,n2,entropy,T,10^-6,10,...
%     10^-3,1);

if export_data==1
    dlmwrite([output_path 'time/iteration_1'],0);
    dlmwrite([output_path 'concentration/iteration_1'],co,'precision',sig_dig);
end

%========================while loop, supposedly=========================

c=co;

c_nodal=get_c_nodal_serendipity(c,nx,ny);
surf(c_nodal)
drawnow
% frame_cap=getframe(gcf);
second_frame=2;

%===================CONTINUATION CODE======================================

co=dlmread([pwd '/Results_ellipse/6000 0.5 0.3 1 wetting_boundary=0.000001/concentration/iteration_267']);
c=dlmread([pwd '/Results_ellipse/6000 0.5 0.3 1 wetting_boundary=0.000001/concentration/iteration_268']);

t1=dlmread([pwd '/Results_ellipse/6000 0.5 0.3 1 wetting_boundary=0.000001/time/iteration_267']);
t2=dlmread([pwd '/Results_ellipse/6000 0.5 0.3 1 wetting_boundary=0.000001/time/iteration_268']);

co=co';
c=c';

time=t2;

dto=t2-t1;
dt=dto;

second_frame=269; %PUT THE NEXT ITERATION NUMBER!

%==========================================================================


disp('Iterating')

for iteration=second_frame:1:50000

    cp=c+((c-co)/dto)*dt;
    coo=co;co=c;c=cp;
    jk=0;
    err=inf;
    erro=inf;
    
    while err>tol_nr
        
        jk=jk+1;
        
        if jk>jk_tol
            disp('Too many Newton Raphson iterations')
            jk=-1;
            break
        end
        
        conc_=get_conc_serendipity(c,weights,ne,nx,nex);
        
        disp('computing terms...')

        [terms_sf,terms_sj]=compute_terms_serendipity(...
            conc_,co,n1,n2,ne,two_chi_n1,dt,weights,nx,nex,...
            diffT);
        
%         [test_terms_sf,test_terms_sj]=test_get_terms(...
%             conc_,co,n1,n2,ne,two_chi_n1,dt,weights,nx,nex,...
%             diffT);

%         return

        disp('computing sf...')
        
        sf=sf_serendipity_2d(nthree,terms_sf,weights,determinants,...
            nx,ny,nex,w,diffT,wTerms);
        
%         test_sf=test_get_sf(nthree,test_terms_sf,weights,determinants,...
%             nx,ny,nex,w,diffT,wTerms);
        
        
%         sum(abs(test_sf-sf))
%         return

%         test_terms_sj=0;
        
        disp('computing sj...')
        sj=sj_serendipity_2d(nnz_,irow,icol,nx,ny,nex,...
            terms_sj,weights,determinants,w,diffT,wTerms,dt,...
            nworkers,splitted_tasks);
        
        
        
        
%         return
        
        
%         warning('hey')
%         return
        
        
        %================================ BC ===================
        disp('applying bc')
        

        %-------------------------------------------------------
        
        
        for ijk=1:1:nnz_
            [node,type]=get_node_type_serendipity(irow(ijk));
            [xth,yth]=get_x_y_th_serendipity(node,nx);
            if xth==1 || xth==nx || yth==1 || yth==ny
                if type==2
                    if irow(ijk)~=icol(ijk)
                        sj(ijk)=0;
                    elseif irow(ijk)==icol(ijk)
                        sj(ijk)=1;
                        sf(irow(ijk))=0;
                    end
                end
            end
        end
        
        disp('done applying bc')
        %========================================================
        
        disp('carrying out matrix division...')
        
        sj=sparse(irow,icol,sj);
        
        c_=sj\-sf';
        
        disp('evaluating error...')
        
        err=sqrt(sum(c_.^2))
        
        if err>=erro
            jk=-1;
            break
            
        elseif err>=err_tol
            warning('Error is too high')
            jk=-1;
            break
        else
            erro=err;
        end
        
        c=c_'+c;
        
%         dlmwrite([pwd '/temp/Run_' num2str(iteration)],c)
        %------------------
%         index=-2;
%         for node=1:1:n
%             index=index+3;
% %             [node,type]=get_node_type_serendipity(irow(ijk));
%             [xth,yth]=get_x_y_th_serendipity(node,nx);
%             if xth==1 || xth==nx || yth==1 || yth==nx
% %                 unit_vector=xy(node,:);
% %                 unit_vector=unit_vector./sqrt(xy(node,:).^2);
% %                 unit_vector
% %                 kitty=xy(node,1)*c(index+1)+xy(node,2)*c(index+2);
% %                 kitty
%                 c(index+1)
%             end
%         end
%         
%         return
        
        %-----------------
        
    end
    
%     warning('show jk')
%     
%     jk
    
    
    %=================TEST===============--
    
%     c=adjust_numerical_error(c,ci,n,nthree);
%     c_nodal=get_c_nodal_serendipity(c,nx,ny);
%     if min(min(c_nodal))<=e1 || max(max(c_nodal))>=e2
%         jk=-1;
%         warning('bad')
%     end
    
    %==========================================
    
    if jk==-1
        c=co;co=coo;
        dt=dt*dt_down;
        warning('going down now!')
        continue
    end
    
    c=ave_c_adjuster(c,ci,weights,ne,nx,nex,w,determinants,nthree);
    %--------------------------- Check new solution
    c_nodal=get_c_nodal_serendipity(c,nx,ny);
    if max(max(c_nodal))>=1 || min(min(c_nodal))<=0
        warning('Exceeded the limit')
        jk=-1;
    end
%     checker=check_horizontal_zig_zag(c_nodal,nx,ny);
%     if checker==0
%         error('Bad evaluation')
%         jk=-1;
%     end
    if jk==-1
        c=co;co=coo;
        dt=dt*dt_down;
        warning('going down now!')
        continue
    end
    %---------------------------
    
    
    disp('Finished iteration')
    disp(iteration)
    
    if export_data==1
        fprintf('Exporting Data...\n')
        dlmwrite([output_path 'time/iteration_' num2str(iteration)],time+dt,'precision',sig_dig);
        dlmwrite([output_path 'concentration/iteration_' num2str(iteration)],c,'precision',sig_dig);
    end
    
    %-------------------------------------- 
    
%     c=adjust_numerical_error(c,ci,n,nthree);
    
%     
%     summation=0;
%     for gbf=1:3:nthree-2
%         summation=summation+c(gbf);
%     end
%     summation/n-ci
%     
%     error('fsdafasg')
    
    
    for gbf=2:3:nthree-1
        [node,type]=get_node_type_serendipity(gbf);
        [xth,yth]=get_x_y_th_serendipity(node,nx);
        if xth==1 || xth==nx || yth==1 || yth==nx
            if type==2 && c(gbf)~=wetting_boundary
                c(gbf)
                error('hello there')
                
            end
        end
    end
%         if c(gbf)~=0
%             warning('hello there')
%             c(gbf)
%         end
%     end
    %--------------------------------------
%     return
        
%     x_mesh=xy(:,1);
%     y_mesh=xy(:,2);
    
%     x_mesh=reshape(x_mesh,[ny nx]);
%     y_mesh=reshape(y_mesh,[ny nx]);
%     for ix=1:1:nx

% [X,Y] = meshgrid(1:0.5:2,1:2);
% 
% X
% return

    if iteration==second_frame
        x_mesh=zeros(ny,nx);
        y_mesh=zeros(ny,nx);
        index=0;
        for iy=1:1:ny
            for ix=1:1:nx
                index=index+1;
                x_mesh(iy,ix)=xy(index,1);
                y_mesh(iy,ix)=xy(index,2);
            end
        end
    end
    
    
    
%     c_nodal=get_c_nodal_serendipity(c,nx,ny);
    surf(x_mesh,y_mesh,c_nodal,'EdgeColor','none')
    xlabel('x');ylabel('y')

%     surf(linspace(0,1,nx),linspace(0,1,ny),c_nodal)

    b=1/(pi*rx_to_ry_ratio)^0.5;
    a=(rx_to_ry_ratio*b);
    
    r=max([a b]);

    axis([-r r -r r 0 1])
    
    camlight
    
    drawnow
%     frame_cap=getframe(gcf);
    
    
    time=time+dt;
    
    
    if jk<=2
        dto=dt;
        dt=dt*dt_up;
        warning('upping')
    end
    
end
    
end

function c_new=adjust_numerical_error(c,ci,n,nthree)
c_new=c;
summation=0;
for i=1:3:nthree-2
    summation=summation+c(i);
end
deviation=summation/n-ci;
for i=1:3:nthree-2
    c_new(i)=c_new(i)-deviation;
end
end

function splitted_tasks=task_splitter(nworkers,nnz_)
if nworkers==1
    splitted_tasks=NaN;
elseif nworkers>1
    
    remainder=mod(nnz_,nworkers);
    ntasks=(nnz_-remainder)/nworkers;
    splitted_tasks=zeros(nworkers,2);
    index=0;
    for i=1:1:nworkers
        index=index+1;
        if remainder~=0
            splitted_tasks(i,1)=index;
            index=index+ntasks;
            splitted_tasks(i,2)=index;
            remainder=remainder-1;
        elseif remainder==0
            splitted_tasks(i,1)=index;
            index=index+ntasks-1;
            splitted_tasks(i,2)=index;
        end
    end
end
end

function c_nodal=get_c_nodal_serendipity(c,nx,ny)
c_nodal=zeros(ny,nx);
index=-2;
for ix=1:1:nx
    for iy=1:1:ny
        index=index+3;
        c_nodal(iy,ix)=c(index);
    end
end
end

function sj=sj_serendipity_2d(nnz_,irow,icol,nx,ny,nex,...
    terms_sj,weights,determinants,w,diffT,wTerms,dt,...
    nworkers,splitted_tasks)


sj=zeros(1,nnz_);

if nworkers==1
    for isj=1:1:nnz_
        sj(isj)=get_sj_one_element_serendipity(...
            analyze_relationship(irow(isj),icol(isj),nx,ny,nex),...
            terms_sj,weights,determinants,w,diffT,wTerms,dt);
    end
elseif nworkers>1
    temp_output=zeros(nworkers,splitted_tasks(1,2)-splitted_tasks(1,1)+1);
    
    
%     warning('TURN IT INTO PARFOR')
    parfor core=1:1:nworkers
%         temp_output_1_core=zeros(nworkers,splitted_tasks(core,2)-splitted_tasks(core,1)+1);
%         index_ends=splitted_tasks(core,:);
%         for isj=index_ends(1):1:index_ends(2)
%             temp_output_1_core(isj)=get_sj_one_element_serendipity(...
%                 analyze_relationship(irow(isj),icol(isj),nx,ny,nex),...
%                 terms_sj,weights,determinants,w,diffT,wTerms,dt);
%         end


%         try
        temp_output(core,:)=sj_parallel_computation_serendipity(...
            core,splitted_tasks,irow,icol,nx,ny,nex,...
            terms_sj,weights,determinants,w,diffT,wTerms,dt);
        
%         catch
%             
%             take=sj_parallel_computation_serendipity(...
%                 core,splitted_tasks,irow,icol,nx,ny,nex,...
%                 terms_sj,weights,determinants,w,diffT,wTerms,dt);
%             
%             
%             splitted_tasks
%             
%             size(take)
%             size(temp_output)
%             
%             error('stop here already')
%         end
    end
    
    for core=1:1:nworkers
        index_core=0;
        for isj=splitted_tasks(core,1):1:splitted_tasks(core,2)
            index_core=index_core+1;
            sj(isj)=temp_output(core,index_core);
        end
    end
end
end

function sol=sj_parallel_computation_serendipity(...
    core,splitted_tasks,irow,icol,nx,ny,nex,...
    terms_sj,weights,determinants,w,diffT,wTerms,dt)

sol=zeros(1,splitted_tasks(1,2)-splitted_tasks(1,1)+1);

index_ends=splitted_tasks(core,:);
index=0;
for isj=index_ends(1):1:index_ends(2)
    index=index+1;
    sol(index)=get_sj_one_element_serendipity(...
        analyze_relationship(irow(isj),icol(isj),nx,ny,nex),...
        terms_sj,weights,determinants,w,diffT,wTerms,dt);
end
        
end

function sol=get_sj_one_element_serendipity(mx_relationship,...
    terms_sj,weights,determinants,w,diffT,wTerms,dt)

sol=0;

% test_sol=0;

for ie=1:1:4
    if mx_relationship(ie,1)~=0
        
        e=mx_relationship(ie,1);
        ilocal=mx_relationship(ie,2);
        jlocal=mx_relationship(ie,3);
        
        
%         jlocal=mx_relationship(ie,2);
%         ilocal=mx_relationship(ie,3);
        
        for iy=1:1:3
            for ix=1:1:3
                
                %---Feb 28th-------------------------
                
%                 sol=sol+w(ix)*w(iy)*determinants(mx_relationship(ie,1),ix,iy)*(...
%                     ...
%                     weights(e,ilocal,1,ix,iy)*weights(e,jlocal,1,ix,iy)/dt...
%                     ...
%                     -weights(e,ilocal,1,ix,iy)...
%                     *(terms_sj(e,1,ix,iy)*weights(e,jlocal,1,ix,iy)...
%                     +terms_sj(e,2,ix,iy)*weights(e,jlocal,2,ix,iy)...
%                     +terms_sj(e,3,ix,iy)*weights(e,jlocal,3,ix,iy)...
%                     +terms_sj(e,4,ix,iy)*wTerms(e,jlocal,ix,iy))...
%                     ...
%                     +wTerms(e,ilocal,ix,iy)*wTerms(e,jlocal,ix,iy)...
%                     ...
%                     ...
%                     ...
%                     );
%                 
                %---Feb 27th-------------------------
                
                sol=sol+w(ix)*w(iy)*determinants(mx_relationship(ie,1),ix,iy)*(...
                    ...
                    ...
                    weights(e,ilocal,1,ix,iy)*weights(e,jlocal,1,ix,iy)/dt...
                    ...
                    ...
                    +weights(e,jlocal,1,ix,iy)...
                    *(terms_sj(e,1,ix,iy)*weights(e,ilocal,2,ix,iy)...
                    +terms_sj(e,2,ix,iy)*weights(e,ilocal,3,ix,iy))...
                    ...
                    ...
                    +terms_sj(e,3,ix,iy)...
                    *(weights(e,ilocal,2,ix,iy)*weights(e,jlocal,2,ix,iy)...
                    +weights(e,ilocal,3,ix,iy)*weights(e,jlocal,3,ix,iy))...
                    ...
                    ...
                    +wTerms(e,ilocal,ix,iy)*wTerms(e,jlocal,ix,iy)...
                    ...
                    ...
                    );
                
                %---------------------------------------

%                 sol=sol+w(ix)*w(iy)*determinants(mx_relationship(ie,1),ix,iy)*(...
%                     ...
%                     weights(mx_relationship(ie,1),mx_relationship(ie,2),1,ix,iy)...
%                     *weights(mx_relationship(ie,1),mx_relationship(ie,3),1,ix,iy)/dt...
%                     ...
%                     -diffT*weights(mx_relationship(ie,1),mx_relationship(ie,2),1,ix,iy)*(...
%                     weights(mx_relationship(ie,1),mx_relationship(ie,3),1,ix,iy)...
%                     *terms_sj(mx_relationship(ie,1),1,ix,iy)...
%                     +terms_sj(mx_relationship(ie,1),2,ix,iy)*(...
%                     terms_sj(mx_relationship(ie,1),3,ix,iy)...
%                     *weights(mx_relationship(ie,1),mx_relationship(ie,3),2,ix,iy)...
%                     +terms_sj(mx_relationship(ie,1),4,ix,iy)...
%                     *weights(mx_relationship(ie,1),mx_relationship(ie,3),3,ix,iy)...
%                     +weights(mx_relationship(ie,1),mx_relationship(ie,3),1,ix,iy)...
%                     *terms_sj(mx_relationship(ie,1),5,ix,iy)...
%                     ...
%                     )...
%                     +terms_sj(mx_relationship(ie,1),6,ix,iy)...
%                     *wTerms(mx_relationship(ie,1),mx_relationship(ie,3),ix,iy)...
%                     )...
%                     +wTerms(mx_relationship(ie,1),mx_relationship(ie,2),ix,iy)...
%                     *wTerms(mx_relationship(ie,1),mx_relationship(ie,3),ix,iy)...
%                     ...
%                     );
                
                
                %---------FINDING BUG HERE----------------------
                
%                 if sol~=test_sol
%                 
%                     weights(e,ilocal,1,ix,iy)*weights(e,jlocal,1,ix,iy)/dt
%                     
%                     weights(mx_relationship(ie,1),mx_relationship(ie,2),1,ix,iy)...
%                         *weights(mx_relationship(ie,1),mx_relationship(ie,3),1,ix,iy)/dt
% 
% 
%                     test1=w(ix)*w(iy)*determinants(mx_relationship(ie,1),ix,iy)*(...
%                         ...
%                         weights(mx_relationship(ie,1),mx_relationship(ie,2),1,ix,iy)...
%                         *weights(mx_relationship(ie,1),mx_relationship(ie,3),1,ix,iy)/dt...
%                         ...
%                         -diffT*weights(mx_relationship(ie,1),mx_relationship(ie,2),1,ix,iy)*(...
%                         weights(mx_relationship(ie,1),mx_relationship(ie,3),1,ix,iy)...
%                         *terms_sj(mx_relationship(ie,1),1,ix,iy)...
%                         +terms_sj(mx_relationship(ie,1),2,ix,iy)*(...
%                         terms_sj(mx_relationship(ie,1),3,ix,iy)...
%                         *weights(mx_relationship(ie,1),mx_relationship(ie,3),2,ix,iy)...
%                         +terms_sj(mx_relationship(ie,1),4,ix,iy)...
%                         *weights(mx_relationship(ie,1),mx_relationship(ie,3),3,ix,iy)...
%                         +weights(mx_relationship(ie,1),mx_relationship(ie,3),1,ix,iy)...
%                         *terms_sj(mx_relationship(ie,1),5,ix,iy)...
%                         ...
%                         )...
%                         +terms_sj(mx_relationship(ie,1),6,ix,iy)...
%                         *wTerms(mx_relationship(ie,1),mx_relationship(ie,3),ix,iy)...
%                         )...
%                         +wTerms(mx_relationship(ie,1),mx_relationship(ie,2),ix,iy)...
%                         *wTerms(mx_relationship(ie,1),mx_relationship(ie,3),ix,iy)...
%                         ...
%                         );
%                     
%                     test2=w(ix)*w(iy)*determinants(mx_relationship(ie,1),ix,iy)*(...
%                         ...
%                         weights(e,ilocal,1,ix,iy)*weights(e,jlocal,1,ix,iy)/dt...
%                         ...
%                         -weights(e,ilocal,ix,iy)...
%                         *(test_terms_sj(e,1,ix,iy)*weights(e,jlocal,1,ix,iy)...
%                         +test_terms_sj(e,2,ix,iy)*weights(e,jlocal,2,ix,iy)...
%                         +test_terms_sj(e,3,ix,iy)*weights(e,jlocal,3,ix,iy)...
%                         +test_terms_sj(e,4,ix,iy)*wTerms(e,jlocal,ix,iy))...
%                         ...
%                         +wTerms(e,ilocal,ix,iy)*wTerms(e,jlocal,ix,iy)...
%                         ...
%                         ...
%                         ...
%                         );
%                     
%                     warning('hello')
%                     weights(e,ilocal,ix,iy)
%                     weights(mx_relationship(ie,1),mx_relationship(ie,2),1,ix,iy)
%                     %-------
% %                     weights(mx_relationship(ie,1),mx_relationship(ie,2),1,ix,iy)...
% %                         *weights(mx_relationship(ie,1),mx_relationship(ie,3),1,ix,iy)/dt
% %                     weights(e,ilocal,1,ix,iy)*weights(e,jlocal,1,ix,iy)/dt
%                 
%                     warning('current difference')
%                     ([test1 test2])
%                     warning('[sol test_sol]')
%                     [sol test_sol]
%                     error('fdafda')
%                     
%                 end
                
                
            end
        end
    end
end


end

function mx_relationship=analyze_relationship(gbf1,gbf2,nx,ny,nex)
mx_relationship=ones(4,3);
[node1,type1]=get_node_type_serendipity(gbf1);
[node2,type2]=get_node_type_serendipity(gbf2);

[xth1,yth1]=get_x_y_th_serendipity(node1,nx);
[xth2,yth2]=get_x_y_th_serendipity(node2,nx);

%---------------------------------------
if xth1==1
    mx_relationship(1,1)=0;mx_relationship(3,1)=0;
elseif xth1==nx
    mx_relationship(2,1)=0;mx_relationship(4,1)=0;
end
if yth1==1
    mx_relationship(1,1)=0;mx_relationship(2,1)=0;
elseif yth1==ny
    mx_relationship(3,1)=0;mx_relationship(4,1)=0;
end
%----------
bottom_left=(yth1-2)*nex+xth1-1;
if mx_relationship(1,1)==1
    if xth2>xth1 || yth2>yth1
        mx_relationship(1,1)=0;
    else
        mx_relationship(1,1)=bottom_left;
        mx_relationship(1,2)=9+type1;
        if xth2<xth1
            if yth2<yth1
                mx_relationship(1,3)=type2;
            elseif yth2==yth1
                mx_relationship(1,3)=6+type2;
            end
        elseif xth2==xth1
            if yth2<yth1
                mx_relationship(1,3)=3+type2;
            elseif yth2==yth1
                mx_relationship(1,3)=9+type2;
            end
        end
    end
end
%----------------------------------------
if mx_relationship(2,1)==1
    if xth2<xth1 || yth2>yth1
        mx_relationship(2,1)=0;
    else
        mx_relationship(2,1)=bottom_left+1;
        mx_relationship(2,2)=6+type1;
        if xth2>xth1
            if yth2<yth1
                mx_relationship(2,3)=3+type2;
            elseif yth2==yth1
                mx_relationship(2,3)=9+type2;
            end
        elseif xth2==xth1
            if yth2<yth1
                mx_relationship(2,3)=type2;
            elseif yth2==yth1
                mx_relationship(2,3)=6+type2;
            end
        end
    end
end
%----------------------------------------------------
if mx_relationship(3,1)==1
    if xth2>xth1 || yth2<yth1
        mx_relationship(3,1)=0;
    else
        mx_relationship(3,1)=bottom_left+nex;
        mx_relationship(3,2)=3+type1;
        if xth2<xth1
            if yth2>yth1
                mx_relationship(3,3)=6+type2;
            elseif yth2==yth1
                mx_relationship(3,3)=type2;
            end
        elseif xth2==xth1
            if yth2>yth1
                mx_relationship(3,3)=9+type2;
            elseif yth2==yth1
                mx_relationship(3,3)=3+type2;
            end
        end
    end
end
if mx_relationship(4,1)==1
    if xth2<xth1 || yth2<yth1
        mx_relationship(4,1)=0;
    else
        mx_relationship(4,1)=bottom_left+nex+1;
        mx_relationship(4,2)=type1;
        if xth2>xth1
            if yth2>yth1
                mx_relationship(4,3)=9+type2;
            elseif yth2==yth1
                mx_relationship(4,3)=3+type2;
            end
        elseif xth2==xth1
            if yth2>yth1
                mx_relationship(4,3)=6+type2;
            elseif yth2==yth1
                mx_relationship(4,3)=type2;
            end
        end
    end
end
%--------------------------------------

end

function sf=sf_serendipity_2d(nthree,terms_sf,weights,determinants,...
    nx,ny,nex,w,diffT,wTerms)

sf=zeros(1,nthree);

for gbf=1:1:nthree
    [node,type]=get_node_type_serendipity(gbf);
    adjacency=get_adjacent_elements(node,nx,ny,nex);
    
    if adjacency(1)~=0
        sf_val=get_sf_one_element_serendipity(...
            adjacency(1),9+type,terms_sf,determinants,...
            weights,diffT,wTerms,w);
    elseif adjacency(1)==0
        sf_val=0;
    end
    if adjacency(2)~=0
        sf_val=sf_val+get_sf_one_element_serendipity(...
            adjacency(2),6+type,terms_sf,determinants,...
            weights,diffT,wTerms,w);
    end
    if adjacency(3)~=0
        sf_val=sf_val+get_sf_one_element_serendipity(...
            adjacency(3),3+type,terms_sf,determinants,...
            weights,diffT,wTerms,w);
    end
    if adjacency(4)~=0
        sf_val=sf_val+get_sf_one_element_serendipity(...
            adjacency(4),type,terms_sf,determinants,...
            weights,diffT,wTerms,w);
    end
    
    sf(gbf)=sf_val;
end

end

function sol=get_sf_one_element_serendipity(...
    e,ilocal,terms_sf,determinants,...
    weights,diffT,wTerms,w)
sol=0;

for iy=1:1:3
    for ix=1:1:3
        
        %---Feb 28th-------------------------
        
%         sol=sol+w(ix)*w(iy)*determinants(e,ix,iy)*(...
%             ...
%             terms_sf(e,1,ix,iy)*weights(e,ilocal,1,ix,iy)...
%             ...
%             -terms_sf(e,2,ix,iy)*weights(e,ilocal,1,ix,iy)...
%             ...
%             +terms_sf(e,3,ix,iy)*wTerms(e,ilocal,ix,iy)...
%             ...
%             );
        
        %---Feb 27th-------------------------
        
        sol=sol+w(ix)*w(iy)*determinants(e,ix,iy)*(...
            ...
            terms_sf(e,1,ix,iy)*weights(e,ilocal,1,ix,iy)...
            ...
            +terms_sf(e,2,ix,iy)*weights(e,ilocal,2,ix,iy)...
            +terms_sf(e,3,ix,iy)*weights(e,ilocal,3,ix,iy)...
            ...
            +terms_sf(e,4,ix,iy)*wTerms(e,ilocal,ix,iy)...
            ...
            );
        
        %---------------------------------------
        
%         sol=sol+w(ix)*w(iy)*determinants(e,ix,iy)*(...
%             ...
%             weights(e,ilocal,1,ix,iy)*(terms_sf(e,1,ix,iy)...
%             -diffT*terms_sf(e,2,ix,iy))...
%             ...
%             +terms_sf(e,3,ix,iy)*wTerms(e,ilocal,ix,iy)...
%             ...
%             );
    end
end

end

function wTerms=generate_wTerms_serendipity(weights,ne)
wTerms=zeros(ne,12,3,3);
for e=1:1:ne
    for iy=1:1:3
        for ix=1:1:3
            for type=1:1:12
                wTerms(e,type,ix,iy)=...
                    weights(e,type,4,ix,iy)+weights(e,type,5,ix,iy);
            end
        end
    end
end
end

function [terms_sf,terms_sj]=compute_terms_serendipity(...
    conc_,co,n1,n2,ne,two_chi_n1,dt,weights,nx,nex,...
    diffT)
% terms_sf=zeros(ne,3,3,3);

terms_sf=zeros(ne,4,3,3);

terms_sj=zeros(ne,3,3,3);

conco_=get_conc_first_order_serendipity(co,weights,ne,nx,nex);

for e=1:1:ne
    for iy=1:1:3
        for ix=1:1:3
            
            %---Feb 28th-------------------------
            
%             c2=1-conc_(e,1,ix,iy);
%             
%             d2f_dc2=(conc_(e,1,ix,iy)*n1)^-1+(c2*n2)^-1-two_chi_n1;
%             d3f_dc3=-(conc_(e,1,ix,iy)^2*n1)^-1+(c2^2*n2)^-1;
%             d4f_dc4=2*((conc_(e,1,ix,iy)^3*n1)^-1+(c2^3*n2)^-1);
%             
%             grad_c_sq=conc_(e,2,ix,iy)^2+conc_(e,3,ix,iy)^2;
%             laplacian=conc_(e,4,ix,iy)+conc_(e,5,ix,iy);
%             
%             
%             terms_sf(e,1,ix,iy)=(conc_(e,1,ix,iy)-conco_(e,ix,iy))/dt;
%             terms_sf(e,2,ix,iy)=...
%                 diffT*(d3f_dc3*grad_c_sq+d2f_dc2*laplacian);
%             terms_sf(e,3,ix,iy)=laplacian;
%             
%             
%             terms_sj(e,1,ix,iy)=diffT*(d4f_dc4*grad_c_sq+d3f_dc3*laplacian);
%             
%             sub=2*diffT*d3f_dc3;
%             
%             terms_sj(e,2,ix,iy)=sub*conc_(e,2,ix,iy);
%             terms_sj(e,3,ix,iy)=sub*conc_(e,3,ix,iy);
%             
%             terms_sj(e,4,ix,iy)=diffT*d2f_dc2;
%             
%             %test
%             test_sf=zeros(1,3);
%             test_sj=zeros(1,4);
%             for ijk=1:1:3
%                 test_sf(ijk)=terms_sf(e,ijk,ix,iy);
%             end
%             for ijk=1:1:4
%                 test_sj(ijk)=terms_sj(e,ijk,ix,iy);
%             end
            
            
            %---Feb 27th-------------------------
            
            c2=1-conc_(e,1,ix,iy);
            
            terms_sf(e,1,ix,iy)=(conc_(e,1,ix,iy)-conco_(e,ix,iy))/dt;
            
            terms_sj(e,3,ix,iy)=diffT...
                *((conc_(e,1,ix,iy)*n1)^-1+(c2*n2)^-1-two_chi_n1);
            
            terms_sf(e,2,ix,iy)=terms_sj(e,3,ix,iy)*conc_(e,2,ix,iy);
            terms_sf(e,3,ix,iy)=terms_sj(e,3,ix,iy)*conc_(e,3,ix,iy);
            
            terms_sf(e,4,ix,iy)=conc_(e,4,ix,iy)+conc_(e,5,ix,iy);
            
            sub=diffT*(-(conc_(e,1,ix,iy)^2*n1)^-1+(c2^2*n2)^-1);
            
            terms_sj(e,1,ix,iy)=sub*conc_(e,2,ix,iy);
            terms_sj(e,2,ix,iy)=sub*conc_(e,3,ix,iy);
            
            
            %------------------------------------------
            
%             c2=1-conc_(e,1,ix,iy);
%             
%             sub1=-(conc_(e,1,ix,iy)^2*n1)^-1+(c2^2*n2)^-1;
%             sub2=conc_(e,2,ix,iy)^2+conc_(e,3,ix,iy)^2;
%             sub3=(conc_(e,1,ix,iy)*n1)^-1+(c2*n2)^-1-two_chi_n1;
%             sub4=conc_(e,4,ix,iy)+conc_(e,5,ix,iy);
%             
%             %---Controversal-------------
%             terms_sf(e,1,ix,iy)=(conc_(e,1,ix,iy)-conco_(e,ix,iy))/dt;
%             %-------------------------------
%             terms_sf(e,2,ix,iy)=sub1*sub2+sub3*sub4;
%             terms_sf(e,3,ix,iy)=sub4;
%             
%             terms_sj(e,1,ix,iy)=...
%                 2*((conc_(e,1,ix,iy)^3*n1)^-1+(c2^3*n2)^-1)...
%                 *sub2;
%             terms_sj(e,2,ix,iy)=sub1;
%             terms_sj(e,3,ix,iy)=2*conc_(e,2,ix,iy);
%             terms_sj(e,4,ix,iy)=2*conc_(e,3,ix,iy);
%             terms_sj(e,5,ix,iy)=sub4;
%             terms_sj(e,6,ix,iy)=sub3;
%             
% %             terms(e,ix,iy,1)=sub1*sub2;
% %             terms(e,ix,iy,2)=sub3*sub4;
% %             terms(e,ix,iy,3)=sub4;
% %             terms(e,ix,iy,4)=...
% %                 2*((conc_(e,1,ix,iy)^3*n1)^-1+(c2^3*n2)^-1)...
% %                 *sub2;
% %             terms(e,ix,iy,5)=sub1;
% %             terms(e,ix,iy,6)=sub3;
            %---------------------------------------
            
            % ------------TEST-------------
            
% %             test1=d4f_dc4*grad_c_sq;
% %             [test1 terms_sj(e,1,ix,iy)]
%             test=diffT*(terms_sj(e,1,ix,iy)+terms_sj(e,2,ix,iy)*terms_sj(e,5,ix,iy));
%             [test test_sj(1)];
%             test=diffT*terms_sj(e,2,ix,iy)*terms_sj(e,3,ix,iy);
%             [test test_sj(2)];
%             test=diffT*terms_sj(e,2,ix,iy)*terms_sj(e,4,ix,iy);
%             [test test_sj(3)];
%             test=diffT*terms_sj(e,6,ix,iy);
%             [test test_sj(4)];
%             test=terms_sf(e,1,ix,iy);
%             [test test_sf(1)];
%             test=diffT*terms_sf(e,2,ix,iy);
%             [test test_sf(2)];
%             test=terms_sf(e,3,ix,iy);
%             [test test_sf(3)];
            
            % ------------------------------
            
        end
    end
end
end

function sol=get_conc_first_order_serendipity(c,weights,ne,nx,nex)
sol=zeros(ne,3,3);
for e=1:1:ne
    gbfs=get_gbfs_elemental(e,nx,nex);
    for iy=1:1:3
        for ix=1:1:3
            for type=1:1:12
                sol(e,ix,iy)=sol(e,ix,iy)...
                    +c(gbfs(type))...
                    *weights(e,type,1,ix,iy);
            end
        end
    end
end
end

function sol=get_conc_serendipity(c,weights,ne,nx,nex)
sol=zeros(ne,5,3,3);
for e=1:1:ne
    gbfs=get_gbfs_elemental(e,nx,nex);
    for iy=1:1:3
        for ix=1:1:3
            %-------------
            for order=1:1:5
                for type=1:1:12
                    sol(e,order,ix,iy)=sol(e,order,ix,iy)...
                        +c(gbfs(type))...
                        *weights(e,type,order,ix,iy);
                end
            end
            %---------------------
        end
    end
end
end

function gbfs=get_gbfs_elemental(e,nx,nex)
gbfs=zeros(1,12);
nodes=get_eNodes_serendipity(e,nx,nex);
gbfs(1:1:6)=nodes(1)*3-2:1:nodes(2)*3;
gbfs(7:1:12)=nodes(3)*3-2:1:nodes(4)*3;
end

function sol=get_adjacent_elements(node,nx,ny,nex)
sol=ones(1,4);
[xth,yth]=get_x_y_th_serendipity(node,nx);
if xth==1
    sol(1)=0;sol(3)=0;
elseif xth==nx
    sol(2)=0;sol(4)=0;
end
if yth==1
    sol(1)=0;sol(2)=0;
elseif yth==ny
    sol(3)=0;sol(4)=0;
end
bottom_left=(yth-2)*nex+xth-1;
if sol(1)~=0
    sol(1)=bottom_left;
end
if sol(2)~=0
    sol(2)=bottom_left+1;
end
if sol(3)~=0
    sol(3)=bottom_left+nex;
end
if sol(4)~=0
    sol(4)=bottom_left+nex+1;
end
end

function [node,type]=get_node_type_serendipity(gbf)
type=sp_mod(gbf,3);
node=(gbf-type)/3+1;
end

function [xth,yth]=get_x_y_th_serendipity(node,nx)
xth=sp_mod(node,nx);
yth=(node-xth)/nx+1;
end

function co=generate_co_serendipity(ci,fluc,nthree)
%THIS HAS NOT BEEN ADJUSTED TO ENSURE AVERAGE CONCENTRATION!
co=zeros(1,nthree);
for i=1:3:nthree-2
    co(i)=ci+fluc*(rand(1,1)*2-1);
end
end

%============================================================

function twist=get_twist(v)
if v(2)>0
    if v(1)>0
        twist=atan(v(2)/v(1));
    elseif v(1)<0
        twist=pi-atan(-v(2)/v(1));
    elseif v(1)==0
        twist=pi/2;
    end
elseif v(2)<0
     if v(1)>0
        twist=2*pi-atan(-v(2)/v(1));
    elseif v(1)<0
        twist=pi+atan(v(2)/v(1));
    elseif v(1)==0
        twist=3*pi/2;
     end
elseif v(2)==0
    if v(1)>0
        twist=0;
    elseif v(1)<0
        twist=pi;
    elseif v(1)==0
        warning('Vector must have a length not equal zero!')
        twist=0;
    end
end
end

function [weights,determinants]=compute_weights_serendipity(...
    ne,xy,nx,nex,gp,rx_to_ry_ratio,nworkers)
weights=zeros(ne,12,5,3,3);
determinants=zeros(ne,3,3);


%=======TEST TWISTER PART 1/2=========
% warning('delete this')
% gp(1)=0;
% gp(3)=1;
%===================================

%-------New info--------------
coords=[0 0
    1 0
    0 1
    1 1];
mx_coef=generate_matrix_serendipity(coords);
mx_constraints=mx_constraints_assembler([0 0 0 0]);
mx_poly_coef=generate_poly_coef(mx_coef,mx_constraints);
mx_P=zeros(3,3,12,6);
for ialpha=1:1:3
    for ibeta=1:1:3
        mx_P(ialpha,ibeta,:,:)=P(gp(ialpha),gp(ibeta),mx_poly_coef);
    end
end
%-----------------------------

[allocation,extra,load_less,load_more]=wise_task_splitter_v2(ne,nworkers);

temp1=zeros(nworkers,load_more,12,5,3,3);
temp2=zeros(nworkers,load_more,3,3);

parfor worker=1:1:nworkers
    if worker<=extra
        [sol1,sol2]=compute_weights_serendipity_assist(...
            allocation(worker,:),load_more,load_more,nx,nex,xy,gp,mx_P);
        temp1(worker,:,:,:,:,:)=sol1;
        temp2(worker,:,:,:)=sol2;
    elseif worker>extra
        [sol1,sol2]=compute_weights_serendipity_assist(...
            allocation(worker,:),load_less,load_more,nx,nex,xy,gp,mx_P);
        temp1(worker,:,:,:,:,:)=sol1;
        temp2(worker,:,:,:)=sol2;
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
        weights(index,:,:,:,:)=temp1(worker,i,:,:,:,:);
        determinants(index,:,:)=temp2(worker,i,:,:);
    end
end

% return
% %-------------------------

% for e=1:1:ne
%     nodes=get_eNodes_serendipity(e,nx,nex);
%     for inode=1:1:4
%         x(inode)=xy(nodes(inode),1);
%         y(inode)=xy(nodes(inode),2);
%     end
%     %================TWIST ADDER HEREEEEEEEEEEE TESTING!======================
%     
%     twists=twist_radially(x,y);
%     
% %     twists=[0 0 0 0];
% %     exth=sp_mod(e,nex);
% %     r_sq=rx_to_ry_ratio^2;
% %     if exth==1
% %         twists(1)=get_twist(-[1 r_sq*xy(nodes(1),2)/xy(nodes(1),1)]);
% %         twists(3)=get_twist(-[1 r_sq*xy(nodes(3),2)/xy(nodes(3),1)]);
% %     elseif exth==nex
% %         twists(2)=get_twist([1 r_sq*xy(nodes(2),2)/xy(nodes(2),1)]);
% %         twists(4)=get_twist([1 r_sq*xy(nodes(4),2)/xy(nodes(4),1)]);
% %     end
% %     eyth=(e-exth)/nex+1;
% %     if eyth==1
% %         if exth==1
% %             twists(2)=get_twist(-[1 r_sq*xy(nodes(2),2)/xy(nodes(2),1)]);
% %         elseif exth==nex
% %             twists(1)=get_twist([1 r_sq*xy(nodes(1),2)/xy(nodes(1),1)]);
% %         else
% %             if xy(nodes(1),1)>=0
% %                 twists(1)=get_twist([1 r_sq*xy(nodes(1),2)/xy(nodes(1),1)]);
% %                 twists(2)=get_twist([1 r_sq*xy(nodes(2),2)/xy(nodes(2),1)]);
% %             else
% %                 twists(1)=get_twist(-[1 r_sq*xy(nodes(1),2)/xy(nodes(1),1)]);
% %                 if xy(nodes(2),1)>=0
% %                     twists(2)=get_twist([1 r_sq*xy(nodes(2),2)/xy(nodes(2),1)]);
% %                 else
% %                     twists(2)=get_twist(-[1 r_sq*xy(nodes(2),2)/xy(nodes(2),1)]);
% %                 end
% %             end
% %         end
% %     elseif eyth==ne/nex
% %         if exth==1
% %             twists(4)=get_twist(-[1 r_sq*xy(nodes(4),2)/xy(nodes(4),1)]);
% %         elseif exth==nex
% %             twists(3)=get_twist([1 r_sq*xy(nodes(3),2)/xy(nodes(3),1)]);
% %         else
% %             if xy(nodes(1),1)>=0
% %                 twists(3)=get_twist([1 r_sq*xy(nodes(3),2)/xy(nodes(3),1)]);
% %                 twists(4)=get_twist([1 r_sq*xy(nodes(4),2)/xy(nodes(4),1)]);
% %             else
% %                 twists(3)=get_twist(-[1 r_sq*xy(nodes(3),2)/xy(nodes(3),1)]);
% %                 if xy(nodes(2),1)>=0
% %                     twists(4)=get_twist([1 r_sq*xy(nodes(4),2)/xy(nodes(4),1)]);
% %                 else
% %                     twists(4)=get_twist(-[1 r_sq*xy(nodes(4),2)/xy(nodes(4),1)]);
% %                 end
% %             end
% %         end
% %     end
%     %==========================================================================
%     %--------------
% %     for i=1:1:4
% %         if twists(i)~=0
% %             index=index+1;
% %             take(index)=twists(i);
% %             xp=xy(nodes(i),1);
% %             yp=xy(nodes(i),2);
% %             plot([xp xp+cos(twists(i))],[yp yp+sin(twists(i))])
% %         end
% %     end
% %     if sum(twists)~=0
% % 
% %         plot([xy(nodes(1),1) xy(nodes(2),1)],[xy(nodes(1),2) xy(nodes(2),2)])
% %         plot([xy(nodes(1),1) xy(nodes(3),1)],[xy(nodes(1),2) xy(nodes(3),2)])
% %         plot([xy(nodes(2),1) xy(nodes(4),1)],[xy(nodes(2),2) xy(nodes(4),2)])
% %         plot([xy(nodes(3),1) xy(nodes(4),1)],[xy(nodes(3),2) xy(nodes(4),2)])
% %     end
%     %==========================================================
%     
%     
%     weights(e,:,:,:,:)=get_weight_global(gp,x,y,twists,mx_P);
%     
%     
%     for iy=1:1:3
%         for ix=1:1:3
% %             weights(e,:,:,ix,iy)=get_weight_global(gp(ix),gp(iy),x,y,twists);
%             %------------------Determinants-----------
%             
%             determinants(e,ix,iy)=det([...
%                 mapper_dx_or_y_dalpha(gp(iy),x),mapper_dx_or_y_dbeta(gp(ix),x)...
%                 ;mapper_dx_or_y_dalpha(gp(iy),y),mapper_dx_or_y_dbeta(gp(ix),y)]);
%             
%             %--------------------
%         end
%     end
%         
%     
%     %=======TEST TWISTER PART 2/2=========
% %     test=zeros(4,3);
% %     v=zeros(4,2);
% %     v(1,:)=[weights(e,2,2,1,1),weights(e,2,3,1,1)];
% %     v(2,:)=[weights(e,5,2,3,1),weights(e,5,3,3,1)];
% %     v(3,:)=[weights(e,8,2,1,3),weights(e,8,3,1,3)];
% %     v(4,:)=[weights(e,11,2,3,3),weights(e,11,3,3,3)];
% %     for inode=1:1:4
% %         test(inode,1)=twists(inode);
% %         test(inode,2)=get_twist(v(inode,:));
% %     end
% %     e
% %     test
%     %====================================
%     
% end
% 
% 
% 
% % hold off
% % grid on
% % grid minor
% % take'
% % error('fdafa')


end

function nodes=get_eNodes_serendipity(e,nx,nex)
nodes=zeros(1,4);
[xth1,yth1]=get_xyth1_serendipity(e,nex);
nodes(1)=(yth1-1)*nx+xth1;
nodes(2)=nodes(1)+1;
nodes(3)=nodes(1)+nx;
nodes(4)=nodes(3)+1;
end

function [irow,icol]=setUpSparse_serendipity(nnz_,nx,ny,nthree)
irow=zeros(1,nnz_);
icol=zeros(1,nnz_);
nx3=nx*3;
index=0;
for gbf1=1:1:3
    for gbf2=1:1:6
        index=index+1;
        irow(index)=gbf1;
        icol(index)=gbf2;
    end
    starter=nx*3;
    for i=1:1:6
        index=index+1;
        irow(index)=gbf1;
        icol(index)=starter+i;
    end
end
for e=2:1:nx-1
    starter=(e-2)*3;
    starter_=starter+3*nx;
    for i=1:1:3
        gbf1=(e-1)*3+i;
        for j=1:1:9
            index=index+1;
            irow(index)=gbf1;
            icol(index)=starter+j;
        end
        for j=1:1:9
            index=index+1;
            irow(index)=gbf1;
            icol(index)=starter_+j;
        end
    end
end
for gbf1=nx3-2:1:nx3
    for gbf2=nx3-5:1:nx3
        index=index+1;
        irow(index)=gbf1;
        icol(index)=gbf2;
    end
    for gbf2=nx3*2-5:1:nx3*2
        index=index+1;
        irow(index)=gbf1;
        icol(index)=gbf2;
    end
end
for yth=2:1:ny-1
    pre_gbf1=(yth-1)*nx3;
    pre_gbf2_lower=(yth-2)*nx3;
    pre_gbf2_mid=pre_gbf2_lower+nx3;
    pre_gbf2_higher=pre_gbf2_mid+nx3;
    for gbf1=pre_gbf1+1:1:pre_gbf1+3
        for i=1:1:6
            index=index+1;
            irow(index)=gbf1;
            icol(index)=pre_gbf2_lower+i;
        end
        for i=1:1:6
            index=index+1;
            irow(index)=gbf1;
            icol(index)=pre_gbf2_mid+i;
        end
        for i=1:1:6
            index=index+1;
            irow(index)=gbf1;
            icol(index)=pre_gbf2_higher+i;
        end
    end
    pre_gbf2_lower=pre_gbf2_lower-3;
    pre_gbf2_mid=pre_gbf2_mid-3;
    pre_gbf2_higher=pre_gbf2_higher-3;
    for xth=2:1:nx-1
        pre_gbf1=pre_gbf1+3;
        pre_gbf2_lower=pre_gbf2_lower+3;
        pre_gbf2_mid=pre_gbf2_mid+3;
        pre_gbf2_higher=pre_gbf2_higher+3;
        for igbf1=1:1:3
            gbf1=pre_gbf1+igbf1;
            for igbf2=1:1:9
                index=index+1;
                irow(index)=gbf1;
                icol(index)=pre_gbf2_lower+igbf2;
            end
            for igbf2=1:1:9
                index=index+1;
                irow(index)=gbf1;
                icol(index)=pre_gbf2_mid+igbf2;
            end
            for igbf2=1:1:9
                index=index+1;
                irow(index)=gbf1;
                icol(index)=pre_gbf2_higher+igbf2;
            end
        end
    end
    pre_gbf1=pre_gbf1+3;
    pre_gbf2_lower=pre_gbf2_lower+3;
    pre_gbf2_mid=pre_gbf2_mid+3;
    pre_gbf2_higher=pre_gbf2_higher+3;
    for igbf1=1:1:3
        gbf1=pre_gbf1+igbf1;
        for igbf2=1:1:6
            index=index+1;
            irow(index)=gbf1;
            icol(index)=pre_gbf2_lower+igbf2;
        end
        for igbf2=1:1:6
            index=index+1;
            irow(index)=gbf1;
            icol(index)=pre_gbf2_mid+igbf2;
        end
        for igbf2=1:1:6
            index=index+1;
            irow(index)=gbf1;
            icol(index)=pre_gbf2_higher+igbf2;
        end
    end
end
pre_gbf1=pre_gbf1+3;
pre_gbf2_lower=pre_gbf2_lower+6;
pre_gbf2_mid=pre_gbf2_mid+6;
for igbf1=1:1:3
    gbf1=pre_gbf1+igbf1;
    for igbf2=1:1:6
        index=index+1;
        irow(index)=gbf1;
        icol(index)=pre_gbf2_lower+igbf2;
    end
    for igbf2=1:1:6
        index=index+1;
        irow(index)=gbf1;
        icol(index)=pre_gbf2_mid+igbf2;
    end
end
pre_gbf2_lower=pre_gbf2_lower-3;
pre_gbf2_mid=pre_gbf2_mid-3;
for xth=2:1:nx-1
    pre_gbf1=pre_gbf1+3;
    pre_gbf2_lower=pre_gbf2_lower+3;
    pre_gbf2_mid=pre_gbf2_mid+3;
    for igbf1=1:1:3
        gbf1=pre_gbf1+igbf1;
        for igbf2=1:1:9
            index=index+1;
            irow(index)=gbf1;
            icol(index)=pre_gbf2_lower+igbf2;
        end
        for igbf2=1:1:9
            index=index+1;
            irow(index)=gbf1;
            icol(index)=pre_gbf2_mid+igbf2;
        end
    end
end
pre_gbf2_lower=pre_gbf2_lower+3;
pre_gbf2_mid=pre_gbf2_mid+3;
for gbf1=nthree-2:1:nthree
    for igbf2=1:1:6
        index=index+1;
        irow(index)=gbf1;
        icol(index)=pre_gbf2_lower+igbf2;
    end
    for igbf2=1:1:6
        index=index+1;
        irow(index)=gbf1;
        icol(index)=pre_gbf2_mid+igbf2;
    end
end
end

function nnz_=get_nnz_serendipity(nx,ny)
nx_=nx-2;ny_=ny-2;
nnz_=144+108*(nx_+ny_)+81*nx_*ny_;
end

%=============================================MAR12TH====================

function [xth1,yth1]=get_xyth1_serendipity(e,nex)
xth1=sp_mod(e,nex);
yth1=(e-xth1)/nex+1;
end

function sol=sp_mod(a,b)
sol=mod(a,b);
if sol==0
    sol=b;
end
end

function xy=map_square_to_circle(nx,ny,n,rx_to_ry_ratio)

%============TEMPORARY============================= REMOVE LATER!
% xy=zeros(n,2);
% x=linspace(0,1,nx);
% y=linspace(0,1,ny);
% node=0;
% for iy=1:1:ny
%     for ix=1:1:nx
%         node=node+1;
%         xy(node,1)=x(ix);
%         xy(node,2)=y(iy);
%     end
% end
% return
%==============================================


xy=zeros(n,2);

x=linspace(-1,1,nx);
y=linspace(-1,1,ny);

node=0;
for iy=1:1:ny
    for ix=1:1:nx
        node=node+1;
        xy(node,1)=x(ix)*sqrt(1-y(iy)^2/2);
        xy(node,2)=y(iy)*sqrt(1-x(ix)^2/2);
    end
end

ry=(rx_to_ry_ratio*pi)^-0.5;
rx=rx_to_ry_ratio*ry;

xy(:,1)=xy(:,1)*rx;
xy(:,2)=xy(:,2)*ry;

end

function sol=get_weight_global(gp,x,y,twists,mx_P)
%THIS IS THE FINAL FUNCTION WITH THE TWIST!!!

sol=zeros(12,5,3,3);

derivatives_local=zeros(5,1);

for ialpha=1:1:3
    for ibeta=1:1:3
        
        dx_dalpha=mapper_dx_or_y_dalpha(gp(ibeta),x);
        dy_dalpha=mapper_dx_or_y_dalpha(gp(ibeta),y);
        dx_dbeta=mapper_dx_or_y_dbeta(gp(ialpha),x);
        dy_dbeta=mapper_dx_or_y_dbeta(gp(ialpha),y);
        d2x_dalphadbeta=mapper_d2x_or_y_dalphadbeta(x);
        d2y_dalphadbeta=mapper_d2x_or_y_dalphadbeta(y);
        
        mx_chain=mx_chain_assembler(dx_dalpha,dy_dalpha,dx_dbeta,dy_dbeta,...
            d2x_dalphadbeta,d2y_dalphadbeta);
        
        ilocal1=-2;
        ilocal2=-1;
        ilocal3=0;

        for node=1:1:4

            ilocal1=ilocal1+3;
            ilocal2=ilocal2+3;
            ilocal3=ilocal3+3;

            sol(ilocal1,1,ialpha,ibeta)=mx_P(ialpha,ibeta,ilocal1,1);

            for i=1:1:5
                derivatives_local(i)=mx_P(ialpha,ibeta,ilocal1,i+1);
            end

            derivatives_global=mx_chain\derivatives_local;

            for i=1:1:4
                sol(ilocal1,i+1,ialpha,ibeta)=derivatives_global(i);
            end

            [der_1st,der_2nd,der_3rd,der_mixed]=twisted_chain(cos(twists(node)),sin(twists(node)),...
                dx_dalpha,dy_dalpha,dx_dbeta,dy_dbeta,d2x_dalphadbeta,d2y_dalphadbeta);

            sol(ilocal2,1,ialpha,ibeta)=...
                mx_P(ialpha,ibeta,ilocal2,1)*der_1st(1)...
                +mx_P(ialpha,ibeta,ilocal3,1)*der_1st(2);

            sol(ilocal3,1,ialpha,ibeta)=...
                mx_P(ialpha,ibeta,ilocal2,1)*der_1st(3)...
                +mx_P(ialpha,ibeta,ilocal3,1)*der_1st(4);

            %---------------------------------------------

            derivatives_local(1)=...
                mx_P(ialpha,ibeta,ilocal2,2)*der_1st(1)...
                +mx_P(ialpha,ibeta,ilocal2,1)*der_2nd(1)...
                +mx_P(ialpha,ibeta,ilocal3,2)*der_1st(2)...
                +mx_P(ialpha,ibeta,ilocal3,1)*der_mixed(1);


            derivatives_local(2)=...
                mx_P(ialpha,ibeta,ilocal2,3)*der_1st(1)...
                +mx_P(ialpha,ibeta,ilocal2,1)*der_mixed(1)...
                +mx_P(ialpha,ibeta,ilocal3,3)*der_1st(2)...
                +mx_P(ialpha,ibeta,ilocal3,1)*der_2nd(2);


            derivatives_local(3)=...
                mx_P(ialpha,ibeta,ilocal2,4)*der_1st(1)...
                +2*mx_P(ialpha,ibeta,ilocal2,2)*der_2nd(1)...
                +mx_P(ialpha,ibeta,ilocal2,1)*der_3rd(1)...
                ...
                +mx_P(ialpha,ibeta,ilocal3,4)*der_1st(2)...
                +2*mx_P(ialpha,ibeta,ilocal3,2)*der_mixed(1)...
                +mx_P(ialpha,ibeta,ilocal3,1)*der_mixed(3);


            derivatives_local(4)=...
                mx_P(ialpha,ibeta,ilocal2,5)*der_1st(1)...
                +2*mx_P(ialpha,ibeta,ilocal2,3)*der_mixed(1)...
                +mx_P(ialpha,ibeta,ilocal2,1)*der_mixed(4)...
                ...
                +mx_P(ialpha,ibeta,ilocal3,5)*der_1st(2)...
                +2*mx_P(ialpha,ibeta,ilocal3,3)*der_2nd(2)...
                +mx_P(ialpha,ibeta,ilocal3,1)*der_3rd(2);


            derivatives_local(5)=...
                mx_P(ialpha,ibeta,ilocal2,6)*der_1st(1)...
                +mx_P(ialpha,ibeta,ilocal2,2)*der_mixed(1)...
                +mx_P(ialpha,ibeta,ilocal2,3)*der_2nd(1)...
                +mx_P(ialpha,ibeta,ilocal2,1)*der_mixed(3)...
                ...
                +mx_P(ialpha,ibeta,ilocal3,6)*der_1st(2)...
                +mx_P(ialpha,ibeta,ilocal3,2)*der_2nd(2)...
                +mx_P(ialpha,ibeta,ilocal3,3)*der_mixed(1)...
                +mx_P(ialpha,ibeta,ilocal3,1)*der_mixed(4);

            derivatives_global=mx_chain\derivatives_local;

            for i=1:1:4
                sol(ilocal2,i+1,ialpha,ibeta)=derivatives_global(i);
            end

            derivatives_local(1)=...
                mx_P(ialpha,ibeta,ilocal2,2)*der_1st(3)...
                +mx_P(ialpha,ibeta,ilocal2,1)*der_2nd(3)...
                +mx_P(ialpha,ibeta,ilocal3,2)*der_1st(4)...
                +mx_P(ialpha,ibeta,ilocal3,1)*der_mixed(2);


            derivatives_local(2)=...
                mx_P(ialpha,ibeta,ilocal2,3)*der_1st(3)...
                +mx_P(ialpha,ibeta,ilocal2,1)*der_mixed(2)...
                +mx_P(ialpha,ibeta,ilocal3,3)*der_1st(4)...
                +mx_P(ialpha,ibeta,ilocal3,1)*der_2nd(4);


            derivatives_local(3)=...
                mx_P(ialpha,ibeta,ilocal2,4)*der_1st(3)...
                +2*mx_P(ialpha,ibeta,ilocal2,2)*der_2nd(3)...
                +mx_P(ialpha,ibeta,ilocal2,1)*der_3rd(3)...
                ...
                +mx_P(ialpha,ibeta,ilocal3,4)*der_1st(4)...
                +2*mx_P(ialpha,ibeta,ilocal3,2)*der_mixed(2)...
                +mx_P(ialpha,ibeta,ilocal3,1)*der_mixed(5);


            derivatives_local(4)=...
                mx_P(ialpha,ibeta,ilocal2,5)*der_1st(3)...
                +2*mx_P(ialpha,ibeta,ilocal2,3)*der_mixed(2)...
                +mx_P(ialpha,ibeta,ilocal2,1)*der_mixed(6)...
                ...
                +mx_P(ialpha,ibeta,ilocal3,5)*der_1st(4)...
                +2*mx_P(ialpha,ibeta,ilocal3,3)*der_2nd(4)...
                +mx_P(ialpha,ibeta,ilocal3,1)*der_3rd(4);


            derivatives_local(5)=...
                mx_P(ialpha,ibeta,ilocal2,6)*der_1st(3)...
                +mx_P(ialpha,ibeta,ilocal2,2)*der_mixed(2)...
                +mx_P(ialpha,ibeta,ilocal2,3)*der_2nd(3)...
                +mx_P(ialpha,ibeta,ilocal2,1)*der_mixed(5)...
                ...
                +mx_P(ialpha,ibeta,ilocal3,6)*der_1st(4)...
                +mx_P(ialpha,ibeta,ilocal3,2)*der_2nd(4)...
                +mx_P(ialpha,ibeta,ilocal3,3)*der_mixed(2)...
                +mx_P(ialpha,ibeta,ilocal3,1)*der_mixed(6);

            derivatives_global=mx_chain\derivatives_local;

            for i=1:1:4
                sol(ilocal3,i+1,ialpha,ibeta)=derivatives_global(i);
            end

        end
        
    end
end

end

function mx=mx_constraints_assembler(twists)
% Twists are in radians
 mx=zeros(12,12);
 first_index=-2;
for node=1:1:4
    first_index=first_index+3;
    mx(first_index,first_index)=1;
    %============================FIX FOR NO TWIST CASE LATER!!!   
%     if twist(node)==0
    if 1==1
        mx_rotation=...
            [cos(twists(node)) sin(twists(node))
            -sin(twists(node)) cos(twists(node))];
        x=mx_rotation\[1;0];
        mx(first_index+1,first_index+1)=x(1);
        mx(first_index+2,first_index+1)=x(2);
        x=mx_rotation\[0;1];
        mx(first_index+1,first_index+2)=x(1);
        mx(first_index+2,first_index+2)=x(2);
    end
end
end

function mx=mx_chain_assembler(dx_dalpha,dy_dalpha,dx_dbeta,dy_dbeta,...
    d2x_dalphadbeta,d2y_dalphadbeta)
mx=zeros(5,5);
mx(1,1)=dx_dalpha;
mx(1,2)=dy_dalpha;
mx(2,1)=dx_dbeta;
mx(2,2)=dy_dbeta;
mx(3,3)=dx_dalpha^2;
mx(3,4)=dy_dalpha^2;
mx(3,5)=2*dx_dalpha*dy_dalpha;
mx(4,3)=dx_dbeta^2;
mx(4,4)=dy_dbeta^2;
mx(4,5)=2*dx_dbeta*dy_dbeta;
mx(5,1)=d2x_dalphadbeta;
mx(5,2)=d2y_dalphadbeta;
mx(5,3)=dx_dalpha*dx_dbeta;
mx(5,4)=dy_dalpha*dy_dbeta;
mx(5,5)=dx_dalpha*dy_dbeta+dx_dbeta*dy_dalpha;
end

function mx=P(alpha,beta,mx_poly_coef)
mx=zeros(12,6);
for i=1:1:12
    mx(i,1)=mx_poly_coef(i,1)+mx_poly_coef(i,2)*alpha+mx_poly_coef(i,3)*beta+mx_poly_coef(i,4)*alpha^2+mx_poly_coef(i,5)*(beta^2)...
        +mx_poly_coef(i,6)*alpha*beta+mx_poly_coef(i,7)*alpha^3+mx_poly_coef(i,8)*beta^3+mx_poly_coef(i,9)*alpha^2*beta...
        +mx_poly_coef(i,10)*alpha*beta^2+mx_poly_coef(i,11)*alpha^3*beta+mx_poly_coef(i,12)*alpha*beta^3;
    mx(i,2)=mx_poly_coef(i,2)+2*mx_poly_coef(i,4)*alpha+mx_poly_coef(i,6)*beta+3*mx_poly_coef(i,7)*(alpha^2)...
        +2*mx_poly_coef(i,9)*alpha*beta+mx_poly_coef(i,10)*beta^2+3*mx_poly_coef(i,11)*alpha^2*beta+mx_poly_coef(i,12)*beta^3;
    mx(i,3)=mx_poly_coef(i,3)+2*mx_poly_coef(i,5)*beta+mx_poly_coef(i,6)*alpha+3*mx_poly_coef(i,8)*(beta^2)...
        +mx_poly_coef(i,9)*alpha^2+2*mx_poly_coef(i,10)*alpha*beta+mx_poly_coef(i,11)*alpha^3+3*mx_poly_coef(i,12)*alpha*beta^2;
    mx(i,4)=2*mx_poly_coef(i,4)+6*mx_poly_coef(i,7)*alpha+2*mx_poly_coef(i,9)*beta+6*mx_poly_coef(i,11)*alpha*beta;
    mx(i,5)=2*mx_poly_coef(i,5)+6*mx_poly_coef(i,8)*beta+2*mx_poly_coef(i,10)*alpha+6*mx_poly_coef(i,12)*alpha*beta;
    mx(i,6)=mx_poly_coef(i,6)+2*mx_poly_coef(i,9)*alpha+2*mx_poly_coef(i,10)*beta+3*mx_poly_coef(i,11)*alpha^2+3*mx_poly_coef(i,12)*beta^2;
end
end

function sol=mapper_dx_or_y_dalpha(beta,x_or_y)
sol=-x_or_y(1)*(1-beta)+x_or_y(2)*(1-beta)-x_or_y(3)*beta+x_or_y(4)*beta;
end

function sol=mapper_dx_or_y_dbeta(alpha,x_or_y)
sol=-x_or_y(1)*(1-alpha)-x_or_y(2)*alpha+x_or_y(3)*(1-alpha)+x_or_y(4)*alpha;
end

function sol=mapper_d2x_or_y_dalphadbeta(x_or_y)
sol=x_or_y(1)-x_or_y(2)-x_or_y(3)+x_or_y(4);
end

function mx_poly_coef=generate_poly_coef(mx_coef,mx_constraints)
mx_poly_coef=zeros(12,12);
for m=1:1:12
    mx_poly_coef(:,m)=mx_coef\mx_constraints(:,m);
end
mx_poly_coef=mx_poly_coef';
end

function mx=generate_matrix_serendipity(coords)
mx=zeros(12,12);
m=0;

for node=1:1:4
    for der_compt=0:1:2
        m=m+1;
        mx(m,:)=Hermitian_poly_coef_serendipity(coords(node,1),coords(node,2),der_compt);
    end
end

end

function coef=Hermitian_poly_coef_serendipity(alpha,beta,der_compt)
if der_compt==0
    coef=[1,alpha,beta,alpha^2,beta^2,alpha*beta,alpha^3,beta^3,...
        alpha^2*beta,alpha*beta^2,alpha^3*beta,alpha*beta^3];
elseif der_compt==1
    coef=[0,1,0,2*alpha,0,beta,3*alpha^2,0,2*alpha*beta,beta^2,...
        3*alpha^2*beta,beta^3];
elseif der_compt==2
    coef=[0,0,1,0,2*beta,alpha,0,3*beta^2,alpha^2,2*alpha*beta,...
        alpha^3,3*alpha*beta^2];
elseif der_compt==3
    coef=[0,0,0,2,0,0,6*alpha,0,2*beta,0,...
        6*alpha*beta,0];
elseif der_compt==4
    coef=[0,0,0,0,2,0,0,6*beta,0,2*alpha,...
        0,6*alpha*beta];
elseif der_compt==5
    coef=[0,0,0,0,0,0,6,0,0,0,...
        6*beta,0];
end
end

function sol=compute_laplacian(alpha,beta,mx_poly_coef,nodal_weights)
d2P_dx2=Hermitian_poly_coef_serendipity(alpha,beta,3);
d2P_dy2=Hermitian_poly_coef_serendipity(alpha,beta,4);
summation=0;
for ilocal=1:1:12
    summation=summation...
        +nodal_weights(ilocal)...
        *(dot(mx_poly_coef(ilocal,:),d2P_dx2)...
        +dot(mx_poly_coef(ilocal,:),d2P_dy2));
end
sol=summation;
end

function sol=get_conc_test(alpha,beta,mx_poly_coef,nodal_weights,der_compt)
P=Hermitian_poly_coef_serendipity(alpha,beta,der_compt);
% d2P_dy2=Hermitian_poly_coef_serendipity(alpha,beta,4);
summation=0;
for ilocal=1:1:12
    summation=summation...
        +nodal_weights(ilocal)...
        *(dot(mx_poly_coef(ilocal,:),P));
end
sol=summation;
end

function sol=generateWeights_test(alpha,beta,der_compt,nodal_weights)
orientations=...
    [0 0
    0 1
    1 0
    1 1];
types=...
    [0 0
    1 0
    0 1
    1 1];
orders=...
    [0 0
    1 0
    0 1
    2 0
    0 2];
sol=0;
index=0;
for iorientation=1:1:4
    for itype=1:1:4
        index=index+1;
        sol=sol+nodal_weights(index)...
            *basis(alpha,orientations(iorientation,1),types(itype,1),orders(der_compt,1))...
            *basis(beta,orientations(iorientation,2),types(itype,2),orders(der_compt,2));
    end
end
end

function sol=generateWeights_test_2(alpha,beta,der_compt,iorientation,itype)
orientations=...
    [0 0
    0 1
    1 0
    1 1];
types=...
    [0 0
    1 0
    0 1
    1 1];
orders=...
    [0 0
    1 0
    0 1
    2 0
    0 2];
sol=0;

sol=sol...
    +basis(alpha,orientations(iorientation,1),types(itype,1),orders(der_compt,1))...
    *basis(beta,orientations(iorientation,2),types(itype,2),orders(der_compt,2));

end

function sol=getting_odd_poly_coef(gp)

sol=zeros(12,12);

mx_var=zeros(24,12);
alphas=[0 1];
betas=[0 1];
der_compts=[0 1 2];
ilocal=0;
for ibeta=1:1:2
    for ialpha=1:1:2
        for ider_compt=1:1:3
            ilocal=ilocal+1;
            mx_var(ilocal,:)=Hermitian_poly_coef_serendipity(...
                alphas(ialpha),betas(ibeta),der_compts(ider_compt));
        end
    end
end

for ialpha=1:1:2
    for igp=1:1:3
        ilocal=ilocal+1;
        mx_var(ilocal,:)=Hermitian_poly_coef_serendipity(...
            alphas(ialpha),gp(igp),1);
    end
end

for ibeta=1:1:2
    for igp=1:1:3
        ilocal=ilocal+1;
        mx_var(ilocal,:)=Hermitian_poly_coef_serendipity(...
            gp(igp),betas(ibeta),2);
    end
end

for ilocal=1:1:12
    constraints=zeros(24,1);
    constraints(ilocal)=1;
    xs=mx_var\constraints;
    sol(ilocal,:)=xs';
end
end

%======= 2019 ===============================

function [der_1st,der_2nd,der_3rd,der_mixed]=twisted_chain(cos_,sin_,...
    dx_dalpha,dy_dalpha,dx_dbeta,dy_dbeta,d2x_dalphadbeta,d2y_dalphadbeta)

der_1st=zeros(1,4);
der_2nd=zeros(1,4);
der_3rd=zeros(1,4);
der_mixed=zeros(1,6);

der_1st(1)=cos_*dx_dalpha+sin_*dy_dalpha;
der_1st(2)=cos_*dx_dbeta+sin_*dy_dbeta;
der_1st(3)=-sin_*dx_dalpha+cos_*dy_dalpha;
der_1st(4)=-sin_*dx_dbeta+cos_*dy_dbeta;

der_mixed(1)=cos_*d2x_dalphadbeta+sin_*d2y_dalphadbeta;
der_mixed(2)=-sin_*d2x_dalphadbeta+cos_*d2y_dalphadbeta;

end

function twists=twist_radially(x,y)
twists=zeros(1,4);
for node=1:1:4
    v=[x(node),y(node)];
    if sum(abs(v))==0
        continue
    else
        twists(node)=get_twist(v);
    end
end
end

function output_path=folder_setUp_curved(nex,ney,rx_to_ry_ratio,diff,...
    ci,T,n1,n2)

output_path=1;
while 1
    if ~exist([pwd '/Results_ellipse/Run_' num2str(output_path) '/'],'dir')
        break
    end
    output_path=output_path+1;
end
output_path=[pwd '/Results_ellipse/Run_' num2str(output_path) '/'];
mkdir(output_path);

mkdir([output_path 'specification/'])
mkdir([output_path 'time/'])
mkdir([output_path 'concentration/'])

dlmwrite([output_path 'specification/nex_ney.dat'],[nex ney]);
dlmwrite([output_path 'specification/rx_to_ry_ratio.dat'],rx_to_ry_ratio);

dlmwrite([output_path 'specification/diff.dat'],diff);
dlmwrite([output_path 'specification/ci.dat'],ci);
dlmwrite([output_path 'specification/T.dat'],T);

dlmwrite([output_path 'specification/n1_n2.dat'],[n1 n2]);

end

function c_new=ave_c_adjuster(c,ci,weights,ne,nx,nex,w,determinants,nthree)
%ONLY WORKS ASSUMING THE DOMAIN HAS AN AREA OF ONE
c_new=c;
conc_=get_conc_first_order_serendipity(c,weights,ne,nx,nex);
summation=0;
for e=1:1:ne
    for ix=1:1:3
        for iy=1:1:3
            summation=summation+w(ix)*w(iy)*determinants(e,ix,iy)...
                *conc_(e,ix,iy);
        end
    end
end

deviation=summation-ci;

for i=1:3:nthree-2
    c_new(i)=c(i)-deviation;
end

%--------
% c=c_new;
% conc_=get_conc_first_order_serendipity(c,weights,ne,nx,nex);
% summation=0;
% for e=1:1:ne
%     for ix=1:1:3
%         for iy=1:1:3
%             summation=summation+w(ix)*w(iy)*determinants(e,ix,iy)...
%                 *conc_(e,ix,iy);
%         end
%     end
% end
% 
% summation

end

function [sol1,sol2]=compute_weights_serendipity_assist(keys,load,load_more,...
    nx,nex,xy,gp,mx_P)
sol1=zeros(load_more,12,5,3,3);
sol2=zeros(load_more,3,3);
x=zeros(1,4);y=zeros(1,4);
for ikey=1:1:load
    nodes=get_eNodes_serendipity(keys(ikey),nx,nex);
    for inode=1:1:4
        x(inode)=xy(nodes(inode),1);
        y(inode)=xy(nodes(inode),2);
    end
    %=====ONLY WORKS FOR CIRCLE FOR NOW===
    twists=twist_radially(x,y);
    %==================================
    sol1(ikey,:,:,:,:)=get_weight_global(gp,x,y,twists,mx_P);
    for iy=1:1:3
        for ix=1:1:3
            sol2(ikey,ix,iy)=det([...
                mapper_dx_or_y_dalpha(gp(iy),x),mapper_dx_or_y_dbeta(gp(ix),x)...
                ;mapper_dx_or_y_dalpha(gp(iy),y),mapper_dx_or_y_dbeta(gp(ix),y)]);
        end
    end
end
end
