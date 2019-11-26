function main_test

clear
clc

% === Start the timer =====================================================

counter_init=tic;

% === User defined variables ==============================================

% --- Define mesh ---------------------------------------------------------

nex=5;
ney=4;
nez=3;

shape=2;

% So far, we have the following available
% 1 - Cube
% 2- Sphere

% The data must have been made ready with the following file

% main_isoparametric_mapping_3d

% --- Polymer properties --------------------------------------------------

n1=1;
n2=10;
 
diff=1800;

%======================
ci=0.6;
T=[0.3];
%======================

include_fluc=1;
ci_fluc = 0.0001;

entropy=1;
T_theta=1;

ci_critical=0; % This bypasses ci specified above
T_critical=0; % This bypasses T specified above

% --- Simulation properties -----------------------------------------------

parallel_computing=1; % Non parallel computing uses traditional method

show_figure=0;
export_data=0; % Yes=1, No=0


dt=5.0*10^-7;

obs_t=2.0e-1;

nr_tol=1.0e-6;
nr_max_iteration=20;
%======================
time_out=3600*2.5;
max_frame=1000;
%======================
fps=10;
dt_down=0.5;
dt_up=1.2;
dt_min=1.0e-150;
dt_ideal=5.0e-10;
bypass_tol=5;

tol=1.0e-6;

predict_with_acceleration=1;

commit_sfsj_terms_to_RAM=1;

% --- Phase diagram properties --------------------------------------------

nT=100;
nInitialGuessPts=100;
guess_accuracy=10^-3;

proceed_request=1;

ifig_phase=1;

% --- Analysis specifications ---------------------------------------------

alpha=200;
beta=8;

% --- Illustration specifications -----------------------------------------

frame_pos=[0 0];
frame_size=[720 720];
ifig_simulation=2;

face_colors={'yellow','blue','red'};
% The order is {Mean Lower Higher}
view_angle=3;
transparency=0.9;
edge_color='none';
ifig_conclusion=3;

domain_EdgeColor='green';
domain_FaceColor='white';
domain_face_transparency=0.1;
domain_edge_transparency=0.8;

% === Validate inputs =====================================================

% VALIDATE THEM HERE

if parallel_computing~=1 || commit_sfsj_terms_to_RAM~=1
    error('Only one method is available. Set parallel_computing and commit_sfsj_terms_to_RAM to 1!')
end

% return

% === Set up variables ====================================================

[nx,ny,nz,n,neight]=setUp_3d(nex,ney,nez);
if ci_critical==1 && T_critical==1
    [ci,T]=identify_critical_point(n1,n2,entropy,T_theta);
elseif ci_critical==1
    [ci,~]=identify_critical_point(n1,n2,entropy,T_theta);
elseif T_critical==1
    [~,T]=identify_critical_point(n1,n2,entropy,T_theta);
end

c_obs=[mean(ci) mean(ci)-ci_fluc mean(ci)+ci_fluc];
ne=nex*ney*nez;
nexney=nex*ney;
dto=inf;
if predict_with_acceleration==1
    dtoo=inf;
end
bypass=0;



if shape==1
    nnz_=nnz_sj_3d(nx,ny,nz);
    [gbfs1,gbfs2,iZeros,iOnes]=sj_mapper_3d(nx,ny,nz);
    fprintf('Importing Data...\n')
    [X,Y,Z,weights]=...
        retrieve_isoparametrically_mapped_weights_3d(...
        shape,nex,ney,nez,parallel_computing);
    fprintf('Computing determinants...\n')
    determinants=isoparametric_determinant_3d(...
        ne,nex,nx,ny,nexney,X,Y,Z,parallel_computing);
    wTerms=generate_wTerms_curved_3d(ne,weights);
elseif shape==2
    [nxny,nx_,ny_,n_nSurf,n_eSurf,neTotal,nTotal,nnz_,neightTotal]=setUp_sph(...
        n,nx,ny,nz,ne,nex,ney,nez);
    [gbfs1,gbfs2]=sj_mapper_sph(n,nx,ny,nz,nnz_,nxny,n_nSurf,nx_,nTotal);
    fprintf('Importing Data...\n')
    [X,Y,Z,weights]=...
        retrieve_isoparametrically_mapped_weights_3d(...
        shape,nex,ney,nez,parallel_computing);
    fprintf('Finished importing Data...\n')
    
    %-------TEST-------------------------
    
    rotation_mx=rotate_sph(n_nSurf,n,X,Y,Z);
    inverse_rotation=inversely_rotate_sph(n_nSurf,n,X,Y,Z);
    
%     test68(n_nSurf,rotation_mx,inverse_rotation)
%     
%     disp('afdasdgag')
%     return
    
    bc_refiner=get_surface_normal_sph(...
        ne,nex,ney,nexney,nx,ny,nz,n,nxny,n_eSurf,n_nSurf,nx_,...
        X,Y,Z,parallel_computing);
    
    normals_=get_nodal_normals(n,n_nSurf,X,Y,Z);
    
    [bc_row_refiner,bc_col_refiner,bc_refiner]=...
    gbfs1_gbfs2_with_bc_sphere(neightTotal,neight,n_nSurf,normals_);
%     disp('hey yo')
%     return
    
    


    
    %-------------------------------------
    
    determinants=isoparametric_determinant_sph(...
        ne,neTotal,nex,ney,n,nx,ny,nz,nxny,nexney,X,Y,Z,...
        n_eSurf,n_nSurf,nx_,parallel_computing);
    wTerms=generate_wTerms_curved_3d(neTotal,weights);
    %---Test-------------------------
%     nnz_=nnz_sj_3d(nx,ny,nz);
%     [gbfs1,gbfs2,iZeros,iOnes]=sj_mapper_3d(nx,ny,nz);
%     nTotal=n;
%     neightTotal=neight;
%     neTotal=ne;
    disp('Surface test start')
    take=dlmread('test_weights_surface.dat');
    
    sf=size(take);
    if sf(1)~=n_eSurf(6)
        error('bad')
    end
    
    weights_surface=zeros(n_eSurf(6),8,8,7,3,3);
    disp('reformatting')
    toc(counter_init)
    for e=1:1:n_eSurf(6)
        index=0;
        eInfo=take(e,:);
        for iorientation=1:1:8
            for itype=1:1:8
                for iorder=1:1:7
                    for idim1=1:1:3
                        for idim2=1:1:3
                            index=index+1;
                            weights_surface(e,iorientation,itype,iorder,idim1,idim2)=...
                                eInfo(index);
                        end
                    end
                end
            end
        end
    end
    wTerms_surface=generate_wTerms_curved_surface_sph(n_eSurf,weights_surface);
    determinants_surface=isoparametric_determinant_surface_sph(...
        ne,nex,ney,n,nx,ny,nz,nxny,nexney,X,Y,Z,n_eSurf,n_nSurf,...
        nx_,parallel_computing);
    
    
    
%     weights_surface=weights_surface(:,:,:,2:1:4,:,:);
    
    %---TEST----------------------------------------
    
%     normals_(1,:)
%     
%     take=[rotation_mx(1,1,1) rotation_mx(1,1,2) rotation_mx(1,1,3)];
%     
%     
%     convert_to_unit_vector(take)
%     
%     
%     return
    
    %----------------------------------------------
    
    
    
    normals=get_surface_normal_sphere(n_eSurf,ne,nex,ney,nexney,...
        nx,ny,nz,n,nxny,n_nSurf,nx_,X,Y,Z);
    
%     size(normals)
%     
%     
%     
%     nodes=get_nodes_of_element_sph(6+ne,ne,nex,ney,nexney,nx,ny,nz,n,nxny,n_eSurf,n_nSurf,nx_);
%     
%     a=normals_(nodes(1)-n,:)
%     b=normals_(nodes(3)-n,:)
%     c=normals_(nodes(5)-n,:)
%     d=normals_(nodes(7)-n,:)
%     take=[normals(6,2,2,1) normals(6,2,2,2) normals(6,2,2,3)]
%     hold on
%     plot3([a(1) 0],[a(2) 0],[a(3) 0])
%     plot3([b(1) 0],[b(2) 0],[b(3) 0])
%     plot3([c(1) 0],[c(2) 0],[c(3) 0])
%     plot3([d(1) 0],[d(2) 0],[d(3) 0])
%     plot3([take(1) 0],[take(2) 0],[take(3) 0],'r*')
%     hold off
%     grid on
%     xlabel('x');ylabel('y');zlabel('z')
%     
%     
%     
%     dot(take,normals_(nodes(1),:))
%     
%     disp('fdasdfds')
%     return
%     hold on
%     for ijk=1:1:n_eSurf(6)
%         for dim1=1:1:3
%             for dim2=1:1:3
%                 xxx=normals(ijk,dim1,dim2,1);
%                 yyy=normals(ijk,dim1,dim2,2);
%                 zzz=normals(ijk,dim1,dim2,3);
%                 plot3([xxx],[yyy],[zzz],'ro')
%             end
%         end
%     end
%     hold 
%     grid on
%     xlabel('x')
%     ylabel('y')
%     zlabel('z')
%     size(normals)
%     
%     return
    %--------------------------------
    
    
    
    
    
end


toc(counter_init)
fig_domain=[min(X) max(X) min(Y) max(Y) min(Z) max(Z)];
% THIS FUNCTION CAN BE PUT IN DIFFERENT MAIN FILE
toc(counter_init)



if show_figure==1
    [x_minus_x,x_minus_y,x_minus_z,x_plus_x,x_plus_y,x_plus_z,...
        y_minus_x,y_minus_y,y_minus_z,y_plus_x,y_plus_y,y_plus_z,...
        z_minus_x,z_minus_y,z_minus_z,z_plus_x,z_plus_y,z_plus_z,...
        x_plot,y_plot,z_plot]=...
        setUp_illustration_curved_structure_3d(nx,ny,nz,X,Y,Z);
end

% --- Break down temperature gradient -------------------------------------

if shape==1
    [diffT,chi]=create_temperature_gradient_3d(...
        parallel_computing,shape,T,diff,entropy,T_theta,nx,ny,nexney,ne,nex,X,Y,Z);
elseif shape==2
    %MUST BE FIXED LATER!!!
    [diffT,chi]=create_temperature_gradient_3d(...
        parallel_computing,shape,T,diff,entropy,T_theta,nx,ny,nexney,neTotal,nex,X,Y,Z);
end

% --- Set up folder for outputting results --------------------------------

if export_data==1 || show_figure==1
    output_path=folder_setUp_3d(export_data);
    
    if export_data==1
        sf=[nex;ney;nez;dt;obs_t;diff;include_fluc;ci_fluc;nr_tol;...
            nr_max_iteration;time_out;max_frame;entropy;T_theta;...
            dt_down;dt_up;dt_min;dt_ideal;bypass_tol;n1;n2];
        dlmwrite([output_path 'specification/spec_general'],sf);
        dlmwrite([output_path 'specification/status'],0);
        dlmwrite([output_path 'specification/spec_T'],T);
        dlmwrite([output_path 'specification/spec_ci'],ci);
    end

% === Phase diagram =======================================================

    fprintf('Generating phase diagram...\n')
    [xspinodal,yspinodal,xbinodal,ybinodal]=spinodal_binodal(...
        n1,n2,entropy,T_theta,nT,tol,parallel_computing,nInitialGuessPts,...
        guess_accuracy);
	if show_figure==1
        figure(ifig_phase)
        illustrate_simulation(xspinodal,yspinodal,xbinodal,ybinodal,T,ci)
        if proceed_request==1
            proceed=proceed_dialog();
            if proceed~=1
                return
            end
        end
	end
end

% === Generate initial condition ==========================================

%THIS FUNCTION IS NOT DESGINED FOR RADIAL GRADIENT!!!   
if shape==1
    co=generate_co_3d(ci,include_fluc,ci_fluc,n,neight,nx,ny,nz);
elseif shape==2
    co=generate_co_3d(ci,include_fluc,ci_fluc,nTotal,neightTotal,nx,ny,nz);
end

% === Initialize the first time step ======================================

time=0.0;
c=co;

if predict_with_acceleration==1
    coo=co;
end

iframe=1;

% --- Plot the first time step --- %

if show_figure==1
    fig=figure(ifig_simulation);
    set(fig,'Position',[frame_pos(1) frame_pos(2) frame_size(1) frame_size(2)])
    clf(fig)
    visual(iframe)=illustrate_and_analyze_curved_structure_3d(...
        c,c_obs,fig_domain,face_colors,edge_color,...
        nx,ny,nz,view_angle,transparency,time,...
        x_minus_x,x_minus_y,x_minus_z,x_plus_x,x_plus_y,x_plus_z,...
        y_minus_x,y_minus_y,y_minus_z,y_plus_x,y_plus_y,y_plus_z,...
        z_minus_x,z_minus_y,z_minus_z,z_plus_x,z_plus_y,z_plus_z,...
        x_plot,y_plot,z_plot,domain_EdgeColor,domain_FaceColor,...
        domain_face_transparency,domain_edge_transparency);
end

% --- Output the first time step --- %

if export_data==1
    dlmwrite([output_path 'time/iteration_1'],0);
    dlmwrite([output_path 'concentration/iteration_1'],c);
end

%------Test----------------------------------------

test1=2;
[node,type]=analyze_gbs_3d(neight+(test1-1)*8+2)
take=c(neight+(test1-1)*8+2:1:neight+(test1-1)*8+4)

take=[]

% disp('finch station')
% return

%-------------------------------------------------

% ==== Enter main loop ====================================================

fprintf('Enter main cycle\n')
while 1
    
% === Update concentration profile ========================================

    if shape==1
        if predict_with_acceleration==1
            cp=predict_c_with_acceleration_3d(c,co,coo,dt,dto,dtoo,...
                parallel_computing,neight);
            cooo=coo;
        elseif predict_with_acceleration==0
            cp=predict_c_3d(c,co,dt,dto,parallel_computing,neight);
        end
    elseif shape==2
        if predict_with_acceleration==1
            cp=predict_c_with_acceleration_3d(c,co,coo,dt,dto,dtoo,...
                parallel_computing,neightTotal);
            cooo=coo;
        elseif predict_with_acceleration==0
            cp=predict_c_3d(c,co,dt,dto,parallel_computing,neightTotal);
        end
    end
    
    coo=co;
    co=c;
    c=cp;
    
% === Rotate first derivative of surface ==================================
%     
%     if shape==2
%         fprintf('Rotating the first derivatives\n')
%         cr=rotate_spherical2Cartesian(c,neight,n_nSurf,inverse_rotation);
%         cor=rotate_spherical2Cartesian(co,neight,n_nSurf,inverse_rotation);
%         
%         
% %         crr=rotate_spherical2Cartesian(cr,neight,n_nSurf,rotation_mx);
% %         
% %         sum(abs(crr-cr))
% %         
% %         disp('systematic')
% %         return
%     end
    
%     for i=1:1:neightTotal
%         if c(i)~=cr(i)
%             [node,type]=analyze_gbs_3d(i);
%             if node<=n
%                 error('bad')
%             end
%             if type~=2 && type~=3 && type~=4
%                 error('bad')
%             else
%                 type
%             end
%             
%         end
%     end

% === Previous concentration (Only for parallel algorithm) ================


    if shape==1
        conco_=get_conc_type_I_curved_3d(ne,co,weights,nexney,nex,nx,ny);
    elseif shape==2
        conco_=get_conc_type_I_curved_sph(neTotal,co,weights,ne,...
                nexney,nex,ney,nx,ny,nz,nxny,n,n_eSurf,n_nSurf,nx_,parallel_computing);
    end
    
% === Apply BC to newly predicted concentration ===========================

    % This should be unnecesssary as terms that need to be zeroed should
    % have already been zeroed.
    
% === Enter Newton-Raphson iterations =====================================

    jk=0;
    err=inf;
    fprintf('Entering Newton Raphson iterations\n')
    
	while err>nr_tol
        
        toc(counter_init)
        
% === Assert reasonable iteration count ===================================

        jk=jk+1;
        if jk>nr_max_iteration
            disp('Too many Newton Raphson iterations')
            jk=-1;
            break
        end
        
% === Compute residual vector and Jacobian matrix =========================

% === Compute using parallelization =======================================

        if shape==1
            fprintf('Computing conc...\n')
            sj=get_conc_curved_3d(ne,c,weights,nexney,nex,nx,ny);
            fprintf('Computing sf and sj terms...\n')
            [terms,sf,sj]=compute_terms_curved_3d(ne,sj,chi,n1,n2,nex);
            fprintf('Computing residual vector...\n')
            sf=compute_sf_with_terms_curved_3d(neight,sf,conco_,nx,ny,nz,nex,ney,weights,dt,determinants,diffT,wTerms,terms);
            fprintf('Computing Jacobian matrix...\n')
            sj=compute_sj_with_terms_curved_3d(gbfs1,gbfs2,nnz_,dt,weights,nex,ney,nx,ny,nz,determinants,diffT,wTerms,terms,sj);

%             mean(sf)
%             mean(sj)
%             return
            
            fprintf('Applying boundary conditions...\n')
            [sf,sj]=bc_3d(sf,sj,iZeros,iOnes,gbfs1);
            
        elseif shape==2
            fprintf('Computing conc...\n')
%             sj=get_conc_curved_3d(ne,c,weights,nexney,nex,nx,ny);

            %--------------------------------------------------------------
            sj=get_conc_curved_sph(neTotal,c,weights,ne,nexney,nex,ney,...
                nx,ny,nz,nxny,n,n_eSurf,n_nSurf,nx_);
            %--------------------------------------------------------------
            
            
            
            
            fprintf('Computing sf and sj terms...\n')
            [terms,sf,sj]=compute_terms_curved_3d(neTotal,sj,chi,n1,n2,nex);
            
            %---------------------------------------
            terms_surface=compute_terms_curved_surface_sph(...
                    ne,nex,ney,nexney,nx,ny,nz,n,nxny,n_eSurf,n_nSurf,nx_,weights_surface,c);
            
            %-----------------------------------------
            
%             fprintf('Computing residual vector...\n')
%             sf=compute_sf_with_terms_curved_3d(neight,sf,conco_,nx,ny,nz,nex,ney,weights,dt,determinants,diffT,wTerms,terms);
%             fprintf('Computing Jacobian matrix...\n')
%             sj=compute_sj_with_terms_curved_3d(gbfs1,gbfs2,nnz_,dt,weights,nex,ney,nx,ny,nz,determinants,diffT,wTerms,terms,sj);
%             fprintf('Applying boundary conditions...\n')
%             [sf,sj]=bc_3d(sf,sj,iZeros,iOnes,gbfs1);
%             
%             disp('Miku')
%             size(sf)
%             size(sj)
%             return
            
            fprintf('Computing residual vector...\n')
            
            sf=test62(...
                neightTotal,sf,conco_,nx,ny,nz,nex,ney,weights,...
                dt,determinants,diffT,wTerms,terms,...
                n,ne,n_eSurf,n_nSurf,nTotal,nx_,...
                determinants_surface,terms_surface,normals,weights_surface);
            

            fprintf('Computing Jacobian matrix...\n')
            sj=test63(...
                gbfs1,gbfs2,nnz_,dt,weights,nex,ney,nx,ny,nz,...
                determinants,diffT,wTerms,terms,...
                ne,nez,n,neight,nexney,n_eSurf,n_nSurf,nx_,sj,...
                wTerms_surface,weights_surface,determinants_surface,normals);

            
            
%             size(sj)
%             nnz(sj)
            
            
            
        end

% --- Sparsing --- %

        fprintf('Sparsing matrix...\n')
        sj=sparse(gbfs1,gbfs2,sj);
%         sj=generate_traditional_mx(gbfs1,gbfs2,sj);
        
%         first_index=neight-7;
%         for ijk=1:1:n_nSurf(6)
%             first_index=first_index+8;
%             
% %             sj(first_index+1,:)=zeros(1,neightTotal);
% %             sj(first_index+1,first_index+1:1:first_index+3)=normals_(ijk,:);
% %             sf(first_index+1)=0;
% %             
% %             sj(first_index+2,:)=zeros(1,neightTotal);
% %             sj(first_index+2,first_index+1:1:first_index+3)=normals_(ijk,:);
% %             sf(first_index+2)=0;
%             
%             sj(first_index+3,:)=zeros(1,neightTotal);
%             sj(first_index+3,first_index+1:1:first_index+3)=normals_(ijk,:);
%             sf(first_index+3)=0;
%         end
%         sjr=sj;
%         
%         sj=test66(sj,neight,n_nSurf,neightTotal,inverse_rotation);
%         
%         sj2=sjr;
%         first_index=neight-7;
%         for ijk=1:1:n_nSurf(6)
%             first_index=first_index+8;
%             
%             for col=neight+1:1:neightTotal
%                 [~,type]=analyze_gbs_3d(col);
%                 if type==2
%                     sj2(first_index+1,col)=...
%                         inverse_rotation(ijk,1,1)*sjr(first_index+1,col)...
%                         +inverse_rotation(ijk,2,1)*sjr(first_index+1,col+1)...
%                         +inverse_rotation(ijk,3,1)*sjr(first_index+1,col+2);
%                 elseif type==3
%                     sj2(first_index+1,col)=...
%                         inverse_rotation(ijk,1,2)*sjr(first_index+1,col-1)...
%                         +inverse_rotation(ijk,2,2)*sjr(first_index+1,col)...
%                         +inverse_rotation(ijk,3,2)*sjr(first_index+1,col+1);
%                 elseif type==4
%                     sj2(first_index+1,col)=...
%                         inverse_rotation(ijk,1,3)*sjr(first_index+1,col-2)...
%                         +inverse_rotation(ijk,2,3)*sjr(first_index+1,col-1)...
%                         +inverse_rotation(ijk,3,3)*sjr(first_index+1,col);
%                 end
%                 
%                 if type==2
%                     sj2(first_index+2,col)=...
%                         inverse_rotation(ijk,1,1)*sjr(first_index+2,col)...
%                         +inverse_rotation(ijk,2,1)*sjr(first_index+2,col+1)...
%                         +inverse_rotation(ijk,3,1)*sjr(first_index+2,col+2);
%                 elseif type==3
%                     sj2(first_index+2,col)=...
%                         inverse_rotation(ijk,1,2)*sjr(first_index+2,col-1)...
%                         +inverse_rotation(ijk,2,2)*sjr(first_index+2,col)...
%                         +inverse_rotation(ijk,3,2)*sjr(first_index+2,col+1);
%                 elseif type==4
%                     sj2(first_index+2,col)=...
%                         inverse_rotation(ijk,1,3)*sjr(first_index+2,col-2)...
%                         +inverse_rotation(ijk,2,3)*sjr(first_index+2,col-1)...
%                         +inverse_rotation(ijk,3,3)*sjr(first_index+2,col);
%                 end
%                 
%                 if type==2
%                     sj2(first_index+3,col)=...
%                         inverse_rotation(ijk,1,1)*sjr(first_index+3,col)...
%                         +inverse_rotation(ijk,2,1)*sjr(first_index+3,col+1)...
%                         +inverse_rotation(ijk,3,1)*sjr(first_index+3,col+2);
%                 elseif type==3
%                     sj2(first_index+3,col)=...
%                         inverse_rotation(ijk,1,2)*sjr(first_index+3,col-1)...
%                         +inverse_rotation(ijk,2,2)*sjr(first_index+3,col)...
%                         +inverse_rotation(ijk,3,2)*sjr(first_index+3,col+1);
%                 elseif type==4
%                     sj2(first_index+3,col)=...
%                         inverse_rotation(ijk,1,3)*sjr(first_index+3,col-2)...
%                         +inverse_rotation(ijk,2,3)*sjr(first_index+3,col-1)...
%                         +inverse_rotation(ijk,3,3)*sjr(first_index+3,col);
%                 end
%                 
%             end
%         end
%         
%         sum(sum(sj2-sj))
%         
%         sj=sj2;
%         return
        
%         for ijk=1:1:n_nSurf(6)
%             row=neight+(ijk-1)*8+2;
%             for col=1:1:neightTotal
%                 sj(row,col)=0;
%             end
%             sj(row,row)=1;
%             sf(row)=0;
%         end
        
%         sj=test66(sj,neight,n_nSurf,neightTotal,inverse_rotation);
        
%         sj=sparse([gbfs1 bc_row_refiner],[gbfs2 bc_col_refiner],[sj bc_refiner]);
%         sj=sparse(...
%             [gbfs1 bc_row_refiner 1],[gbfs2 bc_col_refiner neightTotal+n_nSurf(6)],[sj bc_refiner 0]);
%         sf=[sf zeros(1,n_nSurf(6))];
        %----Test--------------------------------------------
%         
%         sj=generate_traditional_mx([gbfs1 bc_row_refiner 1],[gbfs2 bc_col_refiner neightTotal+n_nSurf(6)],[sj bc_refiner 0]);
%         
%         size_=size(sj);
%         
%         for ijk=neightTotal+1:1:neightTotal+n_nSurf(6)
%             for ijk2=1:1:size_(1)
%                 sj(ijk,ijk2)=1;
%             end
%         end
%         
%         sj
%         det(sj)
        
%---------------------------------------
%         sf=[sf zeros(1,n_nSurf(6))];
%         
% %         sj=generate_traditional_mx([gbfs1 bc_row_refiner],[gbfs2 bc_col_refiner],[sj bc_refiner]);
%         sj=generate_traditional_mx(gbfs1,gbfs2,sj);
%         
%         size(sj)
%         neightTotal
%         n_nSurf(6)+neightTotal
%         
%         rank(sj)
%         disp('Systematic love')
%         return
        %------------------------------------
%         
%         max(bc_row_refiner)
%         
%         size(sj)
%         neightTotal
%         n_nSurf(6)
%         
%         return
%         
%         sjt=zeros(neightTotal,neightTotal);
%         for ijk=1:1:nnz_
%             sjt(gbfs1(ijk),gbfs2(ijk))=sj(gbfs1(ijk),gbfs2(ijk));
%         end
%         
%         sj=test66(sjt,neight,n_nSurf,neightTotal,inverse_rotation);
%         
%         
% %         disp('dafafdfasdf')
% %         return
        
        %------------------------------------------------------
        
%         size(sj)
%             disp('no sir')
%             return
        
% === Carry out matrix division ===========================================

        fprintf('Carrying out matrix division...\n')
        c_=sj\-sf';
        
        
%         length(c_)
%         
%         
%         sf(length(sf))
%         disp('da hell')
%         %Test---------------------------------------
%         
%         
%         test=51;
%         
%         take=neight+(test-1)*8+2:1:neight+(test-1)*8+4;
%         
%         taken=c_(take);
%         
%         
%         normal=normals_(test,:);
%         
%         taken
%         
%         dot(taken,normal)
%         
%         disp('franky')
%         
%         normal
%         
%         rotation_mx(test,1,:)
%         
%         size(normals_)
%         
%         return
%         
% %         size(rotation_mx)
%         sj(neightTotal+test,neight+(test-1)*8+2)
%         sj(neightTotal+test,neight+(test-1)*8+3)
%         sj(neightTotal+test,neight+(test-1)*8+4)
%         
%         i1=neight+(test-1)*8+2;
%         i2=neight+(test-1)*8+3;
%         i3=neight+(test-1)*8+4;
%         
%         for ijk=1:1:neightTotal
%             if sj(neightTotal+test,ijk)~=0
%                 if ijk~=i1 && ijk~=i2 && ijk~=i3
%                     error('fdafasdf')
%                 end
%             end
%         end
%         
%         
%         tester1=-sj*c_;
%         
%         6
%         tester1
%         
%         tester1
%         
%         c_(length(c_))
%         
%         tester1(length(tester1))
%         
%         return
        
        
%         len=length(c_);
%         take=c_(len-7:1:len);
%         
%         take';
%         
%         
%         rotation_mx(n_nSurf(6),1,1);
%         rotation_mx(n_nSurf(6),1,2);
%         rotation_mx(n_nSurf(6),1,3);
%         
%         take2=[rotation_mx(n_nSurf(6),1,1) rotation_mx(n_nSurf(6),1,2) rotation_mx(n_nSurf(6),1,3)];
%         
%         take(2)*rotation_mx(n_nSurf(6),1,1)...
%             +take(3)*rotation_mx(n_nSurf(6),1,2)...
%             +take(4)*rotation_mx(n_nSurf(6),1,3);
%         
%         
%         dir=[X(nTotal) Y(nTotal) Z(nTotal)];
%         
%         totarian=sqrt(sum(dir.^2));
%         
%         disp('dfaasfdadgags')
%         return
        
        %--------------------------------------------
        
        
        
           
% === Update the solution =================================================
        
%         length(sf)
%         size(c)
%         size(c_)
%         
%         neightTotal
%         size(sj)
        c_=c_(1:1:neightTotal);
%         return

        
        fprintf('Assessing current iterations...\n')
        c=c+c_';
        
        %---Test---------------
        live=zeros(1,3);
        disp('-----------------------------')
        node=51;
        gbfs=neight+(node-1)*8+2:1:neight+(node-1)*8+4;
        c(gbfs)
        dot(normals_(node,:),c(gbfs))
        if abs(dot(normals_(node,:),c(gbfs)))<min(abs(c(gbfs)))
            live(1)=1;
        end
        disp('-----------------------------')
        node=1;
        gbfs=neight+(node-1)*8+2:1:neight+(node-1)*8+4;
        c(gbfs)
        dot(normals_(node,:),c(gbfs))
        if abs(dot(normals_(node,:),c(gbfs)))<min(abs(c(gbfs)))
            live(2)=1;
        end
        disp('-----------------------------')
        node=round(n_nSurf(6)*0.75);
        gbfs=neight+(node-1)*8+2:1:neight+(node-1)*8+4;
        c(gbfs)
        dot(normals_(node,:),c(gbfs))
        if abs(dot(normals_(node,:),c(gbfs)))<min(abs(c(gbfs)))
            live(3)=1;
        end
        disp('-----------------------------')
        live
%         
        
        %--------------------
        
% === Evaluate error ======================================================

        err=sqrt(sum(c_.^2.0))
        
% === Assert convergence ==================================================

        if jk>1 && err>=erro
            disp('Newton Raphson is diverging!')
            jk=-1;
            break
        else
            erro=err;
        end
	end
    
% === Exit Newton Raphson iteration =======================================

	fprintf('Exit Newton-Raphson iterations\n')
    
% === Resolve divergence ==================================================

    if jk==-1
        bypass=bypass+1;
        if bypass>bypass_tol
            disp('Too many bypassing!')
            if export_data==1
                termination_type=4;
            end
            break
        else
            dt=dt*dt_down;
            if dt<dt_min
                disp('Time step is too small!')
                if export_data==1
                    termination_type=5;
                end
                break
            else
                c=co;
                co=coo;
                if predict_with_acceleration==1
                    coo=cooo;
                end
                continue
            end
        end
    elseif bypass~=0
        bypass=0;
    end

% === Update time for most recently obtained result =======================

    time=time+dt;
    iframe=iframe+1;

% === Exprt the data if specified =========================================

    if export_data==1
        dlmwrite([output_path 'time/iteration_' num2str(iframe)],time);
        dlmwrite([output_path 'concentration/iteration_' num2str(iframe)],c);
    end

% === Analyze the structure factor ========================================

    if show_figure==1

% === PLOT ================================================================
    
        clf(fig)
        visual(iframe)=illustrate_and_analyze_curved_structure_3d(...
            c,c_obs,fig_domain,face_colors,edge_color,...
            nx,ny,nz,view_angle,transparency,time,...
            x_minus_x,x_minus_y,x_minus_z,x_plus_x,x_plus_y,x_plus_z,...
            y_minus_x,y_minus_y,y_minus_z,y_plus_x,y_plus_y,y_plus_z,...
            z_minus_x,z_minus_y,z_minus_z,z_plus_x,z_plus_y,z_plus_z,...
            x_plot,y_plot,z_plot,domain_EdgeColor,domain_FaceColor,...
            domain_face_transparency,domain_edge_transparency);
        
    end

% === Determine whether to end simulation =================================

    if iframe==max_frame
        disp('Maximum frame number specified has been reached!')
        if export_data==1
            termination_type=1;
        end
        break
    elseif time_out<toc(counter_init)
        disp('Time out!')
        if export_data==1
            termination_type=2;
        end
        break
    elseif time>obs_t
        disp('Observation time goal has been achieved!')
        if export_data==1
            termination_type=3;
        end
        break
    end
    
% === Prepare for next time step ==========================================

    if predict_with_acceleration==1
        dtoo=dto;
    end
    dto=dt;
    if jk>=10
        dt=dt*dt_down;
        if dt<dt_min
            disp('Time step is too small!')
            if export_data==1
                termination_type=5;
            end
            break
        end
    elseif jk<=3 || dt<dt_ideal
        dt=dt*dt_up;
    end
    if time+dt>obs_t
        dt=obs_t-time;
    end
    
end

% === Exit main loop ======================================================

fprintf('Exit main loop\n')
toc(counter_init)

if export_data==1
    if predict_with_acceleration==1
        steps=[dt;dto;dtoo];
    elseif predict_with_acceleration==0
        steps=[dt;dto];
    end
    dlmwrite([output_path 'specification/time_steps'],steps);
    dlmwrite([output_path 'specification/elapsed'],toc(counter_init));
    dlmwrite([output_path 'specification/termination'],termination_type);
    dlmwrite([output_path 'specification/status'],1);
    dlmwrite([output_path 'specification/shape'],shape);
end

% === Export the movie ====================================================

if show_figure==1
    try
        if length(ci)==1 && length(T)==1
            vid_title=sprintf('3D nonlinear Cahn Hilliard,%dx%dx%d,ci=%1.3f,fluc=%1.4f,diff=%d,T=%0.3f,elapsed time=%1.3fs.avi',nex,ney,nez,ci,ci_fluc*include_fluc,diff,T,toc(counter_init));
        elseif length(ci)==2 && length(T)==1
            vid_title=sprintf('3D nonlinear Cahn Hilliard,%dx%dx%d,ci=[%1.2f,%1.2f],fluc=%1.4f,diff=%d,T=%0.3f,elapsed time=%1.3fs.avi',nex,ney,nez,ci(1),ci(2),ci_fluc*include_fluc,diff,T,toc(counter_init));
        elseif length(ci)==1 && length(T)==2
            vid_title=sprintf('3D nonlinear Cahn Hilliard,%dx%dx%d,ci=%1.3f,fluc=%1.4f,diff=%d,T=[%1.2f,%1.2f],elapsed time=%1.3fs.avi',nex,ney,nez,ci,ci_fluc*include_fluc,diff,T(1),T(2),toc(counter_init));
        elseif length(ci)==2 && length(T)==2
            vid_title=sprintf('3D nonlinear Cahn Hilliard,%dx%dx%d,ci=[%1.2f,%1.2f],fluc=%1.4f,diff=%d,T=[%1.2f,%1.2f],elapsed time=%1.3fs.avi',nex,ney,nez,ci(1),ci(2),ci_fluc*include_fluc,diff,T(1),T(2),toc(counter_init));
        end
        vid_title=[output_path vid_title];
        video=VideoWriter(vid_title, 'Uncompressed AVI');
        video.FrameRate=fps;
        open(video)
        writeVideo(video,visual);
        close(video)

% === Conclusion ==========================================================

        fig=figure(ifig_conclusion);
        set(fig,'Position',[frame_pos(1) frame_pos(2) frame_size(1) frame_size(2)])
        movie(gcf,visual,1)
    catch
        warning('Failed to create the video')
    end
end

toc(counter_init)

end