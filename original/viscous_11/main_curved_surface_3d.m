function main_curved_surface_3d

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
 
diff=5000;

%======================
ci=0.45;
T=[0.4];
%======================

include_fluc=1;
ci_fluc = 0.0001;

entropy=1;
T_theta=1;

ci_critical=0; % This bypasses ci specified above
T_critical=0; % This bypasses T specified above

% --- Simulation properties -----------------------------------------------

parallel_computing=1; % Non parallel computing uses traditional method

show_figure=1;
export_data=0; % Yes=1, No=0


dt=1.0*10^-7;

obs_t=2.0e-1;

nr_tol=1.0e-6;
nr_max_iteration=20;
%======================
time_out=3600*0.5;
max_frame=1000;
%======================
fps=10;
dt_down=0.5;
dt_up=1.2;
dt_min=1.0e-150;
dt_ideal=5.0e-10;
bypass_tol=5;

tol=1.0e-6;
sparse_traditional_sj=0;

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

if shape==3
    if parallel_computing~=1 || commit_sfsj_terms_to_RAM~=1
        error('Only one method is available for shape 3')
    end
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



if shape==1 || shape==2
    nnz_=nnz_sj_3d(nx,ny,nz);
    [gbfs1,gbfs2,iZeros,iOnes]=sj_mapper_3d(nx,ny,nz);
    fprintf('Importing Data...\n')
    [X,Y,Z,weights]=...
        retrieve_isoparametrically_mapped_weights_3d(...
        shape,nex,ney,nez,parallel_computing);
    fprintf('Computing determinants...\n')
    determinants=isoparametric_determinant_3d(...
        ne,nex,nx,ny,nexney,X,Y,Z,parallel_computing);
elseif shape==3
    nxny=nx*ny;
    nx_=nx-2;
    ny_=ny-2;
    [n_nSurf,n_eSurf]=get_n_nSurf_n_eSurf(nex,ney,nez,ny,nz,nx_,ny_);
    neTotal=ne+n_eSurf(6);
    nTotal=n+n_nSurf(6);
    nnz_=nnz_sj_sph(nx,ny,nz,n_nSurf);
    [gbfs1,gbfs2]=sj_mapper_sph(n,nx,ny,nz,nnz_,nxny,n_nSurf,nx_,nTotal);
    [n_nSurf,n_eSurf]=get_n_nSurf_n_eSurf(nex,ney,nez,ny,nz,nx_,ny_);
    fprintf('Importing Data...\n')
    [X,Y,Z,weights]=...
        retrieve_isoparametrically_mapped_weights_3d(...
        shape,nex,ney,nez,parallel_computing);
    fprintf('Computing determinants...\n')
    
    take=dlmread('test_weights_surface.dat');
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
    
    nnz(wTerms_surface)

    
%     return
    
%     disp('reafdafda')
%     toc(counter_init)
%     
%     test105=isoparametrically_adjust_weights_surface_sph(...
%         ne,nex,ney,n,nx,ny,nz,nxny,nexney,X,Y,Z,parallel_computing,...
%         n_eSurf,n_nSurf,nx_);
%     
%     disp('done verification')
%     toc(counter_init)
%     
%     difference=abs(test105-weights_surface);
%     
%     sum(sum(sum(sum(sum(sum(sum(sum(sum(sum(difference))))))))))
%     
%     
%     max(max(max(max(max(max(max(max(max(max(difference))))))))))
%     
%     toc(counter_init)
    
    
    determinants_surface=isoparametric_determinant_surface_sph(...
        ne,nex,ney,n,nx,ny,nz,nxny,nexney,X,Y,Z,n_eSurf,n_nSurf,...
        nx_,parallel_computing);
    
    take=determinants_surface(:,1,1);
    nnz(take)
    size(take)
    n_eSurf(6)
    
    
    
    

    disp('----new test--------')
%     return
    
    normals=get_surface_normal_sph(ne,nex,ney,nexney,nx,ny,nz,n,nxny,n_eSurf,n_nSurf,nx_,X,Y,Z,parallel_computing);
    
%     normals(1,:)
%     
%     hold on
%     for i=1:1:n_eSurf(6)
%         plot3([normals(i,1),0],[normals(i,2),0],[normals(i,3),0])
%         
%     end
%     hold off
%     grid on
%     xlabel('x');ylabel('y');zlabel('z')
%     
%     return
    %--------------------------
    
    
    determinants=isoparametric_determinant_sph(...
        ne,neTotal,nex,ney,n,nx,ny,nz,nxny,nexney,X,Y,Z,...
        n_eSurf,n_nSurf,nx_,parallel_computing);
    nTotal=n+n_nSurf(6);
    neightTotal=neight+8*n_nSurf(6);
end


toc(counter_init)
fig_domain=[min(X) max(X) min(Y) max(Y) min(Z) max(Z)];
% THIS FUNCTION CAN BE PUT IN DIFFERENT MAIN FILE
toc(counter_init)
if parallel_computing==1
    if shape==1 || shape==2
        wTerms=generate_wTerms_curved_3d(ne,weights);
    elseif shape==3
        wTerms=generate_wTerms_curved_3d(neTotal,weights);
    end
end


if show_figure==1
    [x_minus_x,x_minus_y,x_minus_z,x_plus_x,x_plus_y,x_plus_z,...
        y_minus_x,y_minus_y,y_minus_z,y_plus_x,y_plus_y,y_plus_z,...
        z_minus_x,z_minus_y,z_minus_z,z_plus_x,z_plus_y,z_plus_z,...
        x_plot,y_plot,z_plot]=...
        setUp_illustration_curved_structure_3d(nx,ny,nz,X,Y,Z);
end

% --- Break down temperature gradient -------------------------------------

if shape==1 || shape==2
    [diffT,chi]=create_temperature_gradient_3d(...
        parallel_computing,shape,T,diff,entropy,T_theta,nx,ny,nexney,ne,nex,X,Y,Z);
elseif shape==3
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
if shape==1 || shape==2
    co=generate_co_3d(ci,include_fluc,ci_fluc,n,neight,nx,ny,nz);
elseif shape==3
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

% ==== Enter main loop ====================================================

fprintf('Enter main cycle\n')
while 1
    
% === Update concentration profile ========================================

    if shape==1 || shape==2
        if predict_with_acceleration==1
            cp=predict_c_with_acceleration_3d(c,co,coo,dt,dto,dtoo,...
                parallel_computing,neight);
            cooo=coo;
        elseif predict_with_acceleration==0
            cp=predict_c_3d(c,co,dt,dto,parallel_computing,neight);
        end
    elseif shape==3
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

% === Previous concentration (Only for parallel algorithm) ================

    %Test-------------
    
    conco__test=get_conc_type_I_curved_3d(ne,co,weights,nexney,nex,nx,ny);
    %---------------------------

    if parallel_computing==1
        if shape==1 || shape==2
            conco_=get_conc_type_I_curved_3d(ne,co,weights,nexney,nex,nx,ny); %-----
        elseif shape==3
            conco_=get_conc_type_I_curved_sph(neTotal,co,weights,ne,...
                nexney,nex,ney,nx,ny,nz,nxny,n,n_eSurf,n_nSurf,nx_);
        end
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

        if parallel_computing==1
            if shape==1 || shape==2
                fprintf('Computing conc...\n')
                sj=get_conc_curved_3d(ne,c,weights,nexney,nex,nx,ny); %-----
                if commit_sfsj_terms_to_RAM==1
                    fprintf('Computing sf and sj terms...\n')
                    [terms,sj]=compute_terms_curved_3d(ne,sj,chi,n1,n2,nex);
                    fprintf('Computing residual vector...\n')
                    sf=compute_sf_with_terms_curved_3d(neight,sj,conco_,nx,ny,nz,nex,ney,weights,dt,determinants,diffT,wTerms,terms); %-----
                    fprintf('Computing Jacobian matrix...\n')
                    sj=compute_sj_with_terms_curved_3d(gbfs1,gbfs2,nnz_,dt,weights,nex,ney,nx,ny,nz,determinants,diffT,wTerms,terms);
                elseif commit_sfsj_terms_to_RAM==0
                    fprintf('Computing residual vector...\n')
                    sf=compute_sf_curved_3d(neight,sj,conco_,nx,ny,nz,nex,ney,weights,dt,n1,n2,determinants,diffT,chi,wTerms);
                    fprintf('Computing Jacobian matrix...\n')
                    sj=compute_sj_curved_3d(gbfs1,gbfs2,nnz_,sj,dt,weights,nex,ney,nx,ny,nz,n1,n2,determinants,diffT,chi,wTerms);
                end
            elseif shape==3
                fprintf('Computing conc...\n')
                sj=get_conc_curved_sph(neTotal,c,weights,ne,nexney,nex,ney,nx,ny,nz,nxny,n,n_eSurf,n_nSurf,nx_);
                fprintf('Computing sf and sj terms...\n')
                [terms,sj]=compute_terms_curved_3d(neTotal,sj,chi,n1,n2,nex);
                
%                 terms_surface=compute_terms_curved_surface_sph(n_eSurf,conc_);

                terms_surface=compute_terms_curved_surface_sph(...
                    ne,nex,ney,nexney,nx,ny,nz,n,nxny,n_eSurf,n_nSurf,nx_,weights_surface,c);
                
                
                nnz(terms_surface)
                
%                 return
                
                fprintf('Computing residual vector...\n')
                
                
                sf=compute_sf_with_terms_curved_sph(...
                    neightTotal,sj,conco_,nx,ny,nz,nex,ney,weights,...
                    dt,determinants,diffT,wTerms,terms,...
                    n,ne,n_eSurf,n_nSurf,nTotal,nx_,...
                    determinants_surface,terms_surface,normals,weights_surface);
                
                
                
                
%                 return
                
                fprintf('Computing Jacobian matrix...\n')
                sj=compute_sj_with_terms_curved_sph(...
                    gbfs1,gbfs2,nnz_,dt,weights,nex,ney,nx,ny,nz,...
                    determinants,diffT,wTerms,terms,...
                    ne,nez,n,neight,nexney,n_eSurf,n_nSurf,nx_,...
                    wTerms_surface,weights_surface,determinants_surface,normals);
                
                
                
%                 disp('division')
%                 
%                 sj=sparse(gbfs1,gbfs2,sj);
%                 
%                 hey=sj\-sf'
%                 
%                 
%                 for i=1:1:8
%                     hey(157*8+i)
%                 end
%                 
%                 
%                 disp('================END HERE===================')
%                 return
                
                %---------------------------------------------------------
                
%                 sj_real=sparse(gbfs1,gbfs2,sj);
%                 sf_real=sf;
%                 
%                 disp('-------Start the test-----------')
%                 
%                 conc_=conco__test
%                 
%                 nnz_=nnz_sj_3d(nx,ny,nz);
%                 [gbfs1,gbfs2,iZeros,iOnes]=sj_mapper_3d(nx,ny,nz);
%                 
%                 sj=get_conc_curved_3d(ne,c,weights,nexney,nex,nx,ny);
%                 [terms,sj]=compute_terms_curved_3d(ne,sj,chi,n1,n2,nex);
%                 sf=compute_sf_with_terms_curved_3d(neight,sj,conco_,nx,ny,nz,nex,ney,weights,dt,determinants,diffT,wTerms,terms);
%                 fprintf('Computing Jacobian matrix...\n')
%                 sj=compute_sj_with_terms_curved_3d(gbfs1,gbfs2,nnz_,dt,weights,nex,ney,nx,ny,nz,determinants,diffT,wTerms,terms);
%                 
%                 sj_fake=sparse(gbfs1,gbfs2,sj);
%                 
%                 test1=zeros(1,nnz_);
%                 
%                 for i=1:1:nnz_
%                     difference=abs(sj_real(gbfs1(i),gbfs2(i))-sj_fake(gbfs1(i),gbfs2(i)));
%                     test1(i)=difference;
%                     if difference>1.0*10^-8
%                         
%                         [gbfs1(i),gbfs2(i)]
%                         error('bad')
%                     end
%                 end
%                 
%                 test2=zeros(1,neight);
%                 for i=1:1:neight
%                     test2(i)=abs(sf_real(i)-sf(i));
%                     
%                 end
%                 
%                 sum(test2)
%                 
%                 sum(test1)
%                 
%                 
%                 
%                 
% %                 disp('Nico nico ni')
% %                 
% %                 sum(sj)
% %                 
% %                 sj=sparse(gbfs1,gbfs2,sj);
% %                 
% %                 c=sj\sf'
% 
% 
% %                 
%                 disp('----------Current progress ends---------------')
%                 return
            end
            
% --- Apply BC --- %

%             [sf,sj]=bc_3d(sf,sj,iZeros,iOnes,gbfs1);

            %---Test---------------------------------
            
            
            
            
            %---------------------------------
            
% --- Sparsing --- %
            
            fprintf('Sparsing matrix...\n')
            sj=sparse(gbfs1,gbfs2,sj);
            
            
            %----Test---------------------\
            
%             sj1=zeros(neightTotal,neightTotal);
%             
%             for ijk=1:1:nnz_
%                 sj1(gbfs1(ijk),gbfs2(ijk))=sj(ijk);
%             end
%             
%             sj=sj1;
%             
%             how=0;
%             
%             for ijk=1:1:n_nSurf(6)
%                 how=how+1;
%                 [sf,sj]=raw_boundary_applier(sf,sj,n*8+(ijk-1)*8+2,neightTotal);
%                 [sf,sj]=raw_boundary_applier(sf,sj,n*8+(ijk-1)*8+5,neightTotal);
%                 [sf,sj]=raw_boundary_applier(sf,sj,n*8+(ijk-1)*8+6,neightTotal);
%                 [sf,sj]=raw_boundary_applier(sf,sj,n*8+(ijk-1)*8+8,neightTotal);
%                 
%                 [sf,sj]=raw_boundary_applier(sf,sj,n*8+(ijk-1)*8+3,neightTotal);
%                 [sf,sj]=raw_boundary_applier(sf,sj,n*8+(ijk-1)*8+4,neightTotal);
%                 [sf,sj]=raw_boundary_applier(sf,sj,n*8+(ijk-1)*8+7,neightTotal);
%             end
%             
%             how
% 
% %             nnz_
% %             length(gbfs1)
% %             
% %             max(gbfs1)
% %             max(gbfs2)
% %             size(sf)
% %             size(sj)
%             
%             
% %             return
%             
% %             sj_test=zeros(neightTotal,neightTotal);
% %             
% %             for i=1:1:nnz_
% %                 sj_test(gbfs1(i),gbfs2(i))=sj(gbfs1(i),gbfs2(i));
% %             end
% %             
% %             for irow=neight+1:1:neightTotal
% % %                 if mod(irow,8)~=1
% %                 if 1==1
% %                     for icol=1:1:neightTotal
% %                         if irow~=icol
% %                             sj_test(irow,icol)=1;
% %                         else
% %                             sj_test(irow,icol)=0;
% %                         end
% %                     end
% %                     sf(irow)=0;
% %                 end
% %             end
% %             sj=sj_test;
%             %-----------------------------
            
% === Compute using traditional method, wo parallelization ================            

        elseif parallel_computing==0
            
            fprintf('Computing residual vector and Jacobian matrix...\n')
            [sf,sj]=compute_sfsj_traditional_curved_3d(neight,ne,nexney,nex,nx,ny,...
                weights,c,co,determinants,dt,n1,n2,diffT,chi);
            
% --- Apply BC --- %

            %---------TESTING--------------------
            fprintf('Applying bc...\n')
            if shape==1
                [sf,sj]=bc_3d(sf,sj,iZeros,iOnes,gbfs1); %THIS IS WRONG!
            elseif shape==2
                [sf,sj]=bc_spherical_3d(sf,sj,iOnes,iZeros,surface_derivatives,gbfs1); %THIS DOES NOT WORK FOR TRADITIOANL METHOD!
            end
            %---------------------------------
            
% --- Sparsing --- %
            
            if sparse_traditional_sj==1
                fprintf('Sparsing matrix...\n')
                sj=sparse(sj);
            end
        end
        
        
% === Carry out matrix division ===========================================

        fprintf('Carrying out matrix division...\n')
        c_=sj\-sf';
           
% === Update the solution =================================================
        
        fprintf('Assessing current iterations...\n')
        c=c+c_';
        
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