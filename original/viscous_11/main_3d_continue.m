function main_3d_continue

% GUI is integrated

clear
clc

% === Start the timer =====================================================

counter_init=tic;

% === User defined variables ==============================================

% --- Simulation properties -----------------------------------------------

parallel_computing=1; % Non parallel computing uses traditional method
show_figure=1;
export_data=1; % Yes=1, No=0
tol=1.0e-6;
sparse_traditional_sj=0;
predict_with_acceleration=1;
commit_sfsj_terms_to_RAM=1;
fps=10;

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

face_colors={'blue','yellow','red'};
view_angle=3;
transparency=0.9;
edge_color='none';
ifig_conclusion=3;

% === User specifies simulation number ====================================

while 1
    run_number=inputdlg('Please enter simulation run number to be continued:','Select simulation');
    if isempty(run_number)
        termination_dialog()
        return
    end
    run_number=run_number{1};
    if ~isempty(run_number)
        run_number=str2double(run_number);
        if isempty(run_number)
            err_dialog_3d('Entry is not a number!')
            continue
        end
    else
        err_dialog_3d('Nothing was entered')
        continue
    end
    validation=validate_user_input_3d(run_number);
    if validation==0 || mod(run_number,1)~=0
        err_dialog_3d('Please enter only scalar integer number except infinity and NaN. Computation within the box is not allowed.')
        continue
    end
    break
end

% === Load data ===========================================================

[nex,ney,nez,obs_t,diff,include_fluc,ci_fluc,nr_tol,...
    nr_max_iteration,time_out,max_frame,entropy,T_theta,...
    dt_down,dt_up,dt_min,dt_ideal,bypass_tol,n1,n2,...
    elapsed,ci,T,dts,c,co,coo,time,...
    last_iteration,output_path,cancellation]=import_data_3d(run_number);

if cancellation==1
    termination_dialog()
    return
end

dt=dts(1);
dto=dts(2);
if predict_with_acceleration==1
    if length(dts)==3
        dtoo=dts(3);
    elseif length(dts)==2
        dtoo=inf;
    end
end

% === Validate inputs =====================================================

% VALIDATE THEM HERE

% === Set up variables ====================================================

time_goal=time_out-elapsed;

[nx,ny,nz,~,neight]=setUp_3d(nex,ney,nez);

[gbfs1,gbfs2,iZeros,iOnes]=sj_mapper_3d(nx,ny,nz);
nnz_=nnz_sj_3d(nx,ny,nz);

dx=1/nex;
dy=1/ney;
dz=1/nez;

weights=weight_adjuster_3d(generateWeights_3d(),dx,dy,dz);

ne=nex*ney*nez;
nexney=nex*ney;

dxyz=dx*dy*dz;

bypass=0;

c_obs=[mean(ci)-ci_fluc mean(ci) mean(ci)+ci_fluc];

x_coord=linspace(0,1,nx);
y_coord=linspace(0,1,ny);
z_coord=linspace(0,1,nz);

if parallel_computing==1 %This can be implemented on traditional algorithm as well
    wTerms=generate_wTerms_3d(weights);
end

if export_data==1
    if exist([output_path 'continuation/'],'dir')
        button = questdlg('There is a simulation in progress, which was probably disrupted. Continuing this simulation will delete all these data. Continue?','Warning','Yes','No','No');
        switch button
            case 'Yes'
                try
                    rmdir([output_path 'continuation/'],'s')
                catch
                    err_dialog_3d('Folder could not be removed! Please check the admin privilege.')
                    return
                end
            case 'No'
                termination_dialog()
                return
        end
    end
    mkdir([output_path 'continuation/'])
    original_output_path=output_path;
    output_path=[output_path 'continuation/'];
    mkdir([output_path 'concentration'])
    mkdir([output_path 'time'])
end

% --- Break down temperature gradient -------------------------------------

[diffT,chi]=T_gradient_3d(T,diff,nex,entropy,T_theta);

% === Phase diagram =======================================================

if export_data==1 || show_figure==1

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

% === Initialize the first time step ======================================

iframe=last_iteration;

% --- Plot the first time step --- %

if show_figure==1
    fig=figure(ifig_simulation);
    set(fig,'Position',[frame_pos(1) frame_pos(2) frame_size(1) frame_size(2)])
    clf(fig)
    visual(iframe)=illustrate_and_analyze_3d(c,c_obs,x_coord,y_coord,z_coord,...
        face_colors,edge_color,nx,ny,nz,view_angle,transparency,time);
end

% --- Output the first time step --- %

% This iteration has already been output

% ==== Enter main loop ====================================================

fprintf('Enter main cycle\n')
while 1
    
% === Update concentration profile ========================================

    if predict_with_acceleration==1
        cp=predict_c_with_acceleration_3d(c,co,coo,dt,dto,dtoo,...
            parallel_computing,neight);
        cooo=coo;
    elseif predict_with_acceleration==0
        cp=predict_c_3d(c,co,dt,dto,parallel_computing,neight);
    end

    coo=co;
    co=c;
    c=cp;
    
% === Previous concentration (Only for parallel algorithm) ================

    if parallel_computing==1
        conco_=get_conc_type_I_3d(ne,co,weights,nexney,nex,nx,ny);
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
        fprintf('Remaining time: %d h\n',(time_goal-toc(counter_init))/3600)
        
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
            fprintf('Computing conc...\n')
            sj=get_conc_3d(ne,c,weights,nexney,nex,nx,ny);
            if commit_sfsj_terms_to_RAM==1
                fprintf('Computing sf and sj terms...\n')
                [terms,sj]=compute_terms_3d(ne,sj,chi,n1,n2,nex);
                fprintf('Computing residual vector...\n')
                sf=compute_sf_with_terms_3d(neight,sj,conco_,nx,ny,nz,nex,ney,weights,dt,dxyz,diffT,wTerms,terms);
                fprintf('Computing Jacobian matrix...\n')
                sj=compute_sj_with_terms_3d(gbfs1,gbfs2,nnz_,dt,weights,nex,ney,nx,ny,nz,dxyz,diffT,wTerms,terms);
            elseif commit_sfsj_terms_to_RAM==0
                fprintf('Computing residual vector...\n')
                sf=compute_sf_3d(neight,sj,conco_,nx,ny,nz,nex,ney,weights,dt,n1,n2,dxyz,diffT,chi,wTerms);
                fprintf('Computing Jacobian matrix...\n')
                sj=compute_sj_3d(gbfs1,gbfs2,nnz_,sj,dt,weights,nex,ney,nx,ny,nz,n1,n2,dxyz,diffT,chi,wTerms);
            end
            
% --- Apply BC --- %

            fprintf('Applying bc...\n')
            [sf,sj]=bc_3d(sf,sj,iZeros,iOnes,gbfs1);
            
% --- Sparsing --- %
            
            fprintf('Sparsing matrix...\n')
            sj=sparse(gbfs1,gbfs2,sj);
            
% === Compute using traditional method, wo parallelization ================            

        elseif parallel_computing==0
            
            fprintf('Computing residual vector and Jacobian matrix...\n')
            [sf,sj]=compute_sfsj_traditional_3d(neight,ne,nexney,nex,nx,ny,...
                weights,c,co,dxyz,dt,n1,n2,diffT,chi);
            
            %------------------------------------- TESTING
%             disp('Start duh test')
%             toc(counter_init)
%             wTerms=generate_wTerms_3d(weights);
%             conco_=get_conc_type_I_3d(ne,co,weights,nexney,nex,nx,ny);
%             conc_=get_conc_3d(ne,c,weights,nexney,nex,nx,ny);
%             
%             %--------------No term testing--------------------------
% %             sf2=compute_sf_3d(neight,conc_,conco_,nx,ny,nz,nex,ney,weights,dt,n1,n2,dxyz,diffT,chi,wTerms);
% %             toc(counter_init)
% %             sj2=compute_sj_3d(gbfs1,gbfs2,nnz_,conc_,dt,weights,nex,ney,nx,ny,nz,n1,n2,dxyz,diffT,chi,wTerms);
% %             sj2=sparse(gbfs1,gbfs2,sj2);
%             %---------------------------------------------------------------
% 
%             %-------------- Term testing -----------------------------
%             [terms,conc_]=compute_terms_3d(ne,conc_,chi,n1,n2,nex);
%             toc(counter_init)
%             sf2=compute_sf_with_terms_3d(neight,conc_,conco_,nx,ny,nz,nex,ney,weights,dt,dxyz,diffT,wTerms,terms);
%             %-
%             toc(counter_init)
%             fprintf('Computing sj\n')
%             sj2=compute_sj_with_terms_3d(gbfs1,gbfs2,nnz_,dt,weights,nex,ney,nx,ny,nz,dxyz,diffT,wTerms,terms);
%             fprintf('Done\n')
%             toc(counter_init)
%             %-
%             sj2=sparse(gbfs1,gbfs2,sj2);
%             %------------------------------------------------------------
%             
%             %--------------sj properties testing-------------------------
% %             terms=compute_terms_3d(ne,conc_,chi,n1,n2,nex);
% %             toc(counter_init)
% %             sj_properties=characterize_sjth_3d(nnz_,gbfs1,gbfs2,nex,ney,nx,ny,nz,diffT);
% %             toc(counter_init)
% %             sf2=compute_sf_with_terms_3d(neight,conc_,conco_,nx,ny,nz,nex,ney,weights,dt,dxyz,diffT,wTerms,terms);
% %             %-
% %             toc(counter_init)
% %             fprintf('Computing sj\n')
% % %             sj2=test26(nnz_,dt,weights,dxyz,diffT,wTerms,terms,sj_properties);
% %             sj2=test27(gbfs1,gbfs2,nnz_,dt,weights,nex,ney,nx,ny,nz,dxyz,diffT,wTerms,terms,sj_properties);
% %             
% %             fprintf('Done\n')
% %             toc(counter_init)
% %             %-
% %             nnz(sj_properties)
% %             nnz_*48
% %             size(gbfs1)
% %             size(terms)
% %             nnz(terms)
% %             sj2=sparse(gbfs1,gbfs2,sj2);
%             %------------------------------------------------------------
%             
%             toc(counter_init)
%             sf_diff=sum(abs(sf-sf2))
%             
%             sj_diff=0;
%             for ijk=1:1:nnz_
%                 sj_diff=sj_diff+abs(sj(gbfs1(ijk),gbfs2(ijk))-sj2(gbfs1(ijk),gbfs2(ijk)));
%                 if sj2(gbfs1(ijk),gbfs2(ijk))==0 || sj(gbfs1(ijk),gbfs2(ijk))==0
%                     ijk
%                     sj2(gbfs1(ijk),gbfs2(ijk))
%                     sj(gbfs1(ijk),gbfs2(ijk))
%                     error('Zero is found!')
%                 end
%             end
%             sj_diff
%             
%             disp('Finland')
%             return
            %-------------------------------------
            
% --- Apply BC --- %

            fprintf('Applying bc...\n')
            [sf,sj]=bc_traditional_3d(sf,sj,gbfs1,gbfs2,iZeros,iOnes);
            
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
        visual(iframe)=illustrate_and_analyze_3d(c,c_obs,x_coord,y_coord,z_coord,...
            face_colors,edge_color,nx,ny,nz,view_angle,transparency,time);
    end

% === Determine whether to end simulation =================================

    if iframe==max_frame
        disp('Maximum frame number specified has been reached!')
        if export_data==1
            termination_type=1;
        end
        break
    elseif time_goal<toc(counter_init)
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
    
    timeFolder=[original_output_path 'time/'];
    concFolder=[original_output_path 'concentration/'];
    if parallel_computing==1
        parfor iteration=last_iteration+1:1:iframe
            timeFile=[output_path 'time/iteration_' num2str(iteration)];
            concFile=[output_path 'concentration/iteration_' num2str(iteration)];
            movefile(timeFile,timeFolder)
            movefile(concFile,concFolder)
        end
    elseif parallel_computing==0
        for iteration=last_iteration+1:1:iframe
            timeFile=[output_path 'time/iteration_' num2str(iteration)];
            concFile=[output_path 'concentration/iteration_' num2str(iteration)];
            movefile(timeFile,timeFolder)
            movefile(concFile,concFolder)
        end
    end
    
    sf=[nex;ney;nez;dt;obs_t;diff;include_fluc;ci_fluc;nr_tol;...
        nr_max_iteration;time_out;max_frame;entropy;T_theta;...
        dt_down;dt_up;dt_min;dt_ideal;bypass_tol;n1;n2];
    dlmwrite([original_output_path 'specification/spec_general'],sf);
    dlmwrite([original_output_path 'specification/time_steps'],steps);
    dlmwrite([original_output_path 'specification/elapsed'],toc(counter_init)+elapsed);
    dlmwrite([original_output_path 'specification/termination'],termination_type);
    rmdir(output_path,'s')
end

% === Export the movie ====================================================

if show_figure==1
    try
        
        for iframe_previous=1:1:last_iteration-1
            clf(fig)
            c_previous=dlmread([original_output_path 'concentration/' 'iteration_' num2str(iframe_previous)]);
            time_previous=dlmread([original_output_path 'time/' 'iteration_' num2str(iframe_previous)]);
            visual(iframe_previous)=illustrate_and_analyze_3d(c_previous,c_obs,x_coord,y_coord,z_coord,...
                face_colors,edge_color,nx,ny,nz,view_angle,transparency,time_previous);
        end
        
        if length(ci)==1 && length(T)==1
            vid_title=sprintf('3D nonlinear Cahn Hilliard,%dx%dx%d,ci=%1.3f,fluc=%1.4f,diff=%d,T=%0.3f,elapsed time=%1.3fs.avi',nex,ney,nez,ci,ci_fluc*include_fluc,diff,T,toc(counter_init)+elapsed);
        elseif length(ci)==2 && length(T)==1
            vid_title=sprintf('3D nonlinear Cahn Hilliard,%dx%dx%d,ci=[%1.2f,%1.2f],fluc=%1.4f,diff=%d,T=%0.3f,elapsed time=%1.3fs.avi',nex,ney,nez,ci(1),ci(2),ci_fluc*include_fluc,diff,T,toc(counter_init)+elapsed);
        elseif length(ci)==1 && length(T)==2
            vid_title=sprintf('3D nonlinear Cahn Hilliard,%dx%dx%d,ci=%1.3f,fluc=%1.4f,diff=%d,T=[%1.2f,%1.2f],elapsed time=%1.3fs.avi',nex,ney,nez,ci,ci_fluc*include_fluc,diff,T(1),T(2),toc(counter_init)+elapsed);
        elseif length(ci)==2 && length(T)==2
            vid_title=sprintf('3D nonlinear Cahn Hilliard,%dx%dx%d,ci=[%1.2f,%1.2f],fluc=%1.4f,diff=%d,T=[%1.2f,%1.2f],elapsed time=%1.3fs.avi',nex,ney,nez,ci(1),ci(2),ci_fluc*include_fluc,diff,T(1),T(2),toc(counter_init)+elapsed);
        end

        vid_title=[original_output_path vid_title];
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