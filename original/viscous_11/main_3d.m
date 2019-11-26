function main_3d

clear
clc

% === Start the timer =====================================================

main_start=tic;

% === User defined variables ==============================================

% --- Define mesh ---------------------------------------------------------

nex=6;
ney=6;
nez=6;

% --- Polymer properties --------------------------------------------------

n1=1;
n2=10;

%======================
diff=1800
ci_list=[...
    0.4 0.4
    0.6 0.6
    0.45 0.45
    0.45 0.45
    0.45 0.45];
T_list=[...
    0.55 0.55
    0.3 0.3
    0.4 0.4
    0.5 0.55
    0.6 0.6];

nSimulations=1;
%======================

include_fluc=1;
ci_fluc = 0.001*1;

% ci_fluc=0.3

entropy=1;
T_theta=1;

ci_critical=0; % This bypasses ci specified above
T_critical=0; % This bypasses T specified above

% --- Simulation properties -----------------------------------------------

parallel_computing=1; % Non parallel computing uses traditional method

show_figure=1;
export_data=0; % Yes=1, No=0


dt=1.0*10^-6;

obs_t=2.0e-1;

nr_tol=1.0e-6;
nr_max_iteration=20;
%======================
time_out=3600*(7/60);
max_frame=1000;
%======================
fps=10;
dt_down=0.5;
dt_up=1.2;
dt_min=1.0e-15;
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

proceed_request=0;

ifig_phase=1;

% --- Analysis specifications ---------------------------------------------

alpha=200;
beta=8;

% --- Illustration specifications -----------------------------------------

frame_pos=[0 0];
% frame_size=[720 720];

frame_size=[1280 720];

ifig_simulation=2;

face_colors={'blue','yellow','red'};
view_angle=3;
transparency=0.9;
edge_color='none';
ifig_conclusion=3;

% === Validate inputs =====================================================

% VALIDATE THEM HERE

% === Set up variables ====================================================

[nx,ny,nz,n,neight]=setUp_3d(nex,ney,nez);

[gbfs1,gbfs2,iZeros,iOnes]=sj_mapper_3d(nx,ny,nz);
nnz_=nnz_sj_3d(nx,ny,nz);

dx=1/nex;
dy=1/ney;
dz=1/nez;

weights=weight_adjuster_3d(generateWeights_3d(),dx,dy,dz);

ne=nex*ney*nez;
nexney=nex*ney;

dxyz=dx*dy*dz;

dto=inf;
if predict_with_acceleration==1
    dtoo=inf;
end

bypass=0;

x_coord=linspace(0,1,nx);
y_coord=linspace(0,1,ny);
z_coord=linspace(0,1,nz);

if parallel_computing==1 %This can be implemented on traditional algorithm as well
    wTerms=generate_wTerms_3d(weights);
end

for iSimulation=1:1:nSimulations
    
    counter_init=tic;
    
    ci=ci_list(iSimulation,:);
    T=T_list(iSimulation,:);
    if T(2)==T(1)
        T=T(1);
    end
    if ci(2)==ci(1)
        ci=ci(1);
    end
    
    if iSimulation==1
        if ci_critical==1 && T_critical==1
            [ci,T]=identify_critical_point(n1,n2,entropy,T_theta);
            [nrows,ncols]=size(ci_list);
            ci_list=ones(nrows,ncols)*ci;
            [nrows,ncols]=size(T_list);
            T_list=ones(nrows,ncols)*T;
        elseif ci_critical==1
            [ci,~]=identify_critical_point(n1,n2,entropy,T_theta);
            [nrows,ncols]=size(ci_list);
            ci_list=ones(nrows,ncols)*ci;
        elseif T_critical==1
            [~,T]=identify_critical_point(n1,n2,entropy,T_theta);
            [nrows,ncols]=size(T_list);
            T_list=ones(nrows,ncols)*T;
        end
    end
    
    c_obs=[mean(ci)-ci_fluc mean(ci) mean(ci)+ci_fluc];

    % --- Break down temperature gradient -------------------------------------

    [diffT,chi]=T_gradient_3d(T,diff,nex,entropy,T_theta);

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

        if iSimulation==1
            fprintf('Generating phase diagram...\n')
            [xspinodal,yspinodal,xbinodal,ybinodal]=spinodal_binodal(...
                n1,n2,entropy,T_theta,nT,tol,parallel_computing,nInitialGuessPts,...
                guess_accuracy);
            if show_figure==1
                figure(ifig_phase)
                illustrate_simulations(xspinodal,yspinodal,xbinodal,ybinodal,T_list,ci_list,nSimulations)
                if proceed_request==1
                    proceed=proceed_dialog();
                    if proceed~=1
                        return
                    end
                end
            end
        end
    end

% === Generate initial condition ==========================================

%     co=generate_co_3d(ci,include_fluc,ci_fluc,n,neight,nx,ny,nz);
    
    co=ic_sphere(ci,include_fluc,ci_fluc*10,n,0.5,neight,nx,ny,nz,nex,ney,nez);

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
        visual=illustrate_and_analyze_3d(c,c_obs,x_coord,y_coord,z_coord,...
            face_colors,edge_color,nx,ny,nz,view_angle,transparency,time);
    end
    
%     return

    % --- Output the first time step --- %

    if export_data==1
        dlmwrite([output_path 'time/iteration_1'],0);
        dlmwrite([output_path 'concentration/iteration_1'],c);
    end

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
        
        %----TEST--------------------------------------
%         iStopper=0;
        %------------------------------------------------

        while err>nr_tol
            
            %----TEST--------------------------------------
%             iStopper=iStopper+1;
            %------------------------------------------------

            fprintf('Current simulation is run number %i\n',iSimulation)
            fprintf('Elapsed time for entire program is %f\n',toc(main_start))
            fprintf('Elapsed time for this particular run is %f\n',toc(counter_init))

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
                    [terms,sf,sj]=compute_terms_3d(ne,sj,chi,n1,n2,nex);
                    fprintf('Computing residual vector...\n')
                    sf=compute_sf_with_terms_3d(neight,sf,conco_,nx,ny,nz,nex,ney,weights,dt,dxyz,diffT,wTerms,terms);
                    fprintf('Computing Jacobian matrix...\n')
                    sj=compute_sj_with_terms_3d(gbfs1,gbfs2,nnz_,dt,weights,nex,ney,nx,ny,nz,dxyz,diffT,wTerms,terms,sj);
                    
                    %---TESTING---------------------------------------------------------------------
                    
                    %----------------------------------------------------------------------------
                    
                elseif commit_sfsj_terms_to_RAM==0
                    fprintf('Computing residual vector...\n')
                    sf=compute_sf_3d(neight,sj,conco_,nx,ny,nz,nex,ney,weights,dt,n1,n2,dxyz,diffT,chi,wTerms);
                    fprintf('Computing Jacobian matrix...\n')
                    sj=compute_sj_3d(gbfs1,gbfs2,nnz_,sj,dt,weights,nex,ney,nx,ny,nz,n1,n2,dxyz,diffT,chi,wTerms);
                end

    % --- Apply BC --- %

                fprintf('Applying bc...\n')
%                 warning('NOT APPLYING  BC NOW')
                [sf,sj]=bc_3d(sf,sj,iZeros,iOnes,gbfs1);

    % --- Sparsing --- %

                fprintf('Sparsing matrix...\n')
                sj=sparse(gbfs1,gbfs2,sj);


%--- TEST-----------------
%                 sj=generate_traditional_mx(gbfs1,gbfs2,sj);
%                 
%                 size(sj)
%                 rank(sj)
%                 
%                 disp('yume')
%                 return
                %----------------------------------

% === Compute using traditional method, wo parallelization ================            

            elseif parallel_computing==0

                fprintf('Computing residual vector and Jacobian matrix...\n')
                [sf,sj]=compute_sfsj_traditional_3d(neight,ne,nexney,nex,nx,ny,...
                    weights,c,co,dxyz,dt,n1,n2,diffT,chi);

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
        
        %---Test-------------------------------
        
%         indexer=-7;
%         for ijk=1:1:n
%             indexer=indexer+8;
%             if rand(1,1)>=0.5
%                 c(indexer)=c(indexer)+ci_fluc*include_fluc*rand(1,1);
%             else
%                 c(indexer)=c(indexer)-ci_fluc*include_fluc*rand(1,1);
%             end
%         end
        
        %--------------------------------------

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
        
        %-------------
%         return
        %----------------

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
        dlmwrite([output_path 'specification/shape'],1);
    end

% === Export the movie ====================================================

    if show_figure==1
        try
            if length(ci)==1 && length(T)==1
                vid_title=[output_path sprintf('3D nonlinear Cahn Hilliard,%dx%dx%d,ci=%1.3f,fluc=%1.4f,diff=%d,T=%0.3f,elapsed time=%1.3fs.avi',nex,ney,nez,ci,ci_fluc*include_fluc,diff,T,toc(counter_init))];
            elseif length(ci)==2 && length(T)==1
                vid_title=[output_path sprintf('3D nonlinear Cahn Hilliard,%dx%dx%d,ci=[%1.2f,%1.2f],fluc=%1.4f,diff=%d,T=%0.3f,elapsed time=%1.3fs.avi',nex,ney,nez,ci(1),ci(2),ci_fluc*include_fluc,diff,T,toc(counter_init))];
            elseif length(ci)==1 && length(T)==2
                vid_title=[output_path sprintf('3D nonlinear Cahn Hilliard,%dx%dx%d,ci=%1.3f,fluc=%1.4f,diff=%d,T=[%1.2f,%1.2f],elapsed time=%1.3fs.avi',nex,ney,nez,ci,ci_fluc*include_fluc,diff,T(1),T(2),toc(counter_init))];
            elseif length(ci)==2 && length(T)==2
                vid_title=[output_path sprintf('3D nonlinear Cahn Hilliard,%dx%dx%d,ci=[%1.2f,%1.2f],fluc=%1.4f,diff=%d,T=[%1.2f,%1.2f],elapsed time=%1.3fs.avi',nex,ney,nez,ci(1),ci(2),ci_fluc*include_fluc,diff,T(1),T(2),toc(counter_init))];
            end
            create_movie_3d(visual,vid_title,fps)

% === Conclusion ==========================================================

            fig=figure(ifig_conclusion);
            set(fig,'Position',[frame_pos(1) frame_pos(2) frame_size(1) frame_size(2)])
            movie(gcf,visual,1)
        catch
            warning('Failed to create the video')
        end
    end
    
    fprintf('End of one simulation\n')
    toc(counter_init)
end

fprintf('End of entire main program\n')
toc(main_start)

end