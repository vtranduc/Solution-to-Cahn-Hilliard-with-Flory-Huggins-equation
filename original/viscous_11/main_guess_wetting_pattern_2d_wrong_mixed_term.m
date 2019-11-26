function main_guess_wetting_pattern_2d_wrong_mixed_term

% BC for wetting is WRONG!
% HERE, BC FOR THE MIXED THEM IS LOOSELY DEFINED AS ZERO
% ONLY THE FIRST DERIVATIVE HAS BEEN CORRECTED
% NONETHELESS, IT STILL GIVES A GOOD GUESS FOR 3D IMAGE!

if nargin==0
    clear
end
clc
close all

% --- Start the timer -----------------------------------------------------

counter_init=tic;

% --- User defined variables ----------------------------------------------

% --- Specify the mesh ---

nex=31;
ney=31;

% --- Polymer properties ---

if ~exist('diff','var')
    diff=5000;
end
if ~exist('T','var')
%     T=[0.3 0.31];
    T=0.35;
end
if ~exist('ci','var')
    ci=0.3;
end
if ~exist('n2','var')
    n2=1; %Must be 1 or larger!
end

%---------------------------

wetting=1;

hl=0.0*ones(1,ney+1);
hb=0*ones(1,nex+1);
hr=0*ones(1,ney+1);
ht=0*ones(1,nex+1);

gl=0.0*ones(1,ney+1);
gb=0*ones(1,nex+1);
gr=-100*ones(1,ney+1);
gt=0.0*ones(1,nex+1);

%---------------------

%============================
thermophoresis=0;
distribution_type=4;
% distribution_type==0 Linear temperature gradeint
% distribution_type==3 Cosine wave throughout entire square domain
% distribution_type==4 Temperature focussed at the center
%============================

n1=1; %Must be equal to 1!
include_fluc=1;
ci_fluc = 10^-6; %Initial fluctuation
entropy=1;
T_theta=1;

% --- Simulation properties ---

%======================
time_out=3600*(60/60)*1;
dt=1.0e-10;
% If too big, simulation will stop in first time step
% General rule: The more it is inside unstable region, the smaller it
% should be
proceed_request=0; % Yes=1, No=0
%======================
max_frame=100000;
nr_tol=1.0e-6;
nr_max_iteration=5;
nr_up=2;
dt_down=0.5;
dt_up=1.1;
dt_min=1.0e-20;
dt_min_=dt_min;
dt_ideal=5.0e-10;
bypass_tol=5;
xlen=1;ylen=1;
obs_t=10000000000000000000000000000000000000000000000000; %Observation time
%Declare amount of workers to be used
nworkers=0;
% nworkers=0 if you wanna use maximum number of cores
% nworkers=-1 if you wanna specify it as user inpput
% nworkers=1 if no parallel computing is desired
% nworkers>1 for parallel computing with specified number of workers
communication_reliever=1;
%Phase diagram specifications
nT=100;
tol=1.0e-6;
nInitialGuessPts=100;
guess_accuracy=10^-3;

% --- Analysis and exporting methods ---

alpha=200;


%=================================
beta=8;
% beta=100;
%=================================


limS=10^-12;
limSstretcher=2;
fc=20; %Put fc=0 if you wanna use Nyquist frequency

%==========================================================================
%==========================================================================
%==========================================================================
%==========================================================================
%==========================================================================
%==========================================================================
%==========================================================================
figure_analysis=3; 
%figure_analysis=0 displays no figure
%figure_analysis=1 for full analysis
%figure_analysis=2 for just surface plot of concentration
%figure_analysis=3 for basic fourier analysis
export_figure=0;
export_data=0; % Yes=1, No=0
simulation_continuation=0;
% simulation_continuation=0 if this is a new simulation
% simulation_continuation=1 to continue previously completed simulation.
% This will automatically export data
view_simulation_only=0;
% simulation_continuation=1 This will only create a new video for finished
% simulation. Any video in the folder will be deleted and replaced.
%==========================================================================
%==========================================================================
%==========================================================================
%==========================================================================
%==========================================================================
%==========================================================================
%==========================================================================


if export_figure==1
    fig_output_dim=[0 0 6 4]; %In PaperUnits Inches
    export_fig=1000;
    dpi=500;
    axis_font='CMU Serif';
    font_size=9;
    export_spec=1;
    k_specific=0; %Zero will automatically compute characteristic frequency
elseif export_figure==0
    export_spec=NaN;
    k_specific=NaN;
end

video_method=2;
%video_method=1 for specifying fps manually
%video_method=2 for specifying the duration of the video
if video_method==1
    fps=10;duration=NaN;
elseif video_method==2
    fps=NaN;duration=10;
end

if export_data==1
    sig_dig=17;
end

if simulation_continuation==1
    simulation_folder=[pwd '/Results/Run_' num2str(46) '/'];
end

if view_simulation_only==1
    simulation_folder=[pwd '/Results/Run_' num2str(1) '/'];
    
%     simulation_folder=[pwd '/Results/10000 0.3 0.32 1/'];
    if figure_analysis==0
        warning('No figure has been specified. We will use figure_analysis==3')
        figure_analysis=3;
    end
end

% Specifications for each type of analysis
if figure_analysis==1
    contour_lgd_font_size=7.5;
    nx_grad=20;ny_grad=20;
    frame_size=[1920 1080]*1.0;
    frame_pos=[0 0];
elseif figure_analysis==2
    contour_lgd_font_size=NaN;
    nx_grad=NaN;ny_grad=NaN;
    frame_size=[720 720];
    frame_pos=[100 100];
elseif figure_analysis==3
    contour_lgd_font_size=NaN;
    nx_grad=NaN;ny_grad=NaN;
    frame_size=[1280 720];
    frame_pos=[100 100];
end

% --- Validate inputs -----------------------------------------------------

%Do validadation of all variables here

% --- Reload previous simulation if specified -----------------------------

if view_simulation_only==1
    if simulation_continuation==0
        simulation_continuation=1;
    end
    if export_data==1
        export_data=0;
    end
end

if simulation_continuation==1
    fprintf('Loading data...\n')
    [goal_reached,modification,output_path,nex,ney,diff,...
    include_fluc,ci_fluc,nr_tol,nr_max_iteration,entropy,T_theta,...
    dt_down,dt_up,dt_min,dt_ideal,bypass_tol,n1,n2,dt,...
    sf,c,co,iframe,T,ci,time]=reload_2d(...
    simulation_folder,obs_t,dt,max_frame,...
    view_simulation_only);
    if view_simulation_only==0
        if goal_reached==1
            disp('Since the goal has been reached, simulation will terminate now')
            return
        end
        if modification==1
            dlmwrite([output_path 'specification/spec_general'],sf);
        end
        if export_data==0
            export_data=1;
        end
    end
    nframes_last_sim=iframe;
elseif simulation_continuation==0
    nframes_last_sim=0;simulation_folder=NaN;
end

% --- Initialize multi workers as specified -------------------------------

nworkers=cpu_initializer(nworkers);

% --- Phase diagram -------------------------------------------------------

if proceed_request==1 || figure_analysis==1 || figure_analysis==3 || export_figure==1
    fprintf('Generating free energy values...\n')
    %===============T_theta MUST NO LONGER PLAY ANY ROLE IN OUR NONDIMENSIONALIZED SIMULATION!!!
    [xspinodal,yspinodal,xbinodal,ybinodal]=spinodal_binodal(...
        n1,n2,entropy,T_theta,nT,tol,1,nInitialGuessPts,...
        guess_accuracy);
else
    xspinodal=NaN;yspinodal=NaN;xbinodal=NaN;ybinodal=NaN;
end

if proceed_request==1
    figure(100)
    illustrate_simulation(xspinodal,yspinodal,xbinodal,ybinodal,T,ci)
    proceed=input('Proceed? Yes=1+ENTER or Yes=ENTER, No=Any other numerical+ENTER or string+ENTER\n');
    if proceed~=1
        disp('The simulation has been cancelled by the user')
        return
    else
        disp('The simulation is to proceed')
    end
end

% --- Setup other global variables ----------------------------------------

%T_theta does not even make sense anymore!!!

%------------------------------------------



T_distribution_spec=[0.0 0.2];

%-----------FIXXXXXXXXXXXXXXXXXXXXXXXXXX quick

[ne,nx,ny,n,nfour,ll,li,lu,bl,bi,bu,tl,ti,tu,rl,ri,ru,...
    fc,frequency_domain,x_coord,y_coord,T,...
    nnz_,irow,icol,weights,ci_ave,dxdy,...
    Two_chi_n1,coef_T,cp,bypass,dto,...
    nCoreTasks,index_array,...
    grad_T,~,...
    wTerms,iZeros,iOnes,k_specific]=...
    setUpGlobal_2d(nex,ney,fc,alpha,xlen,ylen,...
    ci,T,entropy,diff,n1,n2,...
    communication_reliever,nworkers,...
    dt,distribution_type,T_distribution_spec,...
    export_figure,k_specific);


%======================================Construction
%======================================Construction
%======================================Construction
%======================================Construction
%======================================Construction
%======================================Construction
%======================================Construction
%======================================Construction




% coef_T=diff*T;
if length(T)~=1
    error('just testing here')
end
size(weights);
if wetting==1
    
    sj_assist=wetting_assist(weights);
    
else
    sj_assist=NaN;
    
end

% warning('adjust here')



% splitting=8;
% 
% for i=1:1:round(ny/splitting)
%     hr(i)=0;
%     gr(i)=0;
% end
% 
% for i=round(ny/splitting*2):1:round(ny/splitting*3)
%     hr(i)=0;
%     gr(i)=0;
% end
% 
% for i=round(ny/splitting*4):1:round(ny/splitting*5)
%     hr(i)=0;
%     gr(i)=0;
% end
% 
% for i=round(ny/splitting*6):1:round(ny/splitting*7)
%     hr(i)=0;
%     gr(i)=0;
% end

% hr'
% 
% return


% hl=0.0*ones(1,ney+1);
% hb=0*ones(1,nex+1);
% hr=10*ones(1,ney+1);
% ht=0*ones(1,nex+1);
% 
% gl=0.0*ones(1,ney+1);
% gb=0*ones(1,nex+1);
% gr=-15*ones(1,ney+1);
% gt=0.0*ones(1,nex+1);

if wetting==1
    [iZeros_wetting,iOnes_wetting,ig]=get_bc_keyVals_wetting(...
        nx,ny,nnz_,irow,icol,hl,hb,hr,ht,gl,gb,gr,gt);
end

%======================================Construction
%======================================Construction
%======================================Construction
%======================================Construction
%======================================Construction
%======================================Construction
%======================================Construction

% test=zeros(ney,nex);
% e=0;
% for exth=1:1:nex
%     for eyth=1:1:ney
%         e=e+1;
%         test(eyth,exth)=coef_T(e,2,2,3);
%     end
% end
% 
% coef_T(e,2,2,3);
% 
% surf(test)
% xlabel('x');ylabel('y')
% % 
% return

[terms_assist,to_be_deleted]=assist_global(...
    grad_T,ne,coef_T,entropy,n1,weights,wTerms);

if figure_analysis~=0
    [nx_grad,ny_grad,x_grad,y_grad]=setUpFigure_2d(...
        figure_analysis,nx,ny,nx_grad,ny_grad);
end



if simulation_continuation==0
    if export_data==1 || figure_analysis~=0 || export_figure==1
        output_path=folder_setUp(export_data,export_figure,export_spec,k_specific);
    end
    if export_data==1
        sf=[nex;ney;dt;obs_t;diff;include_fluc;ci_fluc;nr_tol;...
            nr_max_iteration;time_out;max_frame;entropy;T_theta;...
            dt_down;dt_up;dt_min;dt_ideal;bypass_tol;n1;n2];
        dlmwrite([output_path 'specification/spec_general'],sf);
        dlmwrite([output_path 'specification/status'],0);
        dlmwrite([output_path 'specification/spec_T'],T);
        dlmwrite([output_path 'specification/spec_ci'],ci);
    end
end

%--------------Create initial concentration profile----------------------

if simulation_continuation==0
    % Must be careful to apply bc in co.

    co=generate_co_2D(ci,nx,ny,include_fluc,ci_fluc);

%     ci_fluc=-0.05;radius=0.2;
% %     ci_fluc=0.25;radius=0.2;
% %     ci_fluc=0.001;radius=0.4995;
%     co=ic_circle(ci,include_fluc,ci_fluc,n,nx,ny,nfour,radius,nex,ney);
%     
%     
%     test=extract_nodal_weights_2D(co,nx,ny);
%     test=reshape(test,1,n);
%     
%     [min(test) mean(test) max(test) ci_ave]
%     
%     spinodal=get_specific_spinodal_binodal(T,n1,n2,entropy,10^-6)
%     return

    % co=dlmread('hoshizora');

    % co=create_sine_wave_2d(nx,ny,nfour,A,w,x_coord,y_coord,ci);

%     co=dlmread('iteration_198');
    
%     error('dfadfa')

% plot([0 1],[0 1])
% xlabel('\fontname{CommercialScript BT}F')
% return
    
    if wetting==1
        co=bc_wetting(co,ll,li,lu,bl,bi,bu,rl,ri,ru,tl,ti,tu,...
            gl,gb,gr,gt,hl,hb,hr,ht);
    end

%-------------------------------------------------------------------

% --- Initialize the first time step --------------------------------------

    time=0.0;
    c=co;
    iframe=1;
    if export_data==1
        dlmwrite([output_path 'time/iteration_1'],0);
        dlmwrite([output_path 'concentration/iteration_1'],c,'precision',sig_dig);
    end
    
    
    if export_figure==1
        [logS_export,limS_export,E_updated_export,...
            max_c_export,min_c_export,intensity_specific]=plots_output(...
            fig_output_dim,dpi,axis_font,font_size,...
            output_path,export_fig,...
            c,nx,ny,alpha,beta,fc,ci_ave,...
            iframe,time,NaN,limS,xlen,ylen,limSstretcher,...
            xspinodal,yspinodal,xbinodal,ybinodal,T,ci,...
            frequency_domain,...
            coef_T,n1,n2,ne,nex,ney,Two_chi_n1,dxdy,weights,grad_T,... %CHECK LATER====================================
            NaN,NaN,NaN,diff,...
            export_spec,k_specific,NaN,...
            nworkers);
        
    end
end

% --- Analyze and plot the initial condition ------------------------------

if figure_analysis~=0
    
    if simulation_continuation==1
        if figure_analysis==1 || figure_analysis==3
            fprintf('Loading figures...\n')
            [logS,limS_,ac_t,E]=reload_fig_2d(....
                simulation_folder,figure_analysis,iframe,nx,ny,...
                alpha,beta,fc,ci_ave,limS,limSstretcher,nworkers); %=============CHECK LATER===========================================================
        end
        if figure_analysis==2
            logS=0;E=0;ac_t=[0 0 0];limS_=limS;
        elseif figure_analysis==3
            E=0;ac_t=[0 0 0];
        end
    elseif simulation_continuation==0
        logS=0;E=0;ac_t=[0 0 0];limS_=limS;
    end
    
    visual=struct('cdata', cell(1,iframe),'colormap', cell(1,iframe));


    fig=figure(1);
    set(fig,'Position',[frame_pos(1) frame_pos(2) frame_size(1) frame_size(2)])
    minz=0.0;maxz=1.0;
    [visual(iframe),E,ac_t,logS,limS_]=visual_2d_2018(...
        figure_analysis,c,nx,ny,x_coord,y_coord,minz,maxz,time,...
        frequency_domain,ci_ave,limS_,logS,...
        xspinodal,yspinodal,xbinodal,ybinodal,T,ci,...
        iframe,...
        ac_t,alpha,beta,fc,limSstretcher,...
        nx_grad,ny_grad,x_grad,y_grad,contour_lgd_font_size,...
        nworkers);

end

% --- Enter main loop -----------------------------------------------------

if view_simulation_only==0
    
    fprintf('Enter main cycle\n')
    while 1

% --- Update concentration profile ----------------------------------------

        for i=1:1:nfour
            cp(i)=c(i)+dt*((c(i)-co(i))/dto);
        end
        coo=co;co=c;c=cp;

% --- Apply BC to newly predicted concentration ---------------------------

        % This is unnecessary as the values that must be zeroed should have
        % already been zeroed.
        
% === Enter Newton-Raphson iterations =====================================

        jk=0;
        err=inf;
        fprintf('Entering Newton Raphson iterations\n')

        while err>nr_tol

            toc(counter_init)

% --- Assert reasonable iteration count -----------------------------------

            jk=jk+1;
            if jk>nr_max_iteration
                disp('Too many Newton Raphson iterations')
                jk=-1;
                break
            end

% --- Reset residual vector and Jacobian matrix ---------------------------
            
            %=======CONSTRUCTION===========================
            %=======CONSTRUCTION===========================
            %=======CONSTRUCTION===========================
            %=======CONSTRUCTION===========================
            %=======CONSTRUCTION===========================
            %=======CONSTRUCTION===========================
            %=======CONSTRUCTION===========================
            
            if wetting==1
                c=bc_wetting(c,ll,li,lu,bl,bi,bu,rl,ri,ru,tl,ti,tu,...
                    gl,gb,gr,gt,hl,hb,hr,ht);
                
            end
            
            %=======CONSTRUCTION===========================
            %=======CONSTRUCTION===========================
            %=======CONSTRUCTION===========================
            %=======CONSTRUCTION===========================
            %=======CONSTRUCTION===========================
            %=======CONSTRUCTION===========================

            fprintf('Computing conc...\n')
            sj=get_conc(ne,ney,ny,c,weights);
            fprintf('Computing terms...\n')
            [terms_sf,terms_sj]=terms_wetting_2d(sj,ne,n1,n2,Two_chi_n1,...
                grad_T,coef_T,diff,terms_assist,thermophoresis,wetting);

            fprintf('Computing residual vector...\n')
            sf=sf_wetting_2d(ne,nfour,weights,nx,ny,n,co,dt,dxdy,sj,...
                terms_sf,wTerms,...
                grad_T,wetting);
            
            
            
            fprintf('Computing Jacobian matrix...\n')
            sj=sj_wetting_2d(...
                nnz_,irow,icol,ny,n,weights,dt,coef_T,dxdy,...
                terms_sj,wTerms,communication_reliever,nworkers,...
                nCoreTasks,index_array,grad_T,thermophoresis,...
                wetting,sj_assist);
            
            
%             warning('ending editin')
%             return


% --- Apply boundary conditions -------------------------------------------

            fprintf('Applying bc...\n')
            
%             sf_test=sf;
%             sj_test=sj;
            
            
            warning('IMPROPER BC!')
            [sf,sj]=parallel_wetting_bc(sj,iZeros_wetting,iOnes_wetting,...
                sf,ll,li,lu,bl,bi,bu,rl,ri,ru,tl,ti,tu,...
                wetting,ig);
            
            
            
            
%             return
            
%             ig
%             size(ig)
%             ig(24,1)
%             nnz_
% %             size(ig)
%             4*(nx+ny-1)
%             3*(nx+ny)
%             2*(nx+ny)
            
%             ig
%             temp=1:1:20;
%             [iZeros_wetting temp']
%             
%             size(iOnes_wetting)
%             nnz(iOnes_wetting)
% %             iOnes_wetting'
%             2*(nx+ny)
%             
%             
%             size(iZeros_wetting)
            
%             iZeros;
%             
%             
%             iZeros_wetting
%             size(iZeros_wetting)
%             warning('Clock this')
%             
%             size(iOnes)

%             sj_trial=sj_test;
%             sf_trial=sf_test;
% % %             
%             [sf_,sj_]=parallel_bc(sj,iZeros,iOnes,...
%                 sf,ll,li,lu,bl,bi,bu,rl,ri,ru,tl,ti,tu,...
%                 0);
%             
%             for ijk=1:1:nnz_
%                 gbf1=irow(ijk);
%                 [node,class]=classify_gbf(gbf1);
%                 if class==4
%                     continue
%                 end
%                 sj(ijk)=sj_(ijk);
%             end
%             for ijk=1:1:nfour
%                 [node,class]=classify_gbf(ijk);
%                 if class==4
%                     continue
%                 end
%                 sf(ijk)=sf_(ijk);
%             end
%             
%             if sum(abs(sf_test-sf))~=0 || sum(abs(sj_test-sj))~=0
%                 error('fdasf')
%             end
            
% %             
%             sum1=0;
%             sum2=0;
%             
%             for gbf=1:1:nfour
%                 [node,class]=classify_gbf(gbf);
%                 if class==4
%                     continue
%                 end
%                 sf_trial(gbf)=sf_test(gbf);
%             end
%             
%             for ijk=1:1:nnz_
%                 gbf1=irow(ijk);
%                 [node,class]=classify_gbf(gbf1);
%                 if class==4
%                     continue
%                 end
%                 
%                 sj_trial(ijk)=sj_test(ijk);
%                 
%                 sum1=sum1+abs(sj(ijk));
%                 sum2=sum2+abs(sj_test(ijk));
%                 if abs(sj(ijk))~=abs(sj_test(ijk))
%                     warning('see into prob')
%                     sj(ijk)
%                     sj_test(ijk)
%                     irow(ijk)
%                     icol(ijk)
%                     ijk
%                     error('look here')
%                 end
%             end
%             
%             if sum1-sum2~=0
%                 error('bad')
%             end
% %             
% %             sf=sf_test;
% %             sj=sj_test;
% % %             
%             tester=abs(sf)-abs(sf_test);
%             tester';
%             
%             sf=sf_trial;
%             sj=sj_trial;
%             return

% --- Carry out linear matrix division ------------------------------------

            fprintf('Sparsing matrix...\n')
            sj=sparse(irow,icol,sj);
            fprintf('Carrying out matrix division...\n')
            c_=sj\-sf';
            
            
%             test=zeros(1,nfour);
%             index=0;
%             for i=1:1:n
%                 for j=1:1:4
%                     index=index+1;
%                     test(index)=i;
%                 end
%             end
%             c'+c_;
% %             iZeros
% %             [c_ test'*10^-12]
% %             warning('look here')
%             return
            
% --- Evaluate error ------------------------------------------------------

            fprintf('Assessing current iterations...\n')
            err=sqrt(sum(c_.^2.0))

% --- Assert convergence --------------------------------------------------

            if isnan(err)
                disp('Likely matrix is too close to or equal to singular.')
%                 jk=-2;
                jk=-1;
                break
            elseif ~isreal(c_) || (jk>1 && err>=erro)
                disp('Newton Raphson is diverging or giving imaginary!')
                jk=-1;
                break
            else
                erro=err;
            end

% --- Update the solution -------------------------------------------------

            fprintf('Updating iteration...\n')
            c=c+c_';

        end
        

% === Exit Newton Raphson iteration =======================================

        fprintf('Exit Newton-Raphson iterations\n')

% --- Resolve divergence --------------------------------------------------

        if jk==-2
            disp('Convergence is difficult from this point onward. Simulation will terminate!')
            break
        elseif jk==-1
            bypass=bypass+1;
            c=co;
            co=coo;
            if bypass>bypass_tol
                disp('Too many bypassing!')
                break
            else
                dt=dt*dt_down;
                if dt<dt_min_
                    
                    disp('Time step is too small!')
                    break
                else
                    continue
                end
            end
        elseif bypass~=0
            bypass=0;
        end

% --- Adjust the shift in average concentration due to numerical error ---

        fprintf('Adjusting numerical errors...\n')
        [c,range_check]=c_adjuster_2d(c,ci_ave,n,nfour);
        if range_check==0
            disp('Concentration has reached or exceeded the range between 0 and 1!')
            break
        end
        
        %=============CONSTRUCTING================================
        %=============CONSTRUCTING================================
        %=============CONSTRUCTING================================
        %=============CONSTRUCTING================================
        %=============CONSTRUCTING================================
        %=============CONSTRUCTING================================
        
%         c'
%         warning('hey hey hey')
%         return
        
        if wetting==0
            c=bc(c,ll,li,lu,bl,bi,bu,rl,ri,ru,tl,ti,tu);
        elseif wetting==1
            c=bc_wetting(c,ll,li,lu,bl,bi,bu,rl,ri,ru,tl,ti,tu,...
                gl,gb,gr,gt,hl,hb,hr,ht);
        end
        
        
        %=============CONSTRUCTING================================
        %=============CONSTRUCTING================================
        %=============CONSTRUCTING================================
        %=============CONSTRUCTING================================
        %=============CONSTRUCTING================================
        %=============CONSTRUCTING================================
        %=============CONSTRUCTING================================
        %=============CONSTRUCTING================================

% --- Update time for most recently obtained result -----------------------

        time=time+dt;
        iframe=iframe+1;

% --- Export the data if specified ----------------------------------------

        if export_data==1
            fprintf('Exporting Data...\n')
            dlmwrite([output_path 'time/iteration_' num2str(iframe)],time,'precision',sig_dig);
            dlmwrite([output_path 'concentration/iteration_' num2str(iframe)],c,'precision',sig_dig);
        end
        if export_figure==1
            [logS_export,limS_export,E_updated_export,...
                max_c_export,min_c_export,intensity_specific]=plots_output(...
                fig_output_dim,dpi,axis_font,font_size,...
                output_path,export_fig,...
                c,nx,ny,alpha,beta,fc,ci_ave,...
                iframe,time,logS_export,limS_export,xlen,ylen,limSstretcher,...
                NaN,NaN,NaN,NaN,NaN,NaN,...
                frequency_domain,...
                coef_T,n1,n2,ne,nex,ney,Two_chi_n1,dxdy,weights,grad_T,...
                E_updated_export,max_c_export,min_c_export,diff,...
                export_spec,k_specific,intensity_specific,...
                nworkers);
        end

% --- ANALYZE & PLOT ------------------------------------------------------

        if figure_analysis~=0

            fprintf('Analyzing and preparing for display...\n')
            [visual(iframe),E,ac_t,logS,limS_]=visual_2d_2018(...
                figure_analysis,c,nx,ny,x_coord,y_coord,minz,maxz,time,...
                frequency_domain,ci_ave,limS_,logS,...
                xspinodal,yspinodal,xbinodal,ybinodal,T,ci,...
                iframe,...
                ac_t,alpha,beta,fc,limSstretcher,...
                nx_grad,ny_grad,x_grad,y_grad,contour_lgd_font_size,...
                nworkers);

        end
        

% --- Determine whether to end simulation ---------------------------------

        if iframe==max_frame
            disp('Maximum frame number specified has been reached!')
            break
        elseif time_out<toc(counter_init)
            disp('Time out!')
            break
        elseif time>obs_t
            disp('Observation time goal has been achieved!')
            break
        end

% --- Prepare for next time step ------------------------------------------

        dto=dt;
        if jk>=10
            dt=dt*dt_down;
            if dt<dt_min
                disp('Time step is too small!')
                break
            end
        elseif jk<=nr_up || dt<dt_ideal
            dt=dt*dt_up;
        end
        if time+dt>obs_t
            dt=obs_t-time;
            
            if dt<=0
                disp('Maximum observation time has been reached!')
                break
            end
            
        end

    end
end

% --- Exit main loop ------------------------------------------------------

fprintf('Exit main loop\n')
toc(counter_init)

if export_data==1
    dlmwrite([output_path 'specification/elapsed'],toc(counter_init),'precision',sig_dig);
    dlmwrite([output_path 'specification/status'],1);
end

% --- Export the movie ----------------------------------------------------

if figure_analysis~=0
    try
%     if 1==1
        
        fprintf('Creating and exporting the movie...\n')
        if simulation_continuation==0
            simulation_folder=NaN;
        end
        
        movie_creator_2d(nframes_last_sim,visual,ci,T,nex,ney,ci_ave,ci_fluc,include_fluc,diff,...
            counter_init,output_path,video_method,fps,duration,iframe,...
            2,frame_pos,frame_size,minz,maxz,...
            figure_analysis,nx,ny,x_coord,y_coord,...
            frequency_domain,limS,...
            xspinodal,yspinodal,xbinodal,ybinodal,...
            alpha,beta,fc,limSstretcher,...
            nx_grad,ny_grad,x_grad,y_grad,contour_lgd_font_size,...
            simulation_folder,...
            view_simulation_only,...
            nworkers)

    catch
        warning('Failed to create the video')
        figure_analysis=0;
    end
    
    %--------------------------------
    
    
    if view_simulation_only==1 && export_figure==1
        plots_output_all_iterations(...
            fig_output_dim,dpi,axis_font,font_size,...
            simulation_folder,export_fig,...
            nx,ny,alpha,beta,fc,ci_ave,...
            xlen,ylen,limSstretcher,...
            xspinodal,yspinodal,xbinodal,ybinodal,T,ci,...
            frequency_domain,limS,...
            coef_T,n1,n2,ne,nex,ney,Two_chi_n1,dxdy,weights,grad_T,...
            nframes_last_sim,...
            export_spec,diff,k_specific,...
            nworkers);
    end
    
    %--------------------------------
end

% --- Reorganize the files ------------------------------------------------

if simulation_continuation==1 && figure_analysis~=0
    fprintf('Transferring newly obtained data...\n')
    try
        organize_files(simulation_folder,output_path,nframes_last_sim,iframe,view_simulation_only)
    catch
        warning('Failed to organize the files!')
    end
end

% --- Conclusion ----------------------------------------------------------

fprintf('Simulation has concluded!\n')
toc(counter_init)

end