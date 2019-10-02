function main(diff,ci,T,n2,thermophoresis,distribution_type,T_distribution_spec)

if nargin==0
    clear
end
clc
close all

% --- Start the timer -----------------------------------------------------

addpath(genpath('./standard_accessory_2d/'));
counter_init=tic;

% --- User defined variables ----------------------------------------------

% --- Specify the mesh ---

nex=30;
ney=29;

% --- Polymer properties ---

if ~exist('diff','var')
    diff=200000;
end
if ~exist('T','var')
%     T=[0.39999999 0.39999998];
%     T=0.31;
%     T=[0.72 0.71];
    
%     T=[0.72 0.7];
%     T=[0.65 0.6];
%     T=[0.31 0.3]*10^-4;
%     T=[0.31 0.3];
    T=[0.3];
    
end
if ~exist('ci','var')
%     ci=0.5;
    ci=sqrt(10)/(1+sqrt(10));
    
%     ci=[ci-0.1 ci+0.1];
%     
%     ci=[0.8 0.75];
    ci=[0.3];
    
end
if ~exist('n2','var')
%     n2=[1 1.001]; %Must be 1 or larger!
    n2=[1];
end


%-------------------------------------
viscosity=0.1;

mobility=1;



%============================
if ~exist('thermophoresis','var')
    thermophoresis=-0;
    distribution_type=10;
    T_distribution_spec=[3,0];
end

% distribution_type==0 Linear temperature gradeint
% distribution_type==3 Cosine wave throughout entire square domain
% distrubution_type==4 Heat focused at the center
% T_distribution_spec=[smaller_radius, larger_radius]
% distribution_type==5 Sinusoidal wave along the edge of x=0. Such
% amplitude decreases linearly to average temperature at specified position
% T_distribution_spec=[wavenumber, extent_of_wave]

n2_distribution_type=1;
n2_distribution_spec=[3,0.0];

time_dependent_temperature=0;

temperature_type=1;
temperature_spec=1.0*10^-6;

% temperature_type=1 Linear variation of temperature at T(2)
% temperature_spec=time, which is duration in which the temperature
% varies

time_dependent_n2=0;

n2_type=3;
n2_spec=[1.0*10^-7];

%============================

n1=1; %Must be equal to 1!
include_fluc=1;
ci_fluc = 10^-6; %Initial fluctuation
entropy=1;
T_theta=1;

chi_a=-0.5;
chi_b=1;

% 
chi_type=1;

% --- Simulation properties ---

%======================
time_out=3600*(60/60)*6;
dt=1.0e-10;
% If too big, simulation will stop in first time step
% General rule: The more it is inside unstable region, the smaller it
% should be
proceed_request=0; % Yes=1, No=0
%======================
max_frame=1000000;
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
obs_t=1; %Observation time
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
export_data=1; % Yes=1, No=0
simulation_continuation=0;
% simulation_continuation=0 if this is a new simulation
% simulation_continuation=1 to continue previously completed simulation.
% This will automatically export data
if simulation_continuation==1
    simulation_folder=[pwd '/Results/Run_40/'];
end
view_simulation_only=0;
% simulation_continuation=1 This will only create a new video for finished
% simulation. Any video in the folder will be deleted and replaced.
if view_simulation_only==1
    simulation_folder=[pwd '/Results/Run_' num2str(39) '/'];
    if figure_analysis==0
        warning('No figure has been specified. We will use figure_analysis==3')
        figure_analysis=3;
    end
end
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
    export_spec=2;
    k_specific=5; %Zero will automatically compute characteristic frequency
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
% validation_2d(diff,ci,T,n1,n2,time_dependent_temperature,...
%     temperature_type,temperature_spec)

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
    if length(n2)==2
        n2_view=n2(2);
    else
        n2_view=n2;
    end
    [xspinodal,yspinodal,xbinodal,ybinodal]=...
        spinodal_binodal_chi_type_1(...
        n1,n2_view,entropy,T_theta,nT,tol,1,nInitialGuessPts,...
        chi_a,chi_b,...
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

% xspinodal'
% 
% error('fdafda')


% --- Setup other global variables ----------------------------------------

%T_theta does not even make sense anymore!!!

%------------------------------------------


[ne,nx,ny,n,nfour,ll,li,lu,bl,bi,bu,tl,ti,tu,rl,ri,ru,...
    fc,frequency_domain,x_coord,y_coord,T,...
    nnz_,irow,icol,weights,ci_ave,dxdy,...
    chi,coef_T,cp,bypass,dto,...
    ...
    ...
    wTerms,iZeros,iOnes,k_specific,...
    coef_n2,wwTerms]=...
    setUpGlobal_2d(nex,ney,fc,alpha,xlen,ylen,...
    ci,T,entropy,diff,n1,n2,...
    communication_reliever,nworkers,...
    dt,distribution_type,T_distribution_spec,...
    export_figure,k_specific,...
    chi_a,chi_b,chi_type,...
    n2_distribution_type,n2_distribution_spec);


if time_dependent_temperature==1
    varying_duration=1;
end


%-----REMOVE COMMENTS TO VISUALIZE TEMPERATURE!-------

coef_M=zeros(ne,3,3,3);
for e=1:1:ne
    for ix=1:1:3
        for iy=1:1:3
            coef_M(e,ix,iy,1)=mobility;
        end
    end
end

size(coef_T)
coef_T(10,3,3,1)
size(coef_T)
size(wTerms)

sfsj_assist=get_sfsj_assist(ne,weights,wTerms,...
    coef_M, coef_T,viscosity,thermophoresis);

size(sfsj_assist)

warning("stop now please emi")

% 
% test=zeros(ney,nex);
% e=0;
% for exth=1:1:nex
%     for eyth=1:1:ney
%         e=e+1;
%         test(eyth,exth)=coef_T(e,3,2,1);
%     end
% end
% 
% % max(max(test))
% % min(min(test))
% 
% coef_T(e,2,2,3);
% 
% surf(test)
% xlabel('x');ylabel('y')
% % 
% warning('ending')
% return
%------------------------------------------------

% terms_assist=assist_global(...
%     grad_T,ne,coef_T,entropy,n1);

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
            dt_down;dt_up;dt_min;dt_ideal;bypass_tol;n1;n2(1)];
        if length(n2)==2
            sf(21)=n2(2);
        end
        dlmwrite([output_path 'specification/spec_general'],sf);
        dlmwrite([output_path 'specification/status'],0);
        dlmwrite([output_path 'specification/spec_T'],T);
        dlmwrite([output_path 'specification/spec_ci'],ci);
        dlmwrite([output_path 'specification/spec_n2'],n2);
    end
end

%--------------Create initial concentration profile----------------------

if simulation_continuation==0
    % Must be careful to apply bc in co.

%     co=generate_co_2D(ci,nx,ny,include_fluc,ci_fluc);
    
    %----------
%     ci=[0.2 0.3];
    co_type=2;
    co_spec=0;
%     
    co=generate_co_anisotropic_2d(ci,co_type,co_spec,n,nx,ny,nfour,ci_fluc,...
        xlen,ylen);
    
%     error('dfadfa')
%     c_nodal=extract_nodal_weights_2D(co,nx,ny);
%     c_relative=c_nodal-ci_ave*ones(ny,nx);
%     
%     [co_fft,xfft,yfft,~]=fftn_new(c_relative,[ylen xlen],[ny,nx]);
%     
%     z1=reshape(real(co_fft),[1,length(xfft)*length(yfft)]);
%     z2=reshape(imag(co_fft),[1,length(xfft)*length(yfft)]);
%     
%     co_fft_data=[max(abs(z1)),mean(z1),std(z1),max(abs(z2)),mean(z2),std(z2)];
    
%     co_fft_data
%     
%     return
    
    
%-----------------
%     co=dlmread([pwd '/Results/i75/concentration/iteration_50']);
%     
% %     co=dlmread([pwd '/Results/10000 0.3 0.3 1/concentration/iteration_59']);
% % % %     
% % % %     
%     co_noise=generate_co_anisotropic_2d(ci,co_type,co_spec,n,nx,ny,nfour,ci_fluc,...
%         xlen,ylen);
%     for i=1:4:nfour-3
%         co(i)=co(i)+co_noise(i)-mean(ci);
%     end
%-------------------
    
%     co_noise(1)
%     return
    
%     
%     co=generate_random_noise(co,ci_ave,co_fft_data,nx,ny,xlen,ylen);
    
%     co_fft_data
% 
%     length(xfft)
%     length(yfft)
% %     
%     return

% plot([0 1],[0 1])
% xlabel('\fontname{CommercialScript BT}F')
% return

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

% warning('ending editin')
% return

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
    [visual(iframe),E,ac_t,logS,limS_]=visual_2d(...
        figure_analysis,c,nx,ny,x_coord,y_coord,minz,maxz,time,...
        frequency_domain,ci_ave,limS_,logS,...
        xspinodal,yspinodal,xbinodal,ybinodal,T,ci,...
        iframe,...
        ac_t,alpha,beta,fc,limSstretcher,...
        nx_grad,ny_grad,x_grad,y_grad,contour_lgd_font_size,...
        nworkers);
    

end


% return



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
            
% --- Adjust time dependent properties ------------------------------------

            if time_dependent_temperature==1 && varying_duration==1
                fprintf('Adjusting time dependent properties...\n')
                [coef_T,Two_chi_n1,terms_assist]=vary_temperature(...
                    T,temperature_type,temperature_spec,time+dt,...
                    ne,nex,ney,xlen,ylen,entropy,n1,...
                    distribution_type,T_distribution_spec);
            end
            
            if time_dependent_n2==1
                
                [coef_n2,n2_new]=vary_n2(n2_type,n2_spec,n2,time+dt,...
                    n2_distribution_type,ne,n2_distribution_spec,...
                    nex,ney,xlen,ylen);
                warning('heyheu')
                n2_new
            end

% --- Reset residual vector and Jacobian matrix ---------------------------

            fprintf('Computing conc...\n')
            sj=get_conc(ne,ney,ny,c,weights);
            
            fprintf('Computing terms...\n')
            [terms_sf,terms_sj]=terms_2d(sj,ne,n1,chi,...
                coef_T,coef_n2,diff,thermophoresis,...
                nworkers);

            
            fprintf('Computing residual vector...\n')
            sf=sf_2d(ne,nfour,weights,nx,ny,n,co,dt,dxdy,sj,...
                terms_sf,wTerms,viscosity,sfsj_assist);
            
%             error("banken")
%             return
            
            fprintf('Computing Jacobian matrix...\n')
            sj=sj_2d(...
                nnz_,irow,icol,ny,n,weights,dt,dxdy,...
                terms_sj,wTerms,nworkers,viscosity,sfsj_assist);
            
            
            warning('jump jump jump')
%             return

% --- Apply boundary conditions -------------------------------------------

            fprintf('Applying bc...\n')
            [sf,sj]=parallel_bc(sj,iZeros,iOnes,...
                sf,ll,li,lu,bl,bi,bu,rl,ri,ru,tl,ti,tu);
            
            
            
            

% --- Carry out linear matrix division ------------------------------------

            fprintf('Sparsing matrix...\n')
            sj=sparse(irow,icol,sj);
            fprintf('Carrying out matrix division...\n')
            c_=sj\-sf';

% --- Evaluate error ------------------------------------------------------

            fprintf('Assessing current iterations...\n')
            err=sqrt(sum(c_.^2.0))

% --- Assert convergence --------------------------------------------------

            if isnan(err)
                disp('Likely matrix is too close to or equal to singular.')
%                 jk=-2; %Reconsider this later============================
                jk=-1;
                break
            elseif imag(err)~=0 || (jk>1 && err>=erro)
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
        
        [c,range_check]=c_adjuster_2d_new(...
            ne,ney,dxdy,ny,c,weights,ci_ave,xlen,ylen,nfour);
        
        if range_check==0
            disp('Concentration has reached or exceeded the range between 0 and 1!')
            break
        end
        c=bc(c,ll,li,lu,bl,bi,bu,rl,ri,ru,tl,ti,tu);

% --- Update time for most recently obtained result -----------------------

        time=time+dt;
        iframe=iframe+1;
        if time_dependent_temperature==1
            if temperature_type==1 && varying_duration==1 && time>=temperature_spec
                varying_duration=0;
            end
        end

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
            [visual(iframe),E,ac_t,logS,limS_]=visual_2d(...
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

%==================================================================
%==================================================================

