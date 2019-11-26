function main

clear
clc

% --- Start the timer -----------------------------------------------------

counter_init=tic;

% --- User defined variables ----------------------------------------------

nex=100;
ney=100;

dt=1.0e-6; %If too big, simulation will stop in first time step
% General rule: The more it is inside unstable region, the smaller it
% should be

obs_t=10; %Observation time

diff=15000;

%-----------------TEST
A=0.1;
w=2;
%-------------------


% T=0.65;
% T=0.46;
T=0.36;

% ci = 0.6; %Overall initial concentrations
ci=0.3;
% ci=0.7
include_fluc=1;
ci_fluc = 0.001; %Initial fluctuation

nr_tol=1.0e-6;
nr_max_iteration=20;

%======================
time_out=3600*0.5;
max_frame=1000;
%======================

entropy=1;
T_theta=1;

fps=10;

dt_down=0.5;
dt_up=1.2;
dt_min=1.0e-15;
dt_ideal=5.0e-10;
bypass_tol=5;

xlen=1;ylen=1;

alpha=200;
beta=8;

frame_size=[1280 720];
frame_pos=[0 0];

limS=0.0000001; %Guessed value. This will be adjusted as necessary
limShigher=1;
limSlower=0.01;
limSstretcher=10;
limScompressor=0.1;

fc=1000; %Put fc=0 if you wanna use Nyquist frequency

num_str_fac_1d=10;

limS1d=0.01;
limS1dhigher=1;
limS1dlower=0.01;
limSstretcher1d=10;
limScompressor1d=0.1;

nT_phase=1000000000;
T_min_phase=0.001;

n1=1; %Must be equal to 1!
n2=1; %Must be 1 or larger!

proceed_request=0; % Yes=1, No=0
precompute_conc=1; % Yes=1. No=0. Use "No" if the ram is too limited
% INTRODUCE TRADITIONAL METHOD IN WHICH ONLY 1 CORE IS USED, AS
% TRADITIONAL METHOD IS THE MOST EFFECTIVE WITHOUT PARALLEL COMPUTING


compute_with_terms=1;


show_figure=1; % Yes=1. No=0. Selecting yes will also export video.

%Phase diagram specifications
nT=100;
tol=1.0e-6;
% parallel_computing=0;
nInitialGuessPts=100;
guess_accuracy=10^-3;

%Declare amount of workers to be used
nworkers=1;
% nworkers=0 if you wanna use maximum number of cores
% nworkers=-1 if you wanna specify it as user inpput
% nworkers=1 if no parallel computing is desired
% nworkers>1 for parallel computing with specified number of workers

%Decide whether or not to export the data
export_data=0; % Yes=1, No=0

% --- Initialize multi workers if specified -------------------------------

myCluster=parcluster('local');
if nworkers>myCluster.NumWorkers
    warning('nworkers specified exceeds number of workers available.\nThe number of workers will be reduced to maximum available, which is %d.',myCluster.NumWorkers)
elseif nworkers==0
    nworkers=myCluster.NumWorkers;
elseif nworkers==-1
    prompt=['Please specify number of workers to be used. The maximum available is ' num2str(myCluster.NumWorkers) '\n'];
    while 1
        nworkers=input(prompt);
        if length(nworkers)==1 && mod(nworkers,1)==0 && nworkers<=myCluster.NumWorkers && nworkers>=1
            break
        else
            warning('The input is invalid. Please specify a number equal to or lower than %d.',myCluster.NumWorkers)
        end
    end
end
if nworkers~=1
    parallel_computing=1;
    while 1
        try
            delete(gcp('nocreate'))
            parpool('local',nworkers)
            break
        catch
            nworkers=nworkers-1;
            warning('Though resources are physically available, they are not usable at the moment. nworkers will be reduced to %d.',nworkers')
        end
    end
elseif nworkers==1
    parallel_computing=0;
end
if nworkers==1
    fprintf('%d worker will be used for this run.\n',nworkers)
else
    fprintf('%d workers will be used for this run.\n',nworkers)
end

% --- Validate inputs -----------------------------------------------------

%Do validadation of all variables here

% if precompute_conc==1 || precompute_conc==0
%     fprintf('Initializing parallel clusters...')
%     parpool
%     parallel_computing=1;
%     fprintf('Parallel clusters are ready!')
% else
%     parallel_computing=0;
% end

% --- Phase diagram -------------------------------------------------------

%--Test--------------

% nInitialGuessPts=nInitialGuessPts

%-----------------------


fprintf('Generating phase diagram...\n')
[xspinodal,yspinodal,xbinodal,ybinodal]=spinodal_binodal(...
    n1,n2,entropy,T_theta,nT,tol,parallel_computing,nInitialGuessPts,...
    guess_accuracy);
if show_figure==1
    figure(100)
    %--TEST------------
%     subplot(2,2,1)
    %--------------
    illustrate_simulation(xspinodal,yspinodal,xbinodal,ybinodal,T,ci)
end


%-TEST------------------------------

% [cs,f,potential]=free_energy_diagram(n1,n2,T,entropy,T_theta);
% subplot(2,2,2)
% plot(cs,f)
% grid on
% subplot(2,2,3)
% plot(cs,potential)
% grid on
% 
% 
% 
% disp('Finland')
% % return
%-------------------------------


if proceed_request==1
    proceed=input('Proceed? Yes=1+ENTER or Yes=ENTER, No=Any other numerical+ENTER or string+ENTER\n');
    if proceed~=1
        disp('The simulation has been cancelled by the user')
        return
    else
        disp('The simulation is to proceed')
    end
end

% --- Setup other global variables ----------------------------------------

[ne,nx,ny,n,nfour,ll,li,lu,bl,bi,bu,tl,ti,tu,rl,ri,ru,...
    fc,frequency_domain,isf1d,x_coord,y_coord,T,co,...
    nnz_,irow,icol,weights]=...
    setUpGlobal(nex,ney,fc,alpha,num_str_fac_1d,xlen,ylen,...
    T,ci,include_fluc,ci_fluc);

%--------------BYPASSING HERE--------------------------------------------------

% co=ic_circle(ci,include_fluc,ci_fluc,n,nx,ny,nfour,0.4,nex,ney);

% co=dlmread('hoshizora');



co=create_sine_wave_2d(nx,ny,nfour,A,w,x_coord,y_coord,ci);

%-------------------------------------------------------------------

dx=x_coord(2)-x_coord(1);
dy=y_coord(2)-y_coord(1);
weights=weight_adjuster(weights,dx,dy);

if compute_with_terms==1
    wTerms=compute_wTerms_2d(weights);
end

[adjusted_chi,adjusted_diffT]=chi_diffT_adjuster(T,dx,entropy,T_theta,diff);
dxdy=dx*dy;
[iZeros,iOnes]=get_bc_keyVals(nx,ny,nnz_,irow,icol);

if T(1)==T(nx)
    T=T(1);
else
    T=[T(1) T(nx)];
end

cp=zeros(1,nfour);

bypass=0;

dto=inf;

%------Testing--------------------
% co=zeros(1,nfour);
% sf=-3;
% for i=1:1:n
%     sf=sf+4;
%     co(sf)=ci(1);
% end
% co=brownian_noise(co,nx,ny,ci_fluc);
% max(max(co))-min(min(co))
% co';
% % warning('show me here')


% co=dlmread('iteration_1');
% 
% mean(mean(extract_nodal_weights_2D(co,nx,ny)))
% cotest=co;
% itest=-3;
% for hitoshi=1:1:n
%     itest=itest+4;
%     co(itest)=2*ci-co(itest);
% end
% 
% mean(mean(extract_nodal_weights_2D(co,nx,ny)))
% 
% [cotest' co']
% 
% warning('srfd')


%------------------------------

if export_data==1 || show_figure==1
    output_path=folder_setUp();
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


% --- Initialize the first time step --------------------------------------

time=0.0;
c=co;
iframe=1;
if export_data==1
    dlmwrite([output_path 'time/iteration_1'],0);
    dlmwrite([output_path 'concentration/iteration_1'],c);
end

% --- Analyze initial concentration ---------------------------------------

if show_figure==1
    c_nodal=extract_nodal_weights_2D(c,nx,ny);
    ci_ave=mean(mean(c_nodal));

% --- Obtain the structure factor of initial conditions -------------------

% if show_figure==1
    structure_factor=structure_factor_2D(alpha,beta,fc,c_nodal,ci_ave);
    maxS=max(structure_factor);
    logS=[time log(maxS)];
    freq_domain=linspace(0,fc,alpha);
    structure_factor_1d=structure_factor_1D_multiple(c_nodal,...
        freq_domain,isf1d);
    maxS1d_mult=max(structure_factor_1d,[],2);
    maxS1d=max(maxS1d_mult);
    logS1d=log(maxS1d_mult);

% --- Plot the initial condition ------------------------------------------

% if show_figure==1
    fig=figure(1);
    while limS<=maxS
        limS=limS*limSstretcher;
    end
    while limS1d<=maxS1d
        limS1d=limS1d*limSstretcher1d;
    end
    set(fig,'Position',[frame_pos(1) frame_pos(2) frame_size(1) frame_size(2)])
    minz=0.0;maxz=1.0;
    visual(iframe)=visual2D_gradient(c,nx,ny,x_coord,y_coord,minz,maxz,time,...
        frequency_domain,structure_factor,ci_ave,limS,logS,...
        structure_factor_1d,freq_domain,logS1d,limS1d,...
        xspinodal,yspinodal,xbinodal,ybinodal,T,ci);
    
    
%     return
    
end

% --- Enter main loop -----------------------------------------------------
fprintf('Enter main cycle\n')
while 1
    
% --- Update concentration profile ----------------------------------------

    for i=1:1:nfour
        cp(i)=c(i)+dt*((c(i)-co(i))/dto);
    end
    coo=co;
    co=c;
    c=cp;
    
% --- Apply BC to newly predicted concentration ---------------------------

    c=bc(c,ll,li,lu,bl,bi,bu,rl,ri,ru,tl,ti,tu);
    
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
        
        if precompute_conc==1
            fprintf('Computing conc...\n')
            sj=get_conc(ne,ney,ny,c,weights);
            
            %-------------TEST----------------
%             compute_with_terms=0;
%             sj_=sj;
%             terms=compute_terms_2d(sj,ne,n1,n2,adjusted_chi(1,1));
            %--------------------
            
            if compute_with_terms==1
            %--------------------------INCOMPLETE-------------------
                terms=compute_terms_2d(sj,ne,n1,n2,adjusted_chi(1,1)); % MUST MAKE IT AVAILABLE FOR GRADIENT!!!
                fprintf('Computing residual vector...\n')
                sf=test76(...
                    nfour,weights,nx,ny,n,co,n1,n2,dt,dxdy,sj,adjusted_chi,adjusted_diffT,...
                    terms,wTerms); %CHANGE THE NAME
                fprintf('Computing Jacobian matrix...\n')
                sj=test77(nnz_,irow,icol,ny,n,weights,...
                    ney,dt,adjusted_chi,n1,n2,adjusted_diffT,dxdy,sj,...
                    terms,wTerms); %CHANGE THE NAME
            %----------------------------
            elseif compute_with_terms==0
                fprintf('Computing residual vector...\n')
                sf=compute_sf_given_conc(...
                    nfour,weights,nx,ny,n,co,n1,n2,dt,dx,dy,...
                    sj,adjusted_chi,adjusted_diffT); % INPUT DXDY INSTEAD OF DX AND DY
                fprintf('Computing Jacobian matrix...\n')
                sj=compute_sj_given_conc(nnz_,irow,icol,ny,n,weights,...
                    ney,dt,adjusted_chi,n1,n2,adjusted_diffT,dxdy,sj);
            end
        elseif precompute_conc==0
            fprintf('Computing residual vector...\n')
            sf=compute_sf(nfour,weights,nx,ny,n,c,co,n1,n2,dt,dx,dy,...
                adjusted_chi,adjusted_diffT); % INPUT DXDY INSTEAD OF DX AND DY
            fprintf('Computing Jacobian matrix...\n')
            sj=compute_sj(nnz_,irow,icol,ny,n,c,weights,ney,dt,...
                adjusted_chi,n1,n2,adjusted_diffT,dxdy); % 
        end
        
        
        %-------------------
%         toc(counter_init)
%         display('Finland')
        
%         size(conc)
        
%         take=abs(sf-sf_);
%         take'
%         sum(take)
%         take=abs(sj-sj_);
%         
%         sum(take)

%         try
%             ijkl=ijkl+1
%         catch
%             
%             ijkl=0;
%         end
%         
%         fprintf('Computing residual vector...\n')
%         sf_=test76(...
%             nfour,weights,nx,ny,n,co,n1,n2,dt,dx,dy,sj_,adjusted_chi,adjusted_diffT,...
%             terms,wTerms);
%         fprintf('Computing Jacobian matrix...\n')
%         sj_=test77(nnz_,irow,icol,ny,n,weights,...
%             ney,dt,adjusted_chi,n1,n2,adjusted_diffT,dxdy,sj_,...
%             terms,wTerms);
%         
%         sf_=abs(sf-sf_);
%         sj_=abs(sj-sj_);
%         
%         sf_'
%         
%         sum(sf_)
% %         sum(sj_)
%         
%         if ijkl==3
%             return
%         end
%         
%         return

        %--------------------

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

        err=sqrt(sum(c_.^2.0))
        
        %===TEST=======================
        if isnan(err)
            disp('Too close')
            jk=-100;
            break
        end
        %=================================
        
% --- Update the solution -------------------------------------------------
        
        fprintf('Assessing current iterations...\n')
        c=c+c_';
        
% --- Assert convergence --------------------------------------------------

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

% --- Resolve divergence --------------------------------------------------

    %------- TEST -------------------------
    
    if jk==-100
        
        bypass=bypass+1;
        c=co;
        co=coo;
        if bypass>20
            disp('Too many bypassing!')
            break
        else
            dt_up=dt_up*1000;
            dt=dt*dt_up;
%             if bypass==1
%                 dt=0.45;
%             elseif bypass==2
%                 dt=0.8;
%             elseif bypass==3
%                 dt=0.9;
%             elseif bypass==4
%                 dt=0.999;
%             end
            continue
        end
    elseif bypass~=0
        bypass=0;
    end
    
    %---------------------------------------
    
    if jk==-1
        
        bypass=bypass+1;
        c=co;
        co=coo;
        if bypass>bypass_tol
            disp('Too many bypassing!')
            break
        else
            dt=dt*dt_down;
            
            %==TEST===========================
            
%             dlmwrite('hoshizora',c)
%             
%             disp('Exporting')
%             return
            
            %=========================
            
            
            if dt<dt_min
                disp('Time step is too small!')
                break
            else
                continue
            end
        end
    elseif bypass~=0
        bypass=0;
    end
    
    %------Testing--------------------
%     co=zeros(1,nfour);
%     sf=-3;
%     for i=1:1:n
%         sf=sf+4;
%         co(sf)=ci;
%     end
%     co=brownian_noise(co,nx,ny,ci_fluc);
%     max(max(co))-min(min(co))
%     co'
%     warning('show me here')
%     c=brownian_noise(c,nx,ny,0.05);
    %------------------------------

% --- Update time for most recently obtained result -----------------------

    time=time+dt;
    iframe=iframe+1;
    
% --- Exprt the data if specified -----------------------------------------

    if export_data==1
        dlmwrite([output_path 'time/iteration_' num2str(iframe)],time);
        dlmwrite([output_path 'concentration/iteration_' num2str(iframe)],c);
    end
    
% --- Analyze the structure factor ----------------------------------------

    if show_figure==1
        c_nodal=extract_nodal_weights_2D(c,nx,ny);
        structure_factor=structure_factor_2D(alpha,beta,fc,c_nodal,ci_ave);
        maxS=max(structure_factor);
        
        %========================================
        
%         try
%             compare1;
%         catch
%             compare1=inf;
%             compare2=-inf;
%             compare3=inf;
%         end
%         
%         if maxS>compare1
%             warning('ok good')
%             compare1=maxS;
%         else
%             compare1=maxS;
%         end
%         
%         take=min(min(min(c_nodal)));
%         
%         if take<compare2
%             warning('find')
%             compare2=take;
%         else
%             compare2=take;
%         end
%         
%         take=max(max(max(c_nodal)))-take;
%         
%         if take>compare3
%             warning('growth detected!')
%             compare3=take;
%         else
%             compare3=take;
%         end
        
        %========================================
        
        logS(iframe,1:1:2)=[time log(maxS)];
        structure_factor_1d=structure_factor_1D_multiple(c_nodal,...
            freq_domain,isf1d);
        maxS1d_mult=max(structure_factor_1d,[],2);
        maxS1d=max(maxS1d_mult);
        logS1d(1:1:num_str_fac_1d,iframe)=log(maxS1d_mult);
    
% --- PLOT ----------------------------------------------------------------

%     if show_figure==1
        while limS<=maxS
            limS=limS*limSstretcher;
        end
        while limS1d<=maxS1d
            limS1d=limS1d*limSstretcher1d;
        end
        visual(iframe)=visual2D_gradient(c,nx,ny,x_coord,y_coord,minz,maxz,time,...
            frequency_domain,structure_factor,ci_ave,limS,logS,...
            structure_factor_1d,freq_domain,logS1d,limS1d,...
            xspinodal,yspinodal,xbinodal,ybinodal,T,ci);
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
        
        %---------------------------------------% TESTING
    elseif jk==1
        dt_up=dt_up*10;
        dt=dt*dt_up;
        %-------------------------------------
    elseif jk<=3 || dt<dt_ideal
        dt=dt*dt_up;
    end
    if time+dt>obs_t
        dt=obs_t-time;
    end
        
end

% --- Exit main loop ------------------------------------------------------

fprintf('Exit main loop\n')
toc(counter_init)

if export_data==1
    dlmwrite([output_path 'specification/elapsed'],toc(counter_init));
    dlmwrite([output_path 'specification/status'],1);
end

% --- Export the movie ----------------------------------------------------

if show_figure==1
    try
        if length(ci)==1 && length(T)==1
            vid_title=sprintf('2D nonlinear Cahn Hilliard,%dx%d,ci=%1.3f,fluc=%1.4f,diff=%d,T=%0.3f,elapsed time=%1.3fs.avi',nex,ney,ci_ave,ci_fluc*include_fluc,diff,T,toc(counter_init));
        elseif length(ci)==2 && length(T)==1
            vid_title=sprintf('2D nonlinear Cahn Hilliard,%dx%d,ci=[%1.2f,%1.2f],fluc=%1.4f,diff=%d,T=%0.3f,elapsed time=%1.3fs.avi',nex,ney,ci(1),ci(2),ci_fluc*include_fluc,diff,T,toc(counter_init));
        elseif length(ci)==1 && length(T)==2
            vid_title=sprintf('2D nonlinear Cahn Hilliard,%dx%d,ci=%1.3f,fluc=%1.4f,diff=%d,T=[%1.2f,%1.2f],elapsed time=%1.3fs.avi',nex,ney,ci_ave,ci_fluc*include_fluc,diff,T(1),T(2),toc(counter_init));
        elseif length(ci)==2 && length(T)==2
            vid_title=sprintf('2D nonlinear Cahn Hilliard,%dx%d,ci=[%1.2f,%1.2f],fluc=%1.4f,diff=%d,T=[%1.2f,%1.2f],elapsed time=%1.3fs.avi',nex,ney,ci(1),ci(2),ci_fluc*include_fluc,diff,T(1),T(2),toc(counter_init));
        end
        vid_title=[output_path vid_title];
        video=VideoWriter(vid_title, 'Uncompressed AVI');
        video.FrameRate=fps;
        open(video)
        writeVideo(video,visual);
        close(video)

% --- Conclusion ----------------------------------------------------------

        figure(2)
        plot(logS(:,1),logS(:,2),'x')
        grid on
        fig=figure(3);
        set(fig,'Position',[frame_pos(1) frame_pos(2) frame_size(1) frame_size(2)])
        movie(gcf,visual,1)
    catch
        warning('Failed to create the video')
    end
end
toc(counter_init)
end