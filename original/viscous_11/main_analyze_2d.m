function main_analyze_2d(name,distribution_type,T_distribution_spec)

if nargin==0
    clear
end
clc
close all

if ~exist('name','var')
    name='Run_6';
end
if ~exist('distribution_type','var')
    distribution_type=0;
end
if ~exist('T_distribution_spec','var')
    T_distribution_spec=NaN;
end

picked_steps=[1 26 81 83 90 302]; %Make it so only 

profile_6_figs=1; % JUST KEEP THIS EQUAL ONE, OR IT MIGHT NOT WORK!


profile_on=0;
phases_on=0;
% %----Very long time-------
fourier_transform_on=0;
% %-------------------
% view_simulation_on=0;
energy_on=0;
% FT_specific_on=0;
FT_k1_on=0;
FT_k2_on=0;
FT_circular_deliberate=1;
FT_specific=0;
video_create_on=1; % THIS WILL ALSO CREATE LOG PLOT OF K_CIRCULAR, K1, AND K2
zoom_on=0;

if zoom_on==1
    time_range=[0 10^-6];
end


if profile_6_figs==1
    ImageSizeX_profile_6_figs=6;
    ImageSizeY_profile_6_figs=8;

end

ImageSizeX_half_page=6;
ImageSizeY_half_page=4;


%----------Format------------------

% fig_output_dim=[0 0 6 4]; %In PaperUnits Inches
% export_fig=1000;
dpi=150;
axis_font='CMU Serif';
font_size=9;
export_spec=1;
k_specific=0; %Zero will automatically compute characteristic frequency

%--------------------------------

alpha=100;
beta=8;
gamma=8;
fc=20;
nworkers=0;
ifig=1;

%-----------------------------

nworkers=cpu_initializer(nworkers);

folder=[pwd '/Results/' name '/'];

status=dlmread([folder 'specification/status']);
if status==0
    error('This simulation is incomplete')
end

c_folder=[folder 'concentration/iteration_'];
t_folder=[folder 'time/iteration_'];

ci=dlmread([folder 'specification/spec_ci']);

T=dlmread([folder 'specification/spec_T']);


% k_target=get_characteristic_frequency(diff,ci_ave,mean(T),n1,n2,entropy);
if ~exist([folder 'exporting_figures'],'dir')
    mkdir([folder 'exporting_figures'])
end

exporting_folder=[folder 'exporting_figures/'];

if length(picked_steps)~=6
    error('Exactly 6 steps must be chosen')
end

%--------------------------------
spec=dlmread([folder 'specification/spec_general']);

nex=spec(1);
ney=spec(2);
% dt
% obs_t
diff=spec(5);
% include_fluc
% ci_fluc
% nr_tol
% nr_max_iteration
% time_out
% max_frame
entropy=spec(12);
% T_theta
% dt_down
% dt_up
% dt_min
% dt_ideal
% bypass_tol
n1=spec(19);
n2=spec(20);

% --- Assumption ---

xlen=1;
ylen=1;


% %---------------------------------
% 
nx=nex+1;
ny=ney+1;

dx=xlen/nex;
dy=ylen/ney;

dxdy=dx*dy;

ne=nex*ney;

x=linspace(0,xlen,nx);
y=linspace(0,ylen,ny);

ci_ave=mean(ci);

if length(T)==1
    grad_T=0;
    chi=dimensionless_chi(T,entropy);
    Two_chi_n1=2*chi/n1;
    coef_T=diff*T;
elseif length(T)==2
    grad_T=1;
    coef_T=T_characterization(distribution_type,ne,T,...
        T_distribution_spec,nex,ney,xlen,ylen);
    Two_chi_n1=get_two_chi_n1(ne,coef_T,entropy,n1);
end

% size(Two_chi_n1)
% 
% return



weights=weight_adjuster(generateWeights(),dx,dy);

nIterations=determine_nIterations(folder);

[allocation,extra,load_less,load_more]=wise_task_splitter_v2(...
    nIterations,nworkers);

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

if zoom_on==1
    
%     time_range
    
    try
        time_domain=dlmread([exporting_folder 'data_folder/time_domain.dat']);
        F_circular=dlmread([exporting_folder 'data_folder/FT_circular.dat']);
        F_k1=dlmread([exporting_folder 'data_folder/FT_k1.dat']);
        F_k2=dlmread([exporting_folder 'data_folder/FT_k2.dat']);
        energy=dlmread([exporting_folder 'data_folder/energy.dat']);
    catch
        error('Perhaps the files do not exist. Files must be generated first by setting video_create_on=1 and energy_on=1')
    end
    
    
    if time_range(1)>=time_range(2)
        error('Time range must be in increasing order')
    end

    
    if min(time_domain)>time_range(1) || max(time_domain)<time_range(2)
        error('Specified time_range is not within the simulation range')
    end
    
    i1=0;
    while 1==1
        i1=i1+1;
        if time_domain(i1)>=time_range(1)
            break
        end
    end

    for i2=nIterations:-1:i1
        if time_domain(i2)<=time_range(2)
            break
        end
    end
    
    if i1==i2
        error('Selected time_range is too small for any iteration to exist')
    end
    
    %------------------------
    my_fig=figure(ifig);
    energy_homogeneous=compute_homogeneous_energy(...
        ci_ave,diff,coef_T,n1,n2,Two_chi_n1,xlen,ylen,ne,dxdy,grad_T);
    plot(time_domain(i1:1:i2),energy(i1:1:i2)-energy_homogeneous*ones(i2-i1+2,1),'k*')
    xlabel('Time, t');ylabel('Relative energy, \DeltaE')
    grid on
    grid minor
    ytickformat('%,.1f')
    box on
    set(gca,'FontName',axis_font,'FontSize',font_size)
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 ImageSizeX_half_page ImageSizeY_half_page])
    print('-dpng', strcat(exporting_folder, sprintf('energy(t=[%0.3e %0.3e]).png',time_range(1),time_range(2))) , strcat('-r',num2str(dpi)))
    close(my_fig)
    %----------------------------------------------
    
    my_fig=figure(ifig);
    plot(time_domain(i1:1:i2),F_circular(i1:1:i2),'k*')
    xlabel('Time, t');ylabel('log(|| {\fontname{Lucida Calligraphy}F} ||^2)')
    grid on
    grid minor
    ytickformat('%,.1f')
    box on
    set(gca,'FontName',axis_font,'FontSize',font_size)
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 ImageSizeX_half_page ImageSizeY_half_page])
    print('-dpng', strcat(exporting_folder, sprintf('FT_circular(t=[%0.3e %0.3e]).png',time_range(1),time_range(2))) , strcat('-r',num2str(dpi)))
    close(my_fig)
    
    my_fig=figure(ifig);
    plot(time_domain(i1:1:i2),F_k1(i1:1:i2),'k*')
    xlabel('Time, t');ylabel('log(|| {\fontname{Lucida Calligraphy}F} ||^2)')
    grid on
    grid minor
    ytickformat('%,.1f')
    box on
    set(gca,'FontName',axis_font,'FontSize',font_size)
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 ImageSizeX_half_page ImageSizeY_half_page])
    print('-dpng', strcat(exporting_folder, sprintf('FT_k1(t=[%0.3e %0.3e]).png',time_range(1),time_range(2))) , strcat('-r',num2str(dpi)))
    close(my_fig)
    
    my_fig=figure(ifig);
    plot(time_domain(i1:1:i2),F_k2(i1:1:i2),'k*')
    xlabel('Time, t');ylabel('log(|| {\fontname{Lucida Calligraphy}F} ||^2)')
    grid on
    grid minor
    ytickformat('%,.1f')
    box on
    set(gca,'FontName',axis_font,'FontSize',font_size)
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 ImageSizeX_half_page ImageSizeY_half_page])
    print('-dpng', strcat(exporting_folder, sprintf('FT_k2(t=[%0.3e %0.3e]).png',time_range(1),time_range(2))) , strcat('-r',num2str(dpi)))
    close(my_fig)
    
    
end


if FT_specific==1
    
    if ~exist([folder 'exporting_figures/FT'],'dir') || dlmread([folder 'exporting_figures/FT/FT_processing_status'])==0
        error('FT data has not been generated yet. Please generate them first by setting fourier_transform_on==1')
    end
    x_stretcher=1.30;
    left_shifter=0.8;
    right_shifter=1;
    %--------------JUST GUESSING THIS. SHOULD OBTAIN IT FROM ELSEWHERE-----
    fc=20;
    %----------------------------------------------------------------------
    index=0;
    my_fig=figure(ifig);
    %--------
    twopi=2*pi;
    %---------
    
    for iteration=picked_steps
        index=index+1;
        plotter=subplot(3,2,index);
        time=dlmread([t_folder num2str(iteration)]);
        [real_terms,imaginary_terms,resolution]=load_FT_2d(folder,iteration);
        
        k1_domain=linspace(-fc,fc,resolution);
        k2_domain=linspace(-fc,fc,resolution);
        A=real_terms.^2+imaginary_terms.^2;
        
        if mod(resolution,2)==0
            real_resolution=resolution+1;
            
            half_resolution=resolution/2;
            iCenter=half_resolution+1;
            
            A=[A(1:1:half_resolution,1:1:half_resolution) zeros(half_resolution,1) A(1:1:half_resolution,half_resolution+1:1:resolution);
                zeros(1,real_resolution);
                A(half_resolution+1:1:resolution,1:1:half_resolution) zeros(half_resolution,1) A(half_resolution+1:1:resolution,half_resolution+1:1:resolution)];
            k1_domain=[k1_domain(1:1:half_resolution) 0 k1_domain(half_resolution+1:1:resolution)];
            k2_domain=[k2_domain(1:1:half_resolution) 0 k2_domain(half_resolution+1:1:resolution)];
            c=dlmread([c_folder num2str(iteration)]);
            c_nodal=extract_nodal_weights_2D(c,nx,ny);
            c_relative=c_nodal-ci_ave*ones(ny,nx);
            for i=1:1:real_resolution
                A(iCenter,i)=get_mag_2d(c_relative,...
                    k2_domain(i),0,...
                    ny,nx,twopi);
                A(i,iCenter)=get_mag_2d(c_relative,...
                    0,k1_domain(i),...
                    ny,nx,twopi);
            end
        end
        contourf(k1_domain,k2_domain,A','edgecolor','none')
        camlight
        grid on
        grid minor
        colorbar('location','eastoutside');
        colormap hsv
        title(sprintf('t = %d',time),'FontWeight','Normal')
        original=plotter.Position;
        newPos=original;
        newPos(3)=newPos(3)*x_stretcher;
        if mod(index,2)==1
            newPos(1)=newPos(1)*left_shifter;
        else
            newPos(1)=newPos(1)*right_shifter;
        end
        set(plotter,'Position',newPos)
        xlabel('k_1');ylabel('k_2');zlabel('|| {\fontname{Lucida Calligraphy}F} ||^2')
        
    end
    
    set(gca,'FontName',axis_font,'FontSize',font_size)
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 ImageSizeX_profile_6_figs ImageSizeY_profile_6_figs])
    print('-dpng', strcat(exporting_folder,'FT profile.png') , strcat('-r',num2str(dpi)))
    close(my_fig)
    
end

if video_create_on==1
    
    alpha=100;
    beta=4;
    fc=20;
    resolution=[1280 720];
    duration=10;
    fps_max=30;
    quality=50;
    
    frequency_domain=linspace(0,fc,alpha);
    twopi=2*pi;
    half_angle=twopi/(2*beta);
    angle_list=linspace(half_angle,twopi-half_angle,beta);
    
    cos_list=cos(angle_list);
    sin_list=sin(angle_list);
    
    if nnz(cos_list)~=beta || nnz(sin_list)~=beta
        error('Specified beta results in angle that overlap with k1 or k2 or both axes')
    end
    
    [xspinodal,yspinodal,xbinodal,ybinodal]=spinodal_binodal(...
        n1,n2,entropy,1,200,10^-6,0,10,...
        10^-3,0);

    my_fig=figure(ifig);
    set(my_fig,'Position',[0 0 resolution(1) resolution(2)])
    F_circular=zeros(1,nIterations);
    F_k1=zeros(1,nIterations);
    F_k2=zeros(1,nIterations);
    time_domain=zeros(1,nIterations);
    
    video_frames=struct('cdata', cell(1,nIterations),'colormap', cell(1,nIterations));
    
    for iteration=1:1:nIterations
        
        c=dlmread([c_folder num2str(iteration)]);
        time=dlmread([t_folder num2str(iteration)]);
        
        [video_frames(iteration),F_circular,F_k1,F_k2,time_domain]=export_frame(...
            iteration,c,time,x,y,nx,ny,...
            xspinodal,yspinodal,xbinodal,ybinodal,T,ci,...
            xlen,ylen,ci_ave,...
            alpha,beta,nworkers,frequency_domain,cos_list,sin_list,...
            F_circular,F_k1,F_k2,time_domain);
        
    end
    close(my_fig)
    
    video=VideoWriter([folder 'view_simulation.avi']);
    fps=nIterations/duration;
    if fps_max<fps
        fps=fps_max;
    end
    video.FrameRate=fps;
    video.Quality=quality;
    open(video)
    writeVideo(video,video_frames);
    close(video)
    
    %--------------------------output log plot-------------------------
    
    my_fig=figure(ifig);
    plot(time_domain,F_circular,'k*')
    xlabel('Time, t');ylabel('log(|| {\fontname{Lucida Calligraphy}F} ||^2)')
    grid on
    grid minor
    ytickformat('%,.1f')
    box on
    set(gca,'FontName',axis_font,'FontSize',font_size)
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 ImageSizeX_half_page ImageSizeY_half_page])
    print('-dpng', strcat(exporting_folder,'FT_circular.png') , strcat('-r',num2str(dpi)))
    close(my_fig)
    
    my_fig=figure(ifig);
    plot(time_domain,F_k1,'k*')
    xlabel('Time, t');ylabel('log(|| {\fontname{Lucida Calligraphy}F} ||^2)')
    grid on
    grid minor
    ytickformat('%,.1f')
    box on
    set(gca,'FontName',axis_font,'FontSize',font_size)
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 ImageSizeX_half_page ImageSizeY_half_page])
    print('-dpng', strcat(exporting_folder,'FT_k1.png') , strcat('-r',num2str(dpi)))
    close(my_fig)
    
    my_fig=figure(ifig);
    plot(time_domain,F_k2,'k*')
    xlabel('Time, t');ylabel('log(|| {\fontname{Lucida Calligraphy}F} ||^2)')
    grid on
    grid minor
    ytickformat('%,.1f')
    box on
    set(gca,'FontName',axis_font,'FontSize',font_size)
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 ImageSizeX_half_page ImageSizeY_half_page])
    print('-dpng', strcat(exporting_folder,'FT_k2.png') , strcat('-r',num2str(dpi)))
    close(my_fig)
    
    %------------------------------------------------------------------
    
    if ~exist([exporting_folder 'data_folder'],'dir')
        mkdir([exporting_folder 'data_folder'])
        dlmwrite([exporting_folder 'data_folder/time_domain.dat'],time_domain)
    end
    
    dlmwrite([exporting_folder 'data_folder/FT_circular.dat'],F_circular)
    dlmwrite([exporting_folder 'data_folder/FT_k1.dat'],F_k1)
    dlmwrite([exporting_folder 'data_folder/FT_k2.dat'],F_k2)
    
end





%---------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------


if fourier_transform_on==1
    
%     error('fix resolution here')
    f_resolution=128;
    sig_dig=17;
    
    if ~exist([folder 'exporting_figures/FT'],'dir')
        mkdir([folder 'exporting_figures/FT'])
        mkdir([folder 'exporting_figures/FT/real_term'])
        mkdir([folder 'exporting_figures/FT/imaginary_term'])
        dlmwrite([folder 'exporting_figures/FT/FT_processing_status'],0)
        
    end
    if dlmread([folder 'exporting_figures/FT/FT_processing_status'])==0
        
        fprintf('Parallel computing fourier transform now\n')
        export_fourier_transform_data_2d(folder,nworkers,nIterations,...
            fc,nx,ny,xlen,ylen,ci_ave,sig_dig,f_resolution)
        dlmwrite([folder 'exporting_figures/FT/FT_processing_status'],1)
    end
    
end

if FT_k2_on==1
    
    alpha=100;
    fc=20;
    frequency_domain=linspace(0,fc,alpha);
    twopi=2*pi;
    [allocation_k,extra_k,load_less_k,load_more_k]=wise_task_splitter_v2(...
        6,nworkers);
    temp=zeros(nworkers,load_more_k,alpha);
    my_fig=figure(ifig);
    parfor worker=1:1:nworkers
        
        if worker<=extra_k
            temp(worker,:,:)=FT_axial_2_assist(allocation_k(worker,:),...
                load_more_k,load_more_k,...
                picked_steps,ci_ave,c_folder,alpha,nx,ny,...
                frequency_domain,twopi);
        elseif worker>extra_k
            temp(worker,:,:)=FT_axial_2_assist(allocation_k(worker,:),...
                load_less_k,load_more_k,...
                picked_steps,ci_ave,c_folder,alpha,nx,ny,...
                frequency_domain,twopi);
        end
        
    end
    k_array=zeros(6,alpha);
    
    iplot=0;
    for worker=1:1:nworkers
        if worker<=extra_k
            load=load_more_k;
        elseif worker>extra_k
            load=load_less_k;
        end
        for i=1:1:load
            iplot=iplot+1;
            k_array(iplot,:)=temp(worker,i,:);
        end
    end
    for i=1:1:6
        subplot(3,2,i)
        plot(frequency_domain,k_array(i,:),'k')
        grid on
        grid minor
        time=dlmread([t_folder num2str(picked_steps(i))]);
        title(sprintf('t = %d',time),'FontWeight','Normal')
        xlabel('k_2');ylabel('|| {\fontname{Lucida Calligraphy}F} ||^2');
    end
    set(gca,'FontName',axis_font,'FontSize',font_size)
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 ImageSizeX_profile_6_figs ImageSizeY_profile_6_figs])
    print('-dpng', strcat(exporting_folder,'k2.png') , strcat('-r',num2str(dpi)))
    close(my_fig)
end

if FT_k1_on==1
    
    alpha=100;
    fc=20;
    frequency_domain=linspace(0,fc,alpha);
    twopi=2*pi;
    [allocation_k,extra_k,load_less_k,load_more_k]=wise_task_splitter_v2(...
        6,nworkers);
    temp=zeros(nworkers,load_more_k,alpha);
    my_fig=figure(ifig);
    parfor worker=1:1:nworkers
        
        if worker<=extra_k
            temp(worker,:,:)=FT_axial_assist(allocation_k(worker,:),...
                load_more_k,load_more_k,...
                picked_steps,ci_ave,c_folder,alpha,nx,ny,...
                frequency_domain,twopi);
        elseif worker>extra_k
            temp(worker,:,:)=FT_axial_assist(allocation_k(worker,:),...
                load_less_k,load_more_k,...
                picked_steps,ci_ave,c_folder,alpha,nx,ny,...
                frequency_domain,twopi);
        end
        
    end
    k_array=zeros(6,alpha);
    
    iplot=0;
    for worker=1:1:nworkers
        if worker<=extra_k
            load=load_more_k;
        elseif worker>extra_k
            load=load_less_k;
        end
        for i=1:1:load
            iplot=iplot+1;
            k_array(iplot,:)=temp(worker,i,:);
        end
    end
    for i=1:1:6
        subplot(3,2,i)
        plot(frequency_domain,k_array(i,:),'k')
        grid on
        grid minor
        time=dlmread([t_folder num2str(picked_steps(i))]);
        title(sprintf('t = %d',time),'FontWeight','Normal')
        xlabel('k_1');ylabel('|| {\fontname{Lucida Calligraphy}F} ||^2');
    end
    set(gca,'FontName',axis_font,'FontSize',font_size)
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 ImageSizeX_profile_6_figs ImageSizeY_profile_6_figs])
    print('-dpng', strcat(exporting_folder,'k1.png') , strcat('-r',num2str(dpi)))
    close(my_fig)
end

if FT_circular_deliberate==1
    
    %-
    alpha=100;
    beta=4;
    fc=20;
    %-
    twopi=2*pi;
    frequency_domain=linspace(0,fc,alpha);
    
    half_angle=twopi/(2*beta);
    angle_list=linspace(half_angle,twopi-half_angle,beta);
    
    cos_list=cos(angle_list);
    sin_list=sin(angle_list);
    
    if nnz(cos_list)~=beta || nnz(sin_list)~=beta
        error('Specified beta results in angle that overlap with k1 or k2 or both axes')
    end
    F=zeros(1,alpha);
    index=0;
    my_fig=figure(ifig);
    for iteration=picked_steps
        index=index+1;
        subplot(3,2,index)
        c=dlmread([c_folder num2str(iteration)]);
        time=dlmread([t_folder num2str(iteration)]);
        c_nodal=extract_nodal_weights_2D(c,nx,ny);
        c_relative=c_nodal-ci_ave*ones(ny,nx);
        
        [allocation_FT,extra_FT,load_less_FT,load_more_FT]=...
            wise_task_splitter_v2(...
            alpha,nworkers);
        temp=zeros(nworkers,load_more_FT);
        parfor worker=1:1:nworkers
            if worker<=extra_FT
                temp(worker,:)=FT_circular_assist(...
                    allocation_FT(worker,:),load_more_FT,load_more_FT,...
                    frequency_domain,cos_list,sin_list,nx,ny,twopi,beta,...
                    c_relative);
            elseif worker>extra_FT
                temp(worker,:)=FT_circular_assist(...
                    allocation_FT(worker,:),load_less_FT,load_more_FT,...
                    frequency_domain,cos_list,sin_list,nx,ny,twopi,beta,...
                    c_relative);
            end
        end
        index_FT=0;
        for worker=1:1:nworkers
            if worker<=extra_FT
                load=load_more_FT;
            elseif worker>extra_FT
                load=load_less_FT;
            end
            for i=1:1:load
                index_FT=index_FT+1;
                F(index_FT)=temp(worker,i);
            end
        end
        plot(frequency_domain,F,'k')
        grid on
        grid minor
        xlabel('Frequency domain, k')
        ylabel('|| {\fontname{Lucida Calligraphy}F} ||^2')
        title(sprintf('t = %d',time),'FontWeight','Normal')
        title(sprintf('t = %d',time),'FontWeight','Normal')
    end
    set(gca,'FontName',axis_font,'FontSize',font_size)
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 ImageSizeX_profile_6_figs ImageSizeY_profile_6_figs])
    print('-dpng', strcat(exporting_folder,'k_circular.png') , strcat('-r',num2str(dpi)))
    close(my_fig)
end

if energy_on==1
    my_fig=figure(ifig);
    time_domain=zeros(1,nIterations);
    energy=zeros(1,nIterations);
    energy_homogeneous=compute_homogeneous_energy(...
        ci_ave,diff,coef_T,n1,n2,Two_chi_n1,xlen,ylen,ne,dxdy,grad_T);
    
    energy_homogeneous=0;
    
    temp=zeros(nworkers,load_more);
    temp2=zeros(nworkers,load_more);
    
    fprintf('Parallel computing energy...\n')
    parfor worker=1:1:nworkers
        if worker<=extra
            temp(worker,:)=load_energy_assist(allocation(worker,:),...
                load_more,load_more,...
                coef_T,n1,n2,ne,ney,Two_chi_n1,dxdy,...
                ny,weights,diff,grad_T,c_folder);
            temp2(worker,:)=load_time_assist(allocation(worker,:),...
                load_more,load_more,folder);
        elseif worker>extra
            temp(worker,:)=load_energy_assist(allocation(worker,:),...
                load_less,load_more,...
                coef_T,n1,n2,ne,ney,Two_chi_n1,dxdy,...
                ny,weights,diff,grad_T,c_folder);
            temp2(worker,:)=load_time_assist(allocation(worker,:),...
                load_less,load_more,folder);
        end
        
    end
    fprintf('Finished paralllel computing of energy profile\n')
    index=0;
    for worker=1:1:nworkers
        if worker<=extra
            load=load_more;
        elseif worker>extra
            load=load_less;
        end
        for i=1:1:load
            index=index+1;
            time_domain(index)=temp2(worker,i);
            energy(index)=temp(worker,i);
        end
    end
    plot(time_domain,energy-energy_homogeneous*ones(1,nIterations),'k*')
    xlabel('Time, t');ylabel('Relative energy, \DeltaE')
    grid on
    grid minor
    ytickformat('%,.1f')
    box on
    set(gca,'FontName',axis_font,'FontSize',font_size)
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 ImageSizeX_half_page ImageSizeY_half_page])
    print('-dpng', strcat(exporting_folder,'energy.png') , strcat('-r',num2str(dpi)))
    close(my_fig)
    
    if ~exist([exporting_folder 'data_folder'],'dir')
        mkdir([exporting_folder 'data_folder'])
        dlmwrite([exporting_folder 'data_folder/time_domain.dat'],time_domain)
    end
    
    dlmwrite([exporting_folder 'data_folder/energy.dat'],energy)
    
end

%-------------------------------------

if phases_on==1
    index=0;
    my_fig=figure(ifig);
    for iteration=picked_steps
        index=index+1;
        c=dlmread([c_folder num2str(iteration)]);
        time=dlmread([t_folder num2str(iteration)]);
        c_nodal=extract_nodal_weights_2D(c,nx,ny);
        subplot(3,2,index);
        contourf(x,y,c_nodal,[0 ci_ave],'linestyle', 'none');
        colormap([0 0 0;1 1 1])
        xlabel('x');ylabel('y')
        box on
        title(sprintf('t = %d',time),'FontWeight','Normal')
    end
    set(gca,'FontName',axis_font,'FontSize',font_size)
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 ImageSizeX_profile_6_figs ImageSizeY_profile_6_figs])
    print('-dpng', strcat(exporting_folder,'phases.png') , strcat('-r',num2str(dpi)))
    close(my_fig)
end

if profile_on==1
    index=0;
    my_fig=figure(ifig);
    for iteration=picked_steps
        index=index+1;
        c=dlmread([c_folder num2str(iteration)]);
        time=dlmread([t_folder num2str(iteration)]);
        c_nodal=extract_nodal_weights_2D(c,nx,ny);
        subplot(3,2,index);
        surf(x,y,c_nodal,'edgecolor','none')
        xlabel('x');ylabel('y');zlabel('c')
        grid on
        grid minor
        axis([0 xlen 0 ylen 0 1])
        camlight
        title(sprintf('t = %d',time),'FontWeight','Normal')
    end
    set(gca,'FontName',axis_font,'FontSize',font_size)
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 ImageSizeX_profile_6_figs ImageSizeY_profile_6_figs])
    print('-dpng', strcat(exporting_folder,'profile.png') , strcat('-r',num2str(dpi)))
    close(my_fig)
end

end

function [real_terms,imaginary_terms,resolution]=load_FT_2d(folder,iteration)

try
% if 1==1
    real_terms=dlmread([folder 'exporting_figures/FT/real_term/iteration_' num2str(iteration)]);
    imaginary_terms=dlmread([folder 'exporting_figures/FT/imaginary_term/iteration_' num2str(iteration)]);
    [resolution,~]=size(real_terms);
catch
% elseif 1==2
    error('The file loading failed. Please make sure files exist. If not, you must first create them!\nSet fourier_transform_on=1 to obtain the files!')
end

end

function [frame_cap,F_circular,F_k1,F_k2,time_domain]=export_frame(...
    iteration,c,time,x,y,nx,ny,...
    xspinodal,yspinodal,xbinodal,ybinodal,T,ci,...
    xlen,ylen,ci_ave,...
    alpha,beta,nworkers,frequency_domain,cos_list,sin_list,...
    F_circular,F_k1,F_k2,time_domain)

time_domain(iteration)=time;

subplot(2,4,1)
illustrate_simulation(xspinodal,yspinodal,xbinodal,ybinodal,T,ci);

subplot(2,4,2)
c_nodal=extract_nodal_weights_2D(c,nx,ny);
surf(x,y,c_nodal,'edgecolor','none')
xlabel('x');ylabel('y');zlabel('c')
axis([0 xlen 0 ylen 0 1])
colorbar('southoutside')

subplot(2,4,3)
contourf(x,y,c_nodal,[0 ci_ave],'linestyle', 'none');
xlabel('x');ylabel('y')

subplot(2,4,4)
[allocation,extra,load_less,load_more]=wise_task_splitter_v2(...
    alpha,nworkers);
temp=zeros(nworkers,load_more);
% temp1=zeros(nworkers,load_more);
% temp2=zeros(nworers,load_more);
c_relative=c_nodal-ci_ave*ones(ny,nx);
twopi=2*pi;
% warning('replace with parfor here')
parfor worker=1:1:nworkers
    if worker<=extra
        temp(worker,:)=FT_circular_assist(allocation(worker,:),...
            load_more,load_more,...
            frequency_domain,cos_list,sin_list,nx,ny,twopi,beta,...
            c_relative);
%         temp1(worker,:)=FT_axial_assist(keys,load_more,load_more,...
%             picked_steps,ci_ave,c_folder,alpha,nx,ny,...
%             frequency_domain,twopi);
    elseif worker>extra
        temp(worker,:)=FT_circular_assist(allocation(worker,:),...
            load_less,load_more,...
            frequency_domain,cos_list,sin_list,nx,ny,twopi,beta,...
            c_relative);
%         temp1(worker,:)=FT_axial_assist(keys,load_less,load_more,...
%             picked_steps,ci_ave,c_folder,alpha,nx,ny,...
%             frequency_domain,twopi);
    end
end
index_FT=0;
F=zeros(1,alpha);
for worker=1:1:nworkers
    if worker<=extra
        load=load_more;
    elseif worker>extra
        load=load_less;
    end
    for i=1:1:load
        index_FT=index_FT+1;
        F(index_FT)=temp(worker,i);
    end
end
plot(frequency_domain,F)
xlabel('k');ylabel('|| {\fontname{Lucida Calligraphy}F} ||^2')
grid on
grid minor

F_circular(iteration)=log(max(F));

subplot(2,4,5)
for i=1:1:alpha
    F(i)=get_mag_2d(c_relative,...
        frequency_domain(i),0,...
        ny,nx,twopi);
end
plot(frequency_domain,F)
grid on
grid minor
xlabel('k_1');ylabel('|| {\fontname{Lucida Calligraphy}F} ||^2')

F_k1(iteration)=log(max(F));

subplot(2,4,6)
for i=1:1:alpha
    F(i)=get_mag_2d(c_relative,...
        0,frequency_domain(i),...
        ny,nx,twopi);
end
plot(frequency_domain,F)
grid on
grid minor
xlabel('k_2');ylabel('|| {\fontname{Lucida Calligraphy}F} ||^2')

F_k2(iteration)=log(max(F));

subplot(2,4,7)
plot(time_domain(1:1:iteration),F_circular(1:1:iteration),'k*')
xlabel('t');ylabel('log(|| {\fontname{Lucida Calligraphy}F} ||^2)')
grid on
grid minor

subplot(2,4,8)
plot(time_domain(1:1:iteration),F_k1(1:1:iteration),'ko',time_domain(1:1:iteration),F_k2(1:1:iteration),'kx')
xlabel('t');ylabel('log(|| {\fontname{Lucida Calligraphy}F} ||^2)')
grid on
grid minor

str=sprintf('steps=%i, t=%d',iteration,time);
suptitle(str)
frame_cap=getframe(gcf);

end



function export_fourier_transform_data_2d(folder,nworkers,nIterations,...
    fc,nx,ny,xlen,ylen,ci_ave,sig_dig,f_resolution)

% nIterations=1;

% ONLY GOOD FOR EQUAL MESH

if nx~=ny || xlen~=1 || ylen~=1
    error('This is designed specifically for equal spaced cube of volume = 1 unit squared')
end

[allocation,extra,load_less,load_more]=wise_task_splitter_v2(...
    nIterations,nworkers);

% temp=zeros(nworkers,load_more,f_resolution,f_resolution,f_resolution,2);


f_x=linspace(-fc,fc,f_resolution);
f_y=linspace(-fc,fc,f_resolution);
% f_z=linspace(-fc,fc,f_resolution);

twopi=2*pi;

parfor worker=1:1:nworkers %This is parfor
    if worker<=extra
        FT_data_assist_2d(...
            allocation(worker,:),load_more,load_more,nx,ny,f_x,f_y,...
            folder,ci_ave,twopi,f_resolution,sig_dig);
    elseif worker>extra
        FT_data_assist_2d(...
            allocation(worker,:),load_less,load_more,nx,ny,f_x,f_y,...
            folder,ci_ave,twopi,f_resolution,sig_dig);
    end

end

end

function sol=FT_data_assist_2d(...
    keys,load,load_more,nx,ny,...
    f_x,f_y,folder,ci_ave,twopi,...
    f_resolution,sig_dig)

% THIS ONLY WRITES THE FILES NOW, SO IT IS NOT NECESSARY TO OUTPUT SOL
% FIX LATER

% [real,imaginary]=FT_real_imaginary_3d(c,k1,k2,k3,M,N,O,twopi)
sol=zeros(load_more,f_resolution,f_resolution,2);
real_terms=zeros(f_resolution,f_resolution);
imaginary_terms=zeros(f_resolution,f_resolution);

for i=1:1:load
    
    c=dlmread([folder 'concentration/iteration_' num2str(keys(i))]);
    c=extract_nodal_weights_2D(c,nx,ny);
    c=c-ci_ave*ones(ny,nx);
    
    for xth=1:1:f_resolution
        for yth=1:1:f_resolution
                
            [real,imaginary]=FT_real_imaginary_2d(...
                c,f_y(yth),f_x(xth),ny,nx,twopi);

            sol(i,yth,xth,1)=real;
            sol(i,yth,xth,2)=imaginary;
        end
    end
    
    %------------EXPORTING HERE---------------------------
    for xth=1:1:f_resolution
        for yth=1:1:f_resolution
            real_terms(xth,yth)=sol(i,yth,xth,1);
            imaginary_terms(xth,yth)=sol(i,yth,xth,2);
        end
    end
    
    dlmwrite([folder 'exporting_figures/FT/real_term/iteration_' num2str(keys(i))],real_terms,'precision',sig_dig)
    dlmwrite([folder 'exporting_figures/FT/imaginary_term/iteration_' num2str(keys(i))],imaginary_terms,'precision',sig_dig)
    %------------------------------------------------------
    
end

end

function [real,imaginary]=FT_real_imaginary_2d(c,k1,k2,M,N,twopi)
real=0;
imaginary=0;
k1M=twopi*k1/M;
k2N=twopi*k2/N;
for m=1:1:M
    for n=1:1:N
        real=real+c(m,n)*cos(k1M*(m-1)+k2N*(n-1));
        imaginary=imaginary-c(m,n)*sin(k1M*(m-1)+k2N*(n-1));
    end
end
divider=sqrt(M*N);
real=real/divider;
imaginary=imaginary/divider;
end



function sol=FT_axial_2_assist(keys,load,load_more,...
    picked_steps,ci_ave,c_folder,alpha,nx,ny,...
    frequency_domain,twopi)

sol=zeros(load_more,alpha);


for i=1:1:load
    
    c=dlmread([c_folder num2str(picked_steps(keys(i)))]);
    c_nodal=extract_nodal_weights_2D(c,nx,ny);
    c_relative=c_nodal-ci_ave*ones(ny,nx);
    
    for ialpha=1:1:alpha
        sol(i,ialpha)=get_mag_2d(c_relative,...
            frequency_domain(ialpha),0,...
            ny,nx,twopi);
    end
    
end

end

function sol=FT_axial_assist(keys,load,load_more,...
    picked_steps,ci_ave,c_folder,alpha,nx,ny,...
    frequency_domain,twopi)

sol=zeros(load_more,alpha);

for i=1:1:load
    
    c=dlmread([c_folder num2str(picked_steps(keys(i)))]);
    c_nodal=extract_nodal_weights_2D(c,nx,ny);
    c_relative=c_nodal-ci_ave*ones(ny,nx);
    
    for ialpha=1:1:alpha
        sol(i,ialpha)=get_mag_2d(c_relative,...
            0,...
            frequency_domain(ialpha),...
            ny,nx,twopi);
    end
    
end

end

function sol=FT_circular_assist(keys,load,load_more,...
    frequency_domain,cos_list,sin_list,nx,ny,twopi,beta,...
    c_relative)
sol=zeros(1,load_more);
for i=1:1:load
    summation=0;
    for ibeta=1:1:beta
        summation=summation+get_mag_2d(c_relative,...
            frequency_domain(keys(i))*cos_list(ibeta),...
            frequency_domain(keys(i))*sin_list(ibeta),...
            ny,nx,twopi);
    end
    sol(i)=summation/beta;
end

end

function mag=get_mag_2d(c,k1,k2,M,N,twopi)

csum=0;
ssum=0;

k1M=twopi*k1/M;
k2N=twopi*k2/N;

for m=1:1:M
    for n=1:1:N
        csum=csum+c(m,n)*cos(k1M*(m-1)+k2N*(n-1));
        ssum=ssum-c(m,n)*sin(k1M*(m-1)+k2N*(n-1));
    end
end
mag=(csum^2+ssum^2)/(M*N);


% if sqrt(k1^2+k2^2)==6
%     
%     error('afafds')
% end



end

function sol=load_energy_assist(keys,load,load_more,...
    coef_T,n1,n2,ne,ney,Two_chi_n1,dxdy,...
    ny,weights,diff,grad_T,c_folder)

sol=zeros(1,load_more);
for i=1:1:load
    c=dlmread([c_folder num2str(keys(i))]);
    sol(i)=evaluate_total_energy_v2(...
        coef_T,n1,n2,ne,ney,Two_chi_n1,dxdy,...
        ny,c,weights,diff,...
        grad_T);
end

end