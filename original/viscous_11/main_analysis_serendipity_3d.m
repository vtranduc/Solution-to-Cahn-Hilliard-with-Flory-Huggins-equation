function main_analysis_serendipity_3d(name)

if nargin==0
    clear
end
clc
close all

%----------------------------- Specifications

if ~exist('name','var')
    name='5000 0.3 0.3 1 WET(type=1,spec=8,hs=[0 0],gs=[-20 0])'; %Must be 1 or larger!
end

picked_steps=[1 193 244 271 486 1381]; %Make it so only 
profile_6_figs=1; % JUST KEEP THIS EQUAL ONE, OR IT MIGHT NOT WORK!


relative_concentration_on=0;
phases_on=0;
isosurface_on=0;
%----Very long time-------
fourier_transform_on=0;
fourier_transform_specific_iteration_on=0;
%-------------------
energy_on=0;
FT_k1_on=0;
FT_k2_on=0;
FT_k3_on=0;
FT_spherical_deliberate=0;
FT_specific_on=0;
view_simulation_on=1;

fft_filter_on=0;
structure_factor_time_specific_on=0;


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

folder=[pwd '/Results_3d/' name '/'];

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

%--------------------------------
spec=dlmread([folder 'specification/spec_general']);

if length(picked_steps)~=6
    error('There must be exactly 6 steps to be shown')
end
if exist([folder 'exporting_figures/picked_steps.dat'],'file')
    previously_picked=dlmread([folder 'exporting_figures/picked_steps.dat']);
    for i=1:1:6
        if previously_picked(i)~=picked_steps(i)
            error('Conflicted with previously picked steps')
        end
    end
else
    dlmwrite([folder 'exporting_figures/picked_steps.dat'],picked_steps)
end

nex=spec(1);
ney=spec(2);
nez=spec(3);

% dt
% obs_t
diff=spec(6);
% include_fluc
ci_fluc=spec(8);
% tol_nr
% tol_jk
% time_out
% max_frame
entropy=spec(13);
% T_theta
% dt_down
% dt_up
% dt_min
% dt_ideal
% bypass_tol
n1=spec(20);
n2=spec(21);

% --- Assumption ---

xlen=1;
ylen=1;
zlen=1;

%---------------------------------

nx=nex+1;
ny=ney+1;
nz=nez+1;

ne=nex*ney*nez;

x=linspace(0,xlen,nx);
y=linspace(0,ylen,ny);
z=linspace(0,zlen,nz);

ci_ave=mean(ci);

nIterations=determine_nIterations(folder);

[allocation,extra,load_less,load_more]=wise_task_splitter_v2(...
    nIterations,nworkers);


% time_domain=[];
% mag_sq_log=[];

if structure_factor_time_specific_on==1
    
    %-----
    frequencies=[1 2 3 4];
    beta=50;
    gamma=50;
    %---
    
    %==================FIX THIS=============
    nIterations=271;
    [allocation,extra,load_less,load_more]=wise_task_splitter_v2(...
        nIterations,nworkers);
    %---===========================
    
    twopi=2*pi;
    half_theta_angle=twopi/(2*beta);
    theta_list=linspace(half_theta_angle,twopi-half_theta_angle,beta);
    half_gamma_angle=twopi/(2*gamma);
    phi_list=linspace(half_gamma_angle,twopi-half_gamma_angle,gamma);
    phi_list=phi_list(1:1:floor(gamma/2));
    
    npts=length(theta_list)*length(phi_list);
    
    x_list_raw=zeros(1,npts);
    y_list_raw=zeros(1,npts);
    z_list_raw=zeros(1,npts);
    
    index=0;
    for ibeta=1:1:beta
        for igamma=1:1:length(phi_list)
            index=index+1;
            x_list_raw(index)=sin(phi_list(igamma))*cos(theta_list(ibeta));
            y_list_raw(index)=sin(phi_list(igamma))*sin(theta_list(ibeta));
            z_list_raw(index)=cos(phi_list(igamma));
            
        end
    end
    temp=zeros(nworkers,load_more,2);
    n=nx*ny*nz;
    
    my_fig=figure(ifig);
    hold on
    for ifrequency=1:1:4
    %--
        frequency=frequencies(ifrequency);
        
        x_list=x_list_raw*frequency;
        y_list=y_list_raw*frequency;
        z_list=z_list_raw*frequency;
        
        parfor worker=1:1:nworkers
            if worker<=extra

                temp(worker,:,:)=FT_data_time_assist_specific_2d_assist(...
                    allocation(worker,:),...
                    load_more,load_more,folder,...
                    nx,ny,nz,ci_ave,npts,x_list,y_list,z_list,twopi,n);
            elseif worker>extra

                temp(worker,:,:)=FT_data_time_assist_specific_2d_assist(...
                    allocation(worker,:),...
                    load_less,load_more,folder,...
                    nx,ny,nz,ci_ave,npts,x_list,y_list,z_list,twopi,n);
            end
        end
        logged_data=zeros(nIterations,2);
        index=0;
        for worker=1:1:nworkers
            if worker<=extra
                load=load_more;
            elseif worker>extra
                load=load_less;
            end
            for i=1:1:load
                index=index+1;
                logged_data(index,1)=temp(worker,i,1);
                logged_data(index,2)=temp(worker,i,2);
            end
        end
        plot(logged_data(:,2),logged_data(:,1))
    %---
    end
    hold off
    
    grid on;grid minor
    xlabel('Time, t');ylabel('log(|| {\fontname{Lucida Calligraphy}F} ||^2)');
    %=====================================ADJUST THIS EVERY TIME!
    legend('f=1','f=2','f=3','f=4','Location','southeast')
    %====================================
    set(gca,'FontName',axis_font,'FontSize',font_size)
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 ImageSizeX_half_page ImageSizeY_half_page])
    print('-dpng', strcat(exporting_folder,'FT_circular_specific.png') , strcat('-r',num2str(dpi)))
%     close(my_fig)
    
end

if fft_filter_on==1
    
%     if ~exist([folder 'exporting_figures/FT'],'dir') || dlmread([folder 'exporting_figures/FT/FT_processing_status'])==0
%         error('FT data has not been generated yet. Please generate them first by setting fourier_transform_on==1')
%     end
%     
%     iteration=62;
%     twopi=2*pi;

    %=====================
%     iteration=590;
    iteration=271;
    
% %     filter_spec=5;
%     filter_spec=[7 100];
    filters=[0 1.5
        1.5 2.5
        2.5 3.5
        3.5 inf];
    %=======================
    
    filter_type=3; %Keep constant


    X=linspace(0,xlen,nx);
    Y=linspace(0,ylen,ny);
    Z=linspace(0,zlen,nz);

    
    
    time=dlmread([t_folder num2str(iteration)]);

%     %-------------------
    c=dlmread([c_folder num2str(iteration)]);
    c_nodal=extract_abs_result_serendipity(c,nx,ny,nz);
    c_relative=c_nodal-ci_ave*ones(ny,nx,nz);

    [test2,xfft,yfft,zfft]=fftn_new(c_relative,[ylen xlen zlen],[ny,nx,nz]);
    
    
    my_fig=figure(ifig);
    for ifilter=1:1:4
        subplot(2,2,ifilter)
        
        test3=fft_filter_3d(...
            test2,xfft,yfft,zfft,filter_type,filters(ifilter,:));
        
        %-------------------------------------
%         slice(X,Y,Z,real(ifftn(ifftshift(test3))),xlen/2,ylen/2,zlen/2)
%         
%         
%         camlight
%         grid on
%         grid minor
        
        %----
        
        
        face_colors={'red','blue'};
        edge_color='none';
        fig_domain=[min(x) max(x) min(y) max(y) min(z) max(z)];
        view_angle=3;
        %-----------LITERAL-----------
        half_distance=0.001;
        
        
        isosurface_new(X,Y,Z,...
            real(ifftn(ifftshift(test3)))+ci_ave*ones(ny,nx,nz),...
            ci_ave,fig_domain,view_angle,...
            axis_font,font_size,...
            face_colors,edge_color,half_distance)
        
        
%         contourf(X,Y,real(ifftn(ifftshift(test3))),'edgecolor','none')
%         colorbar
        
        %-----------------------------------------
        xlabel('x');ylabel('y');zlabel('z')
        if filters(ifilter,1)==0
            title(sprintf('%d>',filters(ifilter,2)),'FontWeight','Normal')
        elseif filters(ifilter,2)==inf
            title(sprintf('%d<',filters(ifilter,1)),'FontWeight','Normal')
        else
            title(sprintf('(%d,%d)',filters(ifilter,1),filters(ifilter,2)),'FontWeight','Normal')
        end
    end
    
    %----------
%     subplot(2,2,1)
%     title('<1.5')
%     subplot(2,2,2)
%     title('(1.5,2.5)')
%     subplot(2,2,3)
%     title('2.5,3.5')
%     subplot(2,2,4)
%     title('>3.5')
    %-------------------
    
    set(gca,'FontName',axis_font,'FontSize',font_size)
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 ImageSizeX_half_page ImageSizeY_half_page])
    print('-dpng', strcat([exporting_folder,'filter_' num2str(iteration) '.png']) , strcat('-r',num2str(dpi)))
%     close(my_fig)
    
    
end

if phases_on==1
    my_fig=figure(ifig);
    if profile_6_figs
        index=0;
        for iteration=picked_steps
            index=index+1;
            c=dlmread([c_folder num2str(iteration)]);
            time=dlmread([t_folder num2str(iteration)]);
            plotter=subplot(3,2,index);
            obtain_profile(c,nx,ny,nz,x,y,z,xlen,ylen,zlen,...
                ci_ave,ci_fluc,plotter,time);
        end
        set(gca,'FontName',axis_font,'FontSize',font_size)
        set(gcf,'PaperUnits','inches','PaperPosition',[0 0 ImageSizeX_profile_6_figs ImageSizeY_profile_6_figs])
    %     set(gcf,'units','points','position',[0,0,4000,6000])
        print('-dpng', strcat(exporting_folder,'profile.png') , strcat('-r',num2str(dpi)))
    end
    close(my_fig)
end

if isosurface_on==1
    
    face_colors={'red','blue'};
    edge_color='none';
    fig_domain=[min(x) max(x) min(y) max(y) min(z) max(z)];
    view_angle=3;
    %-----------LITERAL-----------
    half_distance=0.001;
    %-----------------------------
    my_fig=figure(ifig);
    if profile_6_figs
        index=0;
        for iteration=picked_steps
            index=index+1;
            c=dlmread([c_folder num2str(iteration)]);
            time=dlmread([t_folder num2str(iteration)]);
            subplot(3,2,index)
            result=extract_abs_result_serendipity(c,nx,ny,nz);
            
            fv = isosurface(x,y,z,result,ci_ave);
            p2 = patch(fv,'visible','off');
            v=p2.Vertices;
            f=p2.Faces;
            n=isonormals(x,y,z,result,p2);
            [nrows,ncols]=size(n);
            v1=zeros(nrows,ncols);
            v2=zeros(nrows,ncols);
            for i=1:1:nrows
                normal_=n(i,:)/(sqrt(n(i,1)^2+n(i,2)^2+n(i,3)^2));
                
                v1(i,:)=v(i,:)+half_distance*normal_;
                v2(i,:)=v(i,:)-half_distance*normal_;
                
            end
            patch('Faces',f,'Vertices',v1,'FaceColor','blue','EdgeColor',edge_color,'FaceColor',char(face_colors(1)))
            patch('Faces',f,'Vertices',v2,'FaceColor','red','EdgeColor',edge_color,'FaceColor',char(face_colors(2)))
            axis(fig_domain)
            view(view_angle)
            camlight 
            lighting gouraud
            box on
            xlabel('x');ylabel('y');zlabel('z')
            title(sprintf('t = %d',time),'FontWeight','Normal')
            hold on
            plot3([0 0],[0 0],[0 1],'k')
            plot3([0 1],[0 0],[1 1],'k')
            plot3([0 0],[0 1],[1 1],'k')
            hold off
            set(gca,'FontName',axis_font,'FontSize',font_size)
            
        end
        set(gcf,'PaperUnits','inches','PaperPosition',[0 0 ImageSizeX_profile_6_figs ImageSizeY_profile_6_figs])
        print('-dpng', strcat(exporting_folder,'isosurface.png') , strcat('-r',num2str(dpi)))
        close(my_fig)
    end
end

if fourier_transform_on==1

    f_resolution=64;
    sig_dig=17;
    
    if ~exist([folder 'exporting_figures/FT'],'dir')
        mkdir([folder 'exporting_figures/FT'])
        mkdir([folder 'exporting_figures/FT/real_term'])
        mkdir([folder 'exporting_figures/FT/imaginary_term'])
        dlmwrite([folder 'exporting_figures/FT/FT_processing_status'],0)
        
    end
    if dlmread([folder 'exporting_figures/FT/FT_processing_status'])==0
        
        fprintf('Parallel computing fourier transform now\n')
        export_fourier_transform_data(folder,nworkers,nIterations,...
            fc,nx,ny,nz,xlen,ylen,zlen,ci_ave,sig_dig,f_resolution)
        
        dlmwrite([folder 'exporting_figures/FT/FT_processing_status'],1)
        
    end
    
end

if fourier_transform_specific_iteration_on==1
    
    f_resolution=64;
    sig_dig=17;
    
    if ~exist([folder 'exporting_figures/FT'],'dir')
        mkdir([folder 'exporting_figures/FT'])
        mkdir([folder 'exporting_figures/FT/real_term'])
        mkdir([folder 'exporting_figures/FT/imaginary_term'])
        dlmwrite([folder 'exporting_figures/FT/FT_processing_status'],0)
        
    end
    if dlmread([folder 'exporting_figures/FT/FT_processing_status'])==0
        
        fprintf('Parallel computing fourier transform now\n')
        export_fourier_transform_data_specific_3d(folder,nworkers,...
            fc,nx,ny,nz,xlen,ylen,zlen,ci_ave,sig_dig,f_resolution,...
            picked_steps)
        
    end
    
%     if dlmread([folder 'exporting_figures/FT/FT_processing_status'])==0
%         
%         fprintf('Parallel computing fourier transform now\n')
%         export_fourier_transform_data_specific_2d(folder,nworkers,...
%             fc,nx,ny,xlen,ylen,ci_ave,sig_dig,f_resolution,picked_steps)
% %         dlmwrite([folder 'exporting_figures/FT/FT_processing_status'],1)
%     end
    
    
end

if view_simulation_on==1
    
    resolution=[1280 720];
    duration=10;
    fps_max=30;
    quality=50;
    
    alpha=64;
    beta=4;
    gamma=4;
    fc=20;
    twopi=2*pi;
    face_colors={'red','blue'};
    edge_color='none';
    fig_domain=[min(x) max(x) min(y) max(y) min(z) max(z)];
    view_angle=3;
    F_spherical=zeros(1,nIterations);
    F_k1=zeros(1,nIterations);
    F_k2=zeros(1,nIterations);
    F_k3=zeros(1,nIterations);
    time_domain=zeros(1,nIterations);
    
    %-----------LITERAL-----------
    half_distance=0.001;
    
    %-----------------------------
    
    my_fig=figure(ifig);
    set(my_fig,'Position',[0 0 resolution(1) resolution(2)])
    
    video_frames=struct('cdata', cell(1,nIterations),'colormap', cell(1,nIterations));
    
    for iteration=1:1:nIterations
        
        c=dlmread([folder 'concentration/iteration_' num2str(iteration)]);
        time=dlmread([folder 'time/iteration_' num2str(iteration)]);
%         [frame_cap,F_spherical,F_k1,F_k2,time_domain]=
        [video_frames(iteration),F_spherical,F_k1,F_k2,F_k3,time_domain]=export_frame(...
            face_colors,edge_color,fig_domain,view_angle,half_distance,...
            iteration,c,time,x,y,z,nx,ny,nz,...
            xlen,ylen,zlen,ci_ave,...
            alpha,beta,gamma,fc,nworkers,...
            F_spherical,F_k1,F_k2,F_k3,time_domain,twopi);
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
    
    %-------------------------
    
    my_fig=figure(ifig);
    plot(time_domain,F_spherical,'k*')
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
    
    my_fig=figure(ifig);
    plot(time_domain,F_k3,'k*')
    xlabel('Time, t');ylabel('log(|| {\fontname{Lucida Calligraphy}F} ||^2)')
    grid on
    grid minor
    ytickformat('%,.1f')
    box on
    set(gca,'FontName',axis_font,'FontSize',font_size)
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 ImageSizeX_half_page ImageSizeY_half_page])
    print('-dpng', strcat(exporting_folder,'FT_k3.png') , strcat('-r',num2str(dpi)))
    close(my_fig)
    
    %-------------------------
    
    if ~exist([exporting_folder 'data_folder'],'dir')
        mkdir([exporting_folder 'data_folder'])
        dlmwrite([exporting_folder 'data_folder/time_domain.dat'],time_domain)
    end
    
    dlmwrite([exporting_folder 'data_folder/FT_spherical.dat'],F_spherical)
    dlmwrite([exporting_folder 'data_folder/FT_k1.dat'],F_k1)
    dlmwrite([exporting_folder 'data_folder/FT_k2.dat'],F_k2)
    dlmwrite([exporting_folder 'data_folder/FT_k3.dat'],F_k3)
    
%     close(my_fig);
    
end

if energy_on==1
    
    if length(T)~=1
        error('This function has only been designed for uniform temperature')
    end
    
    %------Compute energy at uniform condition--------------------
    chi=dimensionless_chi(T,entropy);
    c2=1-ci_ave;
    energy_homogeneous=ci_ave*log(ci_ave)/n1+c2*log(c2)/n2+ci_ave*c2*chi/n1;
    energy_homogeneous=2*diff*T*energy_homogeneous;
    energy_homogeneous=energy_homogeneous*xlen*ylen*zlen;
    
%     return
    
    %--------------------------------------------------------------
    
    my_fig=figure(ifig);
    
    weights=generateWeights_serendipity_3d(xlen,ylen,zlen,nex,ney,nez);
    dxyz=xlen*ylen*zlen/(nex*ney*nez);
    
    
    temp=zeros(nworkers,load_more);
    temp2=zeros(nworkers,load_more);
    
    fprintf('Computing total energy\n')
    
    parfor worker=1:1:nworkers
        if worker<=extra
            temp(worker,:)=energy_assist_serendipity_3d(...
                allocation(worker,:),load_more,load_more,folder,...
                diff,weights,nex,ney,nx,ny,ne,n1,n2,T,entropy,dxyz);
            temp2(worker,:)=load_time_assist(...
                allocation(worker,:),load_more,load_more,folder);
        elseif worker>extra
            temp(worker,:)=energy_assist_serendipity_3d(...
                allocation(worker,:),load_less,load_more,folder,...
                diff,weights,nex,ney,nx,ny,ne,n1,n2,T,entropy,dxyz);
            temp2(worker,:)=load_time_assist(...
                allocation(worker,:),load_less,load_more,folder);
        end
    end
    fprintf('Obtained total energy profile\n')
    energy=zeros(1,nIterations);
    time=zeros(1,nIterations);
    index=0;
    for worker=1:1:nworkers
        if worker<=extra
            load=load_more;
        elseif worker>extra
            load=load_less;
        end
        for i=1:1:load
            index=index+1;
            energy(index)=temp(worker,i);
            time(index)=temp2(worker,i);
        end
    end
    
    plot(time,energy-energy_homogeneous*ones(1,nIterations),'k*')
    axis([0 inf -inf inf])
    xlabel('Time, t');ylabel('Relative energy, \DeltaE')
    grid on
    grid minor
    ytickformat('%,.1f')
    box on
    set(gca,'FontName',axis_font,'FontSize',font_size)
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 ImageSizeX_half_page ImageSizeY_half_page])
    print('-dpng', strcat(folder,'exporting_figures/energy.png') , strcat('-r',num2str(dpi)))
    
    close(my_fig)
    
     if ~exist([exporting_folder 'data_folder'],'dir')
        mkdir([exporting_folder 'data_folder'])
        dlmwrite([exporting_folder 'data_folder/time_domain.dat'],time)
    end
    
    dlmwrite([exporting_folder 'data_folder/energy.dat'],energy)
    
    
end

if relative_concentration_on==1
    
    x_stretcher=1.5;
    left_shifter=0.6;
    right_shifter=1;
    edge_transparency=0.2;
    
    my_fig=figure(ifig);
    index=0;
    for iteration=picked_steps
        index=index+1;
        plotter=subplot(3,2,index);
%         index=index+1;
        c=dlmread([c_folder num2str(iteration)]);
        
        time=dlmread([t_folder num2str(iteration)]);
        result=extract_abs_result_serendipity(c,nx,ny,nz);
        result=result-ci_ave*ones(ny,nx,nz);
        h=slice(x,y,z,result,xlen/2,ylen/2,zlen/2);
        set(h,'edgealpha',edge_transparency)
        hold on
        h=slice([0 xlen/2 xlen],[0 ylen/2 ylen],[0 zlen/2 zlen],zeros(3,3,3),xlen/2,ylen/2,zlen/2);
        set(h,'FaceColor','none','EdgeColor','k','edgealpha',1)
        hold off
%         set(gcf,'edgecolor','none')
        grid on
        grid minor
        colorbar('location','eastoutside');
        caxis([min(min(min(result))) max(max(max(result)))])
        title(sprintf('t = %d',time),'FontWeight','Normal')
        xlabel('x');ylabel('y');zlabel('z')
        
        original=plotter.Position;
        newPos=original;
        newPos(3)=newPos(3)*x_stretcher;
        if mod(index,2)==1
            newPos(1)=newPos(1)*left_shifter;
        else
            newPos(1)=newPos(1)*right_shifter;
        end
        set(plotter,'Position',newPos)
    end
    
    
%     return
    
%     iax = 1; % Or whichever
%     subaxis(4, 6, iax, 'sh', 0.03, 'sv', 0.01, 'padding', 0, 'margin', 0);
    
    set(gca,'FontName',axis_font,'FontSize',font_size)
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 ImageSizeX_profile_6_figs ImageSizeY_profile_6_figs])
    print('-dpng', strcat(exporting_folder,'relative_concentration.png') , strcat('-r',num2str(dpi)))
    close(my_fig)
end

if FT_specific_on==1
    
    x_stretcher=1.5;
    left_shifter=0.6;
    right_shifter=1;
    edge_transparency=0.2;
    my_fig=figure(ifig);
    index=0;
    for iteration=picked_steps
        index=index+1;
        plotter=subplot(3,2,index);
        [real_terms,imaginary_terms,resolution]=load_FT(folder,iteration);
        A=real_terms.^2+imaginary_terms.^2;

        %----------------------------------------------------FC MUST BE
        %READ FROM THE FILE INSTEAD!
        frequency_domain=linspace(-fc,fc,resolution);
        
        
        %--------------------------------------------
        
        if mod(resolution,2)==0
            
            str=sprintf('Computing fourier transform for x,y,z planes for iteration=%i\n',iteration);
            fprintf(str)
%             iCenter=half_resolution+1;
            
            A=zero_inserter(A);
            iCenter=resolution/2+1;
            frequency_domain=[frequency_domain(1:1:iCenter-1) 0 frequency_domain(iCenter:1:resolution)];
            
            c=dlmread([c_folder num2str(iteration)]);
            result=extract_abs_result_serendipity(c,nx,ny,nz);
            result=result-ci_ave*ones(ny,nx,nz);
            
            A=compute_FT_plane(A,nworkers,result,frequency_domain,nx,ny,nz);
            
        end
        
        %-------------------------------------------
        
        h=slice(frequency_domain,frequency_domain,frequency_domain,A,0,0,0);
        set(h,'EdgeColor','none')
        hold on
        h=slice([-fc,0,fc],[-fc,0,fc],[-fc,0,fc],zeros(3,3,3),0,0,0);
        set(h,'FaceColor','none','EdgeColor','k','edgealpha',edge_transparency)
        hold off
        grid on
        grid minor
        colorbar('location','eastoutside');
        colormap hsv
        time=dlmread([t_folder num2str(iteration)]);
        title(sprintf('t = %d',time),'FontWeight','Normal')
        xlabel('k_1');ylabel('k_2');zlabel('k_3')
        original=plotter.Position;
        newPos=original;
        newPos(3)=newPos(3)*x_stretcher;
        if mod(index,2)==1
            newPos(1)=newPos(1)*left_shifter;
        else
            newPos(1)=newPos(1)*right_shifter;
        end
        set(plotter,'Position',newPos)
        %-----------------------------------------------------
    end
    set(gca,'FontName',axis_font,'FontSize',font_size)
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 ImageSizeX_profile_6_figs ImageSizeY_profile_6_figs])
    print('-dpng', strcat(exporting_folder,'amplitude.png') , strcat('-r',num2str(dpi)))
    close(my_fig)
end

%================================================================

if FT_k3_on==1
    resolution=64;
    fc=20;
    my_fig=figure(ifig);
    index=0;
    frequency_domain=linspace(-fc,fc,resolution);
    A=zeros(1,resolution);
    twopi=2*pi;
    for iteration=picked_steps
        index=index+1;
        subplot(3,2,index);
        
        c=dlmread([c_folder num2str(iteration)]);
        
        time=dlmread([t_folder num2str(iteration)]);
        result=extract_abs_result_serendipity(c,nx,ny,nz);
        result=result-ci_ave*ones(ny,nx,nz);
        for j=1:1:resolution
            A(j)=get_mag_3d(result,0,0,...
                frequency_domain(j),...
                ny,nx,nz,twopi);
        end
        plot(frequency_domain,A,'k')
        grid on
        grid minor
        title(sprintf('t = %d',time),'FontWeight','Normal')
        xlabel('k_3');ylabel('|| {\fontname{Lucida Calligraphy}F} ||^2');
    end
    set(gca,'FontName',axis_font,'FontSize',font_size)
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 ImageSizeX_profile_6_figs ImageSizeY_profile_6_figs])
    print('-dpng', strcat(exporting_folder,'k3.png') , strcat('-r',num2str(dpi)))
    close(my_fig)
end

if FT_k1_on==1
    resolution=64;
    fc=20;
    my_fig=figure(ifig);
    index=0;
    frequency_domain=linspace(-fc,fc,resolution);
    A=zeros(1,resolution);
    twopi=2*pi;
    for iteration=picked_steps
        index=index+1;
        subplot(3,2,index);
        
        c=dlmread([c_folder num2str(iteration)]);
        
        time=dlmread([t_folder num2str(iteration)]);
        result=extract_abs_result_serendipity(c,nx,ny,nz);
        result=result-ci_ave*ones(ny,nx,nz);
        for j=1:1:resolution
            A(j)=get_mag_3d(result,0,...
                frequency_domain(j),0,...
                ny,nx,nz,twopi);
        end
        plot(frequency_domain,A,'k')
        grid on
        grid minor
        title(sprintf('t = %d',time),'FontWeight','Normal')
        xlabel('k_1');ylabel('|| {\fontname{Lucida Calligraphy}F} ||^2');
    end
    set(gca,'FontName',axis_font,'FontSize',font_size)
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 ImageSizeX_profile_6_figs ImageSizeY_profile_6_figs])
    print('-dpng', strcat(exporting_folder,'k1.png') , strcat('-r',num2str(dpi)))
    close(my_fig)
end

if FT_k2_on==1
    resolution=64;
    fc=20;
    my_fig=figure(ifig);
    index=0;
    frequency_domain=linspace(-fc,fc,resolution);
    A=zeros(1,resolution);
    twopi=2*pi;
    for iteration=picked_steps
        index=index+1;
        subplot(3,2,index);
        
        c=dlmread([c_folder num2str(iteration)]);
        
        time=dlmread([t_folder num2str(iteration)]);
        result=extract_abs_result_serendipity(c,nx,ny,nz);
        result=result-ci_ave*ones(ny,nx,nz);
        for j=1:1:resolution
            A(j)=get_mag_3d(result,frequency_domain(j),...
                0,0,ny,nx,nz,twopi);
        end
        plot(frequency_domain,A,'k')
        grid on
        grid minor
        title(sprintf('t = %d',time),'FontWeight','Normal')
        xlabel('k_1');ylabel('|| {\fontname{Lucida Calligraphy}F} ||^2');
    end
    set(gca,'FontName',axis_font,'FontSize',font_size)
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 ImageSizeX_profile_6_figs ImageSizeY_profile_6_figs])
    print('-dpng', strcat(exporting_folder,'k2.png') , strcat('-r',num2str(dpi)))
    close(my_fig)
end

%=======================================================================

if FT_spherical_deliberate==1
    alpha=100;
    warning('gamma and beta set to high values')
    beta=50;
    gamma=50;
    fc=20;
    index=0;
    my_fig=figure(ifig);
    for iteration=picked_steps
        index=index+1;
        subplot(3,2,index);
        
        c=dlmread([c_folder num2str(iteration)]);
        
        time=dlmread([t_folder num2str(iteration)]);
        result=extract_abs_result_serendipity(c,nx,ny,nz);
%         result=result-ci_ave*ones(ny,nx,nz);
        [frequency_domain,magnitude]=fourier_analysis_deliberate_3d(...
            result,ci_ave,alpha,beta,gamma,fc,nworkers);
        plot(frequency_domain,magnitude,'k')
        grid on
        grid minor
        title(sprintf('t = %d',time),'FontWeight','Normal')
        xlabel('k');ylabel('|| {\fontname{Lucida Calligraphy}F} ||^2');
    end
    set(gca,'FontName',axis_font,'FontSize',font_size)
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 ImageSizeX_profile_6_figs ImageSizeY_profile_6_figs])
    print('-dpng', strcat(exporting_folder,'k_spherical.png') , strcat('-r',num2str(dpi)))
    close(my_fig)
end



end

%==========================================================================
%==========================================================================
%==========================================================================
%==========================================================================
%==========================================================================
%==========================================================================
%==========================================================================

function export_fourier_transform_data_specific_3d(folder,nworkers,...
    fc,nx,ny,nz,xlen,ylen,zlen,ci_ave,sig_dig,f_resolution,picked_steps)

% nIterations=1;

% ONLY GOOD FOR EQUAL MESH

if nx~=ny || ny~=nz || xlen~=1 || ylen~=1 || zlen~=1
    error('This is designed specifically for equal spaced cube of volume = 1 unit squared')
end


[allocation,extra,load_less,load_more]=wise_task_splitter_v2(...
    6,nworkers);

% temp=zeros(nworkers,load_more,f_resolution,f_resolution,f_resolution,2);


f_x=linspace(-fc,fc,f_resolution);
f_y=linspace(-fc,fc,f_resolution);
f_z=linspace(-fc,fc,f_resolution);

twopi=2*pi;


parfor worker=1:1:nworkers %This is parfor
    if worker<=extra
        FT_data_assist_specific_2d_new(...
            allocation(worker,:),load_more,load_more,nx,ny,nz,f_x,f_y,f_z,...
            folder,ci_ave,twopi,f_resolution,sig_dig,picked_steps);
    elseif worker>extra
        FT_data_assist_specific_2d_new(...
            allocation(worker,:),load_less,load_more,nx,ny,nz,f_x,f_y,f_z,...
            folder,ci_ave,twopi,f_resolution,sig_dig,picked_steps);
    end

end


end

function sol=FT_data_time_assist_specific_2d_assist(keys,load,load_more,...
    folder,nx,ny,nz,ci_ave,npts,x_list,y_list,z_list,twopi,n)

% sol=zeros(load_more,f_resolution,f_resolution,f_resolution,2);
% real_terms=zeros(f_resolution,f_resolution^2);
% imaginary_terms=zeros(f_resolution,f_resolution^2);

sol=zeros(load_more,2);


for i=1:1:load

    c=dlmread([folder 'concentration/iteration_' num2str(keys(i))]);
    c=extract_abs_result_serendipity(c,nx,ny,nz);
    c=c-ci_ave*ones(ny,nx,nz);
    
    %=======================================
%     for xth=1:1:f_resolution
%         for yth=1:1:f_resolution
%              for zth=1:1:f_resolution   
    summation=0;
    for j=1:1:npts
%         [real,imaginary]=FT_real_imaginary_3d(...
%             c,y_list(j),x_list(j),z_list(j),ny,nx,nz,twopi);
        
        summation=summation+get_mag_3d(c,y_list(j),x_list(j),z_list(j),...
            ny,nx,nz,twopi);
        
    end
    
    sol(i,1)=log(summation/npts);
    sol(i,2)=dlmread([folder 'time/iteration_' num2str(keys(i))]);

%     sol(i,yth,xth,zth,1)=real;
%     sol(i,yth,xth,zth,2)=imaginary;
                
%              end
%         end
%     end
    
    %==================================================
    
    %------------EXPORTING HERE---------------------------
%     for zth=1:1:f_resolution
%         index_y=(zth-1)*f_resolution;
%         for xth=1:1:f_resolution
%             for yth=1:1:f_resolution
%                 real_terms(xth,index_y+yth)=sol(i,yth,xth,zth,1);
%                 imaginary_terms(xth,index_y+yth)=sol(i,yth,xth,zth,2);
%             end
%         end
%     end
%     
%     dlmwrite([folder 'exporting_figures/FT/real_term/iteration_' num2str(picked_steps(keys(i)))],real_terms,'precision',sig_dig)
%     dlmwrite([folder 'exporting_figures/FT/imaginary_term/iteration_' num2str(picked_steps(keys(i)))],imaginary_terms,'precision',sig_dig)
    %------------------------------------------------------
    
end

end

function sol=FT_data_assist_specific_2d_new(...
    keys,load,load_more,nx,ny,nz,...
    f_x,f_y,f_z,folder,ci_ave,twopi,...
    f_resolution,sig_dig,picked_steps)

% THIS ONLY WRITES THE FILES NOW, SO IT IS NOT NECESSARY TO OUTPUT SOL
% FIX LATER

% [real,imaginary]=FT_real_imaginary_3d(c,k1,k2,k3,M,N,O,twopi)
sol=zeros(load_more,f_resolution,f_resolution,f_resolution,2);
real_terms=zeros(f_resolution,f_resolution^2);
imaginary_terms=zeros(f_resolution,f_resolution^2);

for i=1:1:load
    
    if exist([folder 'exporting_figures/FT/real_term/iteration_' num2str(picked_steps(keys(i)))],'file')
        return
    end

    c=dlmread([folder 'concentration/iteration_' num2str(picked_steps(keys(i)))]);
    c=extract_abs_result_serendipity(c,nx,ny,nz);
    c=c-ci_ave*ones(ny,nx,nz);
    
    %=======================================
    for xth=1:1:f_resolution
        for yth=1:1:f_resolution
             for zth=1:1:f_resolution   
                [real,imaginary]=FT_real_imaginary_3d(...
                    c,f_y(yth),f_x(xth),f_z(zth),ny,nx,nz,twopi);

                sol(i,yth,xth,zth,1)=real;
                sol(i,yth,xth,zth,2)=imaginary;
             end
        end
    end
    
    %==================================================
    
    %------------EXPORTING HERE---------------------------
    for zth=1:1:f_resolution
        index_y=(zth-1)*f_resolution;
        for xth=1:1:f_resolution
            for yth=1:1:f_resolution
                real_terms(xth,index_y+yth)=sol(i,yth,xth,zth,1);
                imaginary_terms(xth,index_y+yth)=sol(i,yth,xth,zth,2);
            end
        end
    end
    
    dlmwrite([folder 'exporting_figures/FT/real_term/iteration_' num2str(picked_steps(keys(i)))],real_terms,'precision',sig_dig)
    dlmwrite([folder 'exporting_figures/FT/imaginary_term/iteration_' num2str(picked_steps(keys(i)))],imaginary_terms,'precision',sig_dig)
    %------------------------------------------------------
    
end

end



function [frame_cap,F_spherical,F_k1,F_k2,F_k3,time_domain]=export_frame(...
    face_colors,edge_color,fig_domain,view_angle,half_distance,...
    iteration,c,time,x,y,z,nx,ny,nz,...
    xlen,ylen,zlen,ci_ave,...
    alpha,beta,gamma,fc,nworkers,...
    F_spherical,F_k1,F_k2,F_k3,time_domain,twopi)

result=extract_abs_result_serendipity(c,nx,ny,nz);
% c_=result-ci_ave*ones(ny,nx,nz);

subplot(2,4,1)
fv = isosurface(x,y,z,result,ci_ave);
p2 = patch(fv,'visible','off');
v=p2.Vertices;
f=p2.Faces;
n=isonormals(x,y,z,result,p2);
[nrows,ncols]=size(n);
v1=zeros(nrows,ncols);
v2=zeros(nrows,ncols);
for i=1:1:nrows
    normal_=n(i,:)/(sqrt(n(i,1)^2+n(i,2)^2+n(i,3)^2));

    v1(i,:)=v(i,:)+half_distance*normal_;
    v2(i,:)=v(i,:)-half_distance*normal_;

end
patch('Faces',f,'Vertices',v1,'FaceColor','blue','EdgeColor',edge_color,'FaceColor',char(face_colors(1)))
patch('Faces',f,'Vertices',v2,'FaceColor','red','EdgeColor',edge_color,'FaceColor',char(face_colors(2)))
axis(fig_domain)
view(view_angle)
camlight 
lighting gouraud
box on
xlabel('x');ylabel('y');zlabel('z')
% title(sprintf('t = %d',time),'FontWeight','Normal')
hold on
plot3([0 0],[0 0],[0 1],'k')
plot3([0 1],[0 0],[1 1],'k')
plot3([0 0],[0 1],[1 1],'k')
hold off
% set(gca,'FontName',axis_font,'FontSize',font_size)

subplot(2,4,2)
c_relative=result-ci_ave*ones(ny,nx,nz);
mySlice=slice(x,y,z,c_relative,xlen/2,ylen/2,zlen/2);
set(mySlice,'edgecolor','none')
grid on
grid minor
colorbar('location','southoutside');
title(sprintf('t = %d',time),'FontWeight','Normal')
xlabel('x');ylabel('y');zlabel('z')
hold on
h=slice([0 xlen/2 xlen],[0 ylen/2 ylen],[0 zlen/2 zlen],zeros(3,3,3),xlen/2,ylen/2,zlen/2);
set(h,'FaceColor','none','EdgeColor','k','edgealpha',0.5)
hold off
% original=plotter.Position;
% newPos=original;
% newPos(3)=newPos(3)*x_stretcher;
% if mod(index,2)==1
%     newPos(1)=newPos(1)*left_shifter;
% else
%     newPos(1)=newPos(1)*right_shifter;
% end
% set(plotter,'Position',newPos)


subplot(2,4,3)
[frequency_domain,magnitude]=fourier_analysis_deliberate_3d(...
    result,ci_ave,alpha,beta,gamma,fc,nworkers);
plot(frequency_domain,magnitude)
grid on
grid minor
title(sprintf('t = %d',time),'FontWeight','Normal')
xlabel('k');ylabel('|| {\fontname{Lucida Calligraphy}F} ||^2');

subplot(2,4,4)
mag1=zeros(1,alpha);
for j=1:1:alpha
    mag1(j)=get_mag_3d(c_relative,0,...
        frequency_domain(j),0,...
        ny,nx,nz,twopi);
end
plot(frequency_domain,mag1)
grid on
grid minor
title(sprintf('t = %d',time),'FontWeight','Normal')
xlabel('k_1');ylabel('|| {\fontname{Lucida Calligraphy}F} ||^2');

subplot(2,4,5)
mag2=zeros(1,alpha);
for j=1:1:alpha
    mag2(j)=get_mag_3d(c_relative,frequency_domain(j),...
        0,0,...
        ny,nx,nz,twopi);
end
plot(frequency_domain,mag2)
grid on
grid minor
title(sprintf('t = %d',time),'FontWeight','Normal')
xlabel('k_2');ylabel('|| {\fontname{Lucida Calligraphy}F} ||^2');

subplot(2,4,6)
mag3=zeros(1,alpha);
for j=1:1:alpha
    mag3(j)=get_mag_3d(c_relative,0,...
        0,frequency_domain(j),...
        ny,nx,nz,twopi);
end
plot(frequency_domain,mag3)
grid on
grid minor
title(sprintf('t = %d',time),'FontWeight','Normal')
xlabel('k_3');ylabel('|| {\fontname{Lucida Calligraphy}F} ||^2');


subplot(2,4,7)
time_domain(iteration)=time;
F_spherical(iteration)=log(max(magnitude));
subplot(2,4,7)
plot(time_domain(1:1:iteration),F_spherical(1:1:iteration),'k*')
xlabel('t');ylabel('log(|| {\fontname{Lucida Calligraphy}F} ||^2)')
grid on
grid minor

subplot(2,4,8)
F_k1(iteration)=log(max(mag1));
F_k2(iteration)=log(max(mag2));
F_k3(iteration)=log(max(mag3));
plot(time_domain(1:1:iteration),F_k1(1:1:iteration),'ko',time_domain(1:1:iteration),F_k2(1:1:iteration),'kx',time_domain(1:1:iteration),F_k3(1:1:iteration),'k*')
xlabel('t');ylabel('log(|| {\fontname{Lucida Calligraphy}F} ||^2)')
grid on
grid minor

str=sprintf('steps=%i, t=%d',iteration,time);
suptitle(str)
frame_cap=getframe(gcf);

end


function sol=compute_FT_plane(A,nworkers,result,frequency_domain,nx,ny,nz)
sol=A;
[m,~]=size(A);
[allocation,extra,load_less,load_more]=wise_task_splitter_v2(...
    3,nworkers);

temp=zeros(nworkers,load_more,m,m);
twopi=2*pi;

parfor worker=1:1:nworkers
    if worker<=extra
        temp(worker,:,:,:)=compute_FT_plane_assist(allocation(worker,:),...
            load_more,load_more,m,...
            result,frequency_domain,nx,ny,nz,twopi);
    elseif worker>extra
        temp(worker,:,:,:)=compute_FT_plane_assist(allocation(worker,:),...
            load_less,load_more,m,...
            result,frequency_domain,nx,ny,nz,twopi);
    end
end
iCenter=(m-1)/2+1;
index=0;
for worker=1:1:nworkers
    if worker<=extra
        load=load_more;
    elseif worker>extra
        load=load_less;
    end
    for i=1:1:load
        index=index+1;
        if index==1
            for yth=1:1:m
                for zth=1:1:m
                    sol(iCenter,yth,zth)=temp(worker,i,yth,zth);
                end
            end
        elseif index==2
            for xth=1:1:m
                for zth=1:1:m
                    sol(xth,iCenter,zth)=temp(worker,i,xth,zth);
                end
            end
        elseif index==3
            for xth=1:1:m
                for yth=1:1:m
                    sol(xth,yth,iCenter)=temp(worker,i,xth,yth);
                end
            end
        end
    end
    
end


end

function sol=compute_FT_plane_assist(keys,load,load_more,resolution,...
    result,frequency_domain,nx,ny,nz,twopi)

sol=zeros(load_more,resolution,resolution);

for i=1:1:load
    if keys(i)==1
        for yth=1:1:resolution
            for zth=1:1:resolution
                sol(i,yth,zth)=get_mag_3d(result,frequency_domain(yth),...
                    0,frequency_domain(zth),ny,nx,nz,twopi);
            end
        end
    elseif keys(i)==2
        for xth=1:1:resolution
            for zth=1:1:resolution
                sol(i,xth,zth)=get_mag_3d(result,0,...
                    frequency_domain(xth),frequency_domain(zth),ny,nx,nz,twopi);
            end
        end
    elseif keys(i)==3
        for xth=1:1:resolution
            for yth=1:1:resolution
                sol(i,xth,yth)=get_mag_3d(result,frequency_domain(yth),...
                    frequency_domain(xth),0,ny,nx,nz,twopi);
            end
        end
    end
end

end

function sol=zero_inserter(A)

% This function works only if all dimensions are equal and even

[m,~,~]=size(A);
real_resolution=m+1;
half_resolution=m/2;
sol=zeros(real_resolution,real_resolution,real_resolution);
for zth=1:1:half_resolution
    for yth=1:1:half_resolution
        for xth=1:1:half_resolution
            sol(xth,yth,zth)=A(xth,yth,zth);
        end
        for xth=half_resolution+1:1:m
            sol(xth+1,yth,zth)=A(xth,yth,zth);
        end
    end
    for yth=half_resolution+1:1:m
        for xth=1:1:half_resolution
            sol(xth,yth+1,zth)=A(xth,yth,zth);
        end
        for xth=half_resolution+1:1:m
            sol(xth+1,yth+1,zth)=A(xth,yth,zth);
        end
    end
end
for zth=half_resolution+1:1:m
    for yth=1:1:half_resolution
        for xth=1:1:half_resolution
            sol(xth,yth,zth+1)=A(xth,yth,zth);
        end
        for xth=half_resolution+1:1:m
            sol(xth+1,yth,zth+1)=A(xth,yth,zth);
        end
    end
    for yth=half_resolution+1:1:m
        for xth=1:1:half_resolution
            sol(xth,yth+1,zth+1)=A(xth,yth,zth);
        end
        for xth=half_resolution+1:1:m
            sol(xth+1,yth+1,zth+1)=A(xth,yth,zth);
        end
    end
end

end


function [real_terms,imaginary_terms,resolution]=load_FT(folder,iteration)

try
% if 1==1
    real_part=dlmread([folder 'exporting_figures/FT/real_term/iteration_' num2str(iteration)]);
    imaginary_part=dlmread([folder 'exporting_figures/FT/imaginary_term/iteration_' num2str(iteration)]);
    [resolution,~]=size(real_part);
    
    real_terms=zeros(resolution,resolution,resolution);
    imaginary_terms=zeros(resolution,resolution,resolution);
    
    for zth=1:1:resolution
        index_y=(zth-1)*resolution;
        for yth=1:1:resolution
            for xth=1:1:resolution
                real_terms(xth,yth,zth)=real_part(xth,index_y+yth);
                imaginary_terms(xth,yth,zth)=imaginary_part(xth,index_y+yth);
            end
        end
    end
catch
% elseif 1==2
    error('The file loading failed. Please make sure files exist. If not, you must first create them!\nSet fourier_transform_on=1 to obtain the files!')
end

end

function sol=energy_assist_serendipity_3d(...
    keys,load,load_more,folder,...
    diff,weights,nex,ney,nx,ny,ne,n1,n2,T,entropy,dxyz)

sol=zeros(1,load_more);

for i=1:1:load

    c=dlmread([folder 'concentration/iteration_' num2str(keys(i))]);
    sol(i)=evaluate_total_energy_serendipity_3d(...
        c,diff,weights,nex,ney,nx,ny,ne,n1,n2,T,entropy,dxyz);
    
end

end

function export_fourier_transform_data(folder,nworkers,nIterations,...
    fc,nx,ny,nz,xlen,ylen,zlen,ci_ave,sig_dig,f_resolution)

% nIterations=1;

% ONLY GOOD FOR EQUAL MESH

if nx~=ny || ny~=nz || xlen~=1 || ylen~=1 || zlen~=1
    error('This is designed specifically for equal spaced cube of volume = 1 unit squared')
end



[allocation,extra,load_less,load_more]=wise_task_splitter_v2(...
    nIterations,nworkers);

% temp=zeros(nworkers,load_more,f_resolution,f_resolution,f_resolution,2);


f_x=linspace(-fc,fc,f_resolution);
f_y=linspace(-fc,fc,f_resolution);
f_z=linspace(-fc,fc,f_resolution);

twopi=2*pi;


parfor worker=1:1:nworkers %This is parfor
    if worker<=extra
        FT_data_assist(...
            allocation(worker,:),load_more,load_more,nx,ny,nz,f_x,f_y,f_z,...
            folder,ci_ave,twopi,f_resolution,sig_dig);
    elseif worker>extra
        FT_data_assist(...
            allocation(worker,:),load_less,load_more,nx,ny,nz,f_x,f_y,f_z,...
            folder,ci_ave,twopi,f_resolution,sig_dig);
    end

end

end

function sol=FT_data_assist(keys,load,load_more,nx,ny,nz,f_x,f_y,f_z,folder,ci_ave,twopi,f_resolution,sig_dig)

% THIS ONLY WRITES THE FILES NOW, SO IT IS NOT NECESSARY TO OUTPUT SOL
% FIX LATER

% [real,imaginary]=FT_real_imaginary_3d(c,k1,k2,k3,M,N,O,twopi)
sol=zeros(load_more,f_resolution,f_resolution,f_resolution,2);
real_terms=zeros(f_resolution,f_resolution^2);
imaginary_terms=zeros(f_resolution,f_resolution^2);

for i=1:1:load
    
    c=dlmread([folder 'concentration/iteration_' num2str(keys(i))]);
    
%     c=dlmread([folder 'concentration/iteration_' num2str(100)]);
    
    c=extract_abs_result_serendipity(c,nx,ny,nz);
    c=c-ci_ave*ones(ny,nx,nz);
    
    for xth=1:1:f_resolution
        for yth=1:1:f_resolution
            for zth=1:1:f_resolution
                
                [real,imaginary]=FT_real_imaginary_3d(...
                    c,f_y(yth),f_x(xth),f_z(zth),ny,nx,nz,twopi);
                
                sol(i,yth,xth,zth,1)=real;
                sol(i,yth,xth,zth,2)=imaginary;
                
            end
        end
    end
    
    %------------EXPORTING HERE---------------------------
%     index_y=0;
    for zth=1:1:f_resolution
        index_y=(zth-1)*f_resolution;
        for xth=1:1:f_resolution
            for yth=1:1:f_resolution
                real_terms(xth,index_y+yth)=sol(i,yth,xth,zth,1);
                imaginary_terms(xth,index_y+yth)=sol(i,yth,xth,zth,2);
            end
        end
    end
    
    dlmwrite([folder 'exporting_figures/FT/real_term/iteration_' num2str(keys(i))],real_terms,'precision',sig_dig)
    dlmwrite([folder 'exporting_figures/FT/imaginary_term/iteration_' num2str(keys(i))],imaginary_terms,'precision',sig_dig)
    %------------------------------------------------------
    
end

end

function obtain_profile(c,nx,ny,nz,x_coord,y_coord,z_coord,xlen,ylen,zlen,...
    ci_ave,ci_fluc,plotter,time)

result=extract_abs_result_serendipity(c,nx,ny,nz);
s1=slice(x_coord,y_coord,z_coord,result,[0 xlen],[0 ylen],[0 zlen]);
set(s1,'FaceColor','interp')
colormap(plotter,'gray')
caxis([ci_ave-ci_fluc ci_ave+ci_fluc])
xlabel('x');ylabel('y');zlabel('z')
title(sprintf('t=%d',time),'FontWeight','Normal')

end

function isosurface_new(x,y,z,result,...
    ci_ave,fig_domain,view_angle,...
    axis_font,font_size,face_colors,edge_color,...
    half_distance)

fv = isosurface(x,y,z,result,ci_ave);
p2 = patch(fv,'visible','off');
v=p2.Vertices;
f=p2.Faces;
n=isonormals(x,y,z,result,p2);
[nrows,ncols]=size(n);
v1=zeros(nrows,ncols);
v2=zeros(nrows,ncols);
for i=1:1:nrows
    normal_=n(i,:)/(sqrt(n(i,1)^2+n(i,2)^2+n(i,3)^2));

    v1(i,:)=v(i,:)+half_distance*normal_;
    v2(i,:)=v(i,:)-half_distance*normal_;

end
patch('Faces',f,'Vertices',v1,'FaceColor','blue','EdgeColor',edge_color,'FaceColor',char(face_colors(1)))
patch('Faces',f,'Vertices',v2,'FaceColor','red','EdgeColor',edge_color,'FaceColor',char(face_colors(2)))
axis(fig_domain)
view(view_angle)
camlight 
lighting gouraud
box on
xlabel('x');ylabel('y');zlabel('z')
% title(sprintf('t = %d',time),'FontWeight','Normal')
hold on
plot3([0 0],[0 0],[0 1],'k')
plot3([0 1],[0 0],[1 1],'k')
plot3([0 0],[0 1],[1 1],'k')
hold off
set(gca,'FontName',axis_font,'FontSize',font_size)

end