function thesis_5_multiple_simuls

clear
clc

%=========================== FOLDERS ======================================


%--------------------------------------Done

% name1='E2 500000 0.5 0.39 1 (i30 DQ 500000 0.5 [0.39 0.41] 1 DT=0 d4s[0,0.25]  (i51 500000 0.5 0.39 1))';
% name2='F2 TQ 500000 0.5 0.39 1 (i30 DQ 500000 0.5 [0.39 0.41] 1 DT=1000000 d4s[0,0.25]  (i51 500000 0.5 0.39 1))';
% name3='G2 TQ 500000 0.5 0.39 0.41 1 (i30 DQ 500000 0.5 [0.39 0.41] 1 DT=10000000 d4s[0,0.25] (i51 500000 0.5 0.39 1))';
% name4='H2 TQ 500000 0.5 0.39 1 (i30 H1 500000 0.5 [0.39 0.41] 1 DT=100000000 d4s[0,0.25] (i51 500000 0.5 0.39 1))';
% 
% 
% frames=[1 50 682
%     1 50 682
%     1 50 682
%     1 50 682];
% 
% DTs=-[0,10^6,10^7,10^8];
% 
% dimensions=[2 2 2 2];

%--------------------------------------Done

name1='5000 0.5 0.3 1';
name2='hs=[-0.0005 0.00 ]_gs=[-0.0005 0.0005] GOOD';
name3='5000 0.5 0.3 1 WET(type=2,spec=1;hs=[-10^-3 10^-3];gs=[-10^-3 10^-3])';
name4='5000 0.5 0.3 1 wetting(d5s[0.05],h=0.1,g=0)';


frames=[1 70 158
    6 80 88
    45 88 107
    42 73 121];

DTs=-[0,10^6,10^7,10^8];

dimensions=[3 3 3 3];

cases={'Case 3.A','Case 3.B','Case 3.C','Case 3.D'};


%==========================================================================


names={name1,name2,name3,name4};


dpi=150;
axis_font='CMU Serif';
font_size=9;

ImageSizeX_profile_6_figs=6;
ImageSizeY_profile_6_figs=8;

export_folder=[pwd '/figures/'];

if ~exist(export_folder,'dir')
    mkdir(export_folder);
end

fig1=figure(1);
% fig2=figure(2);

%-Literals and assumptions-
xlen=1;ylen=1;
az=360-45;
el=45;
%--------------------------

for isimul=1:1:4
    
    if dimensions(isimul)==2
    
        name=char(names(isimul));
        folder=[pwd '/Results/' name '/'];
        spec=dlmread([folder 'specification/spec_general']);
        nex=spec(1);ney=spec(2);
        nx=nex+1;ny=ney+1;
        x=linspace(0,xlen,nx);y=linspace(0,ylen,ny);
    
        for iframe=1:1:3

            c=dlmread([folder 'concentration/iteration_' num2str(frames(isimul,iframe))]);
            time=dlmread([folder 'time/iteration_' num2str(frames(isimul,iframe))]);
            c_nodal=extract_nodal_weights_2D(c,nx,ny);


            handle=subplot(4,3,(isimul-1)*3+iframe);
            contourf(x,y,c_nodal,'linestyle', 'none');
            colorbar('Location','eastoutside')
            colormap(gca,'jet')
            xlabel('x');ylabel('y')
            box on
            title(sprintf('Contour plot at\nt = %d',time),'FontWeight','Normal')
            pos=get(handle,'Position');
            set(handle,'Position',pos.*[1 1 1.3 0.8])

        end
        
    elseif dimensions(isimul)==3
        
        name=char(names(isimul));
        folder=[pwd '/Results_3d/' name '/'];

%         status=dlmread([folder 'specification/status']);
%         if status==0
%             error('This simulation is incomplete')
%         end
        c_folder=[folder 'concentration/iteration_'];
        t_folder=[folder 'time/iteration_'];
        
        xlen=1;
        ylen=1;
        zlen=1;
        spec=dlmread([folder 'specification/spec_general']);
        
        nex=spec(1);
        ney=spec(2);
        nez=spec(3);

        nx=nex+1;
        ny=ney+1;
        nz=nez+1;
        
        x=linspace(0,xlen,nx);
        y=linspace(0,ylen,ny);
        z=linspace(0,zlen,nz);
        
        ci=dlmread([folder 'specification/spec_ci']);
        ci_ave=mean(ci);
        
        %-----------------Literal
        half_distance=0.001;
        %--------------------
        
        
        face_colors={'red','blue'};
        edge_color='none';
        fig_domain=[min(x) max(x) min(y) max(y) min(z) max(z)];
        view_angle=3;
        
        for iframe=1:1:3
%             index=index+1;
            c=dlmread([c_folder num2str(frames(isimul,iframe))]);
            time=dlmread([t_folder num2str(frames(isimul,iframe))]);
            handle=subplot(4,3,(isimul-1)*3+iframe);
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
            
            if iframe==1
            
                my_case=char(cases(isimul));

                str=sprintf(my_case);
                h1 = text(-0.7, 1.2,str);
                set(h1, 'rotation', 90);
            end
        
            
            
        end
        
%         if iframe==1
%             
%             my_case=char(cases(isimul));
%             
%             str=sprintf(my_case);
%             h1 = text(-0.0, 0.0,str);
%             set(h1, 'rotation', 0);
%         end
        
        
        
        
    else
        
        error('No such dimension')
        
    end

    
end

set(gca,'FontName',axis_font,'FontSize',font_size)
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 ImageSizeX_profile_6_figs ImageSizeY_profile_6_figs])
print('-dpng', strcat(export_folder,'profile.png') , strcat('-r',num2str(dpi)))


end


%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%-----------------------END OF MAIN PROGRAM--------------------------------
%--------------------------------------------------------------------------
%----THE REST MAY OR MAY NOT BE ESSENTIAL THAT WERE COPIED AND PASTED------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------



function test125_modified(name,distribution_type,T_distribution_spec)

if nargin==0
    clear
end
clc
close all

if ~exist('name','var')
    name='100000 [0.25 0.35] 0.35 1';
end
if ~exist('distribution_type','var')
    distribution_type=0;
end
if ~exist('T_distribution_spec','var')
    T_distribution_spec=NaN;
end

picked_steps=[1 139 156 265 486 1528]; %Make it so only 

profile_6_figs=1; % JUST KEEP THIS EQUAL ONE, OR IT MIGHT NOT WORK!


profile_on=0;
phases_on=0;
% %----Very long time-------
fourier_transform_on=0;
fourier_transform_new_on=0;
fourier_transform_specific_iteration_on=0;
% %-------------------
% view_simulation_on=0;
energy_on=0;
% FT_specific_on=0;
FT_k1_on=0;
FT_k2_on=0;
FT_circular_deliberate=0;
FT_circular_deliberate_specific=0;
FT_specific=0;
video_create_on=0; % THIS WILL ALSO CREATE LOG PLOT OF K_CIRCULAR, K1, AND K2
zoom_on=0;
potential_on=0;

%----------------
fft_filter_on=0;

test2_on=0;

test3_on=0;
%---------------

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

try
    n2=spec(20);
catch
%     warning('The value of n2 could not be read')
%     n2=input('Manually specify n2 here\n')
    
    n2=dlmread([folder 'specification/spec_n2']);
    if length(n2)~=1
        warning('There is n2 variation')
        n2=input('Manually specify n2 here\n');
    end
    
%     n2=1;
end

% warning('FIX HERE')
% n2=1;

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





return

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

if test3_on==1
    
    
    %---Get the instance of zero energy----
    
    %-Pick new condition or current one-
%     n22=1;
%     T2=0.35;
%     grad_T2=0;
%     coef_T2=diff*T2;
%     Two_chi_n12=2*dimensionless_chi(T2,entropy)/n1;
%     %-
% %     coef_T2=coef_T;
% %     n22=n2;
% %     Two_chi_n12=Two_chi_n1;
% %     grad_T2=grad_T;
%     %-----------------------------------------
%     
%     
%     energy_homogeneous=compute_homogeneous_energy(...
%         ci_ave,diff,coef_T2,n1,n22,Two_chi_n12,xlen,ylen,ne,dxdy,grad_T2);
%     
%     
%     energy=dlmread([folder 'exporting_figures/data_folder/energy.dat']);
%     time=dlmread([folder 'exporting_figures/data_folder/time_domain.dat']);
%     
% 
%     relative_energy=energy-energy_homogeneous*ones(1,nIterations);
%     
%     
%     relative_energy(1)
%     
%     plot(time,relative_energy)
%     grid on
%     grid minor
%     for i=1:1:nIterations
%         if relative_energy(i)<=0
%             fprintf('Relative energy is less than or equal to 0 at\niteration=%i\n',i)
%             break
%         end
%     end
    %----Break down structure factor---------------------
%     iteration=590;
%     
% %     iteration=1
%     
%     filter_spec=[5 7];
%     
%     c=dlmread([c_folder num2str(iteration)]);
%     c_nodal=extract_nodal_weights_2D(c,nx,ny);
%     c_relative=c_nodal-ci_ave*ones(ny,nx);
%     
%     
%     subplot(2,2,1)
%     surf(c_nodal,'edgecolor','none')
%     camlight
%     
%     [test2,xfft,yfft,~]=fftn_new(c_relative,[ylen xlen],[ny,nx]);
%     
%     
%     filter_type=3;
%     
%     test3=fft_filter_2d(...
%         test2,xfft,yfft,filter_type,filter_spec);
%     subplot(2,2,2)
%     c_test3=real(ifftn(ifftshift(test3)))+ci_ave*ones(ny,nx);
%     surf(c_test3,'edgecolor','none')
%     camlight
%     title('Within range')
% 
%     filter_type=1;
%     test4=fft_filter_2d(...
%         test2,xfft,yfft,filter_type,filter_spec(1));
%     subplot(2,2,3)
%     c_test4=real(ifftn(ifftshift(test4)))+ci_ave*ones(ny,nx);
%     surf(c_test4,'edgecolor','none')
%     camlight
%     title('Low frequency noise')
%     
%     filter_type=2;
%     test5=fft_filter_2d(...
%         test2,xfft,yfft,filter_type,filter_spec(2));
%     subplot(2,2,4)
%     c_test5=real(ifftn(ifftshift(test5)))+ci_ave*ones(ny,nx);
%     surf(c_test5,'edgecolor','none')
% %     surf(c_test5+c_test4+c_test3-2*ci_ave*ones(ny,nx),'edgecolor','none')
%     camlight
%     title('High frequency noise')
%     
%     e1=evaluate_total_energy_v2(...
%         coef_T,n1,n2,ne,ney,Two_chi_n1,dxdy,...
%         ny,c,weights,diff,...
%         grad_T)-energy_homogeneous;
%     
%     e2=evaluate_total_energy_v2(...
%         coef_T,n1,n2,ne,ney,Two_chi_n1,dxdy,...
%         ny,c_nodal_to_c(c_test3,nx,ny,nex,ney,dx,dy),...
%         weights,diff,grad_T)-energy_homogeneous;
%     
%     e3=evaluate_total_energy_v2(...
%         coef_T,n1,n2,ne,ney,Two_chi_n1,dxdy,...
%         ny,c_nodal_to_c(c_test4,nx,ny,nex,ney,dx,dy),...
%         weights,diff,grad_T)-energy_homogeneous;
%     
%     e4=evaluate_total_energy_v2(...
%         coef_T,n1,n2,ne,ney,Two_chi_n1,dxdy,...
%         ny,c_nodal_to_c(c_test5,nx,ny,nex,ney,dx,dy),...
%         weights,diff,grad_T)-energy_homogeneous;
%     
%     disp('show here')
    
    %----Line of structure factor---------------------

%     iteration=486;
%     c=dlmread([c_folder num2str(iteration)]);
%     c_nodal=extract_nodal_weights_2D(c,nx,ny);
%     c_relative=c_nodal-ci_ave*ones(ny,nx);
%     trend=char_f_trend_x_dir(c_relative,nx,ny,ylen);
%     
%     
%     
%     plot(linspace(0,xlen,nx),trend,'*')
% 
%     grid on;grid minor

    %----Average concentration at a line-----
    
%     -------
    nx_target=1;
    %------
%     nIterations=2;
%     [allocation,extra,load_less,load_more]=wise_task_splitter_v2(...
%         nIterations,nworkers);
    %-------------
    trend=zeros(nIterations,2);
    for iteration=1:1:nIterations
        c=dlmread([c_folder num2str(iteration)]);
        c_nodal=extract_nodal_weights_2D(c,nx,ny);
        time=dlmread([t_folder num2str(iteration)]);
        summation=0;
        for ijk=1:1:ny
            summation=summation+c_nodal(ijk,nx_target);
        end
        trend(iteration,1)=summation/ny;
        trend(iteration,2)=time;
    end
    my_fig=figure(ifig);
    plot(trend(:,2),trend(:,1),'k')
    xlabel('Time, t');ylabel('Concentration, c')
    grid on;grid minor
    set(gca,'FontName',axis_font,'FontSize',font_size)
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 ImageSizeX_half_page ImageSizeY_half_page])
    print('-dpng', strcat(exporting_folder,'c_last_one.png') , strcat('-r',num2str(dpi)))
    close(my_fig)
    
end

if test2_on==1
    
    iteration=50;
    
    time=dlmread([t_folder num2str(iteration)]);
    twopi=2*pi;

%     %-------------------
    c=dlmread([c_folder num2str(iteration)]);
    c_nodal=extract_nodal_weights_2D(c,nx,ny);
    c_relative=c_nodal-ci_ave*ones(ny,nx);
    
    [real_terms,imaginary_terms,resolution,k1_domain,k2_domain]=load_FT_w_axial_2d(...
        folder,iteration,fc,twopi,nx,ny,ci_ave);
    
    
%     phase=get_phase_specific(real_terms,imaginary_terms,resolution);

    phase=angle(real_terms+sqrt(-1)*imaginary_terms);
    
    
    resolution
    
%     k=2^12;
%     phase=abs(fftn(phase));
%     
%     max(max(phase))
%     min(min(phase))
%     
%     phase
    warning('afdaf')
    max(max(phase))
    phase(1,1)
    
    contourf(phase','edgecolor','none')
    
%     axis([0 0.001 0 0.001])
    
    max(max(phase))
    
    [a,b]=max(phase);
    [c,d]=max(a);
    
    [b(d) d]
    
    
    
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
    iteration=486;
    
% %     filter_spec=5;
%     filter_spec=[7 100];
    filters=[0 1.5
        1.5 3
        3 8
        8 inf];
    
    filters=[0 3
        3 8
        8 12
        12 inf];
    %=======================
    
    filter_type=3; %Keep constant


    X=linspace(0,xlen,nx);
    Y=linspace(0,ylen,ny);

    
    
    time=dlmread([t_folder num2str(iteration)]);

%     %-------------------
    c=dlmread([c_folder num2str(iteration)]);
    c_nodal=extract_nodal_weights_2D(c,nx,ny);
    c_relative=c_nodal-ci_ave*ones(ny,nx);
%     subplot(2,2,1)
%     surf(c_relative,'edgecolor','none')
%     camlight
    [test2,xfft,yfft,~]=fftn_new(c_relative,[ylen xlen],[ny,nx]);
%     test2=fft_filter_2d(...
%         test2,xfft,yfft,filter_type,filter_spec);
%     subplot(2,2,2)
%     contourf(xfft,yfft,(abs(test2)),'edgecolor','none')
%     axis([-fc fc -fc fc])
%     subplot(2,2,3)
%     surf(real(ifftn(ifftshift(test2))),'edgecolor','none')
%     camlight
    
    my_fig=figure(ifig);
    
     %---EXTRA----- Variance
     trend=zeros(4,nx);
     
     %---------------------
     
    for ifilter=1:1:4
        subplot(2,2,ifilter)
        test3=fft_filter_2d(...
            test2,xfft,yfft,filter_type,filters(ifilter,:));
        
        %-------------------------------------
        surf(X,Y,real(ifftn(ifftshift(test3))),'edgecolor','none')
        camlight
        grid on
        grid minor
        
        %---EXTRA----- Variance
        my_data=real(ifftn(ifftshift(test3)));
        for xth=1:1:nx
            summation=0;
            for ijk=1:1:ny
                summation=summation+my_data(ijk,xth);
            end
            average=summation/ny;
            summation=0;
            for ijk=1:1:ny
                summation=summation+(my_data(ijk,xth)-average)^2;
            end
            trend(ifilter,xth)=summation/ny;
        end
        %----------------------
        
        
%         contourf(X,Y,real(ifftn(ifftshift(test3))),'edgecolor','none')
%         colorbar
        
        %-----------------------------------------
        xlabel('x');ylabel('y');zlabel('c_r')
        if filters(ifilter,1)==0
            title(sprintf('%d>',filters(ifilter,2)),'FontWeight','Normal')
        elseif filters(ifilter,2)==inf
            title(sprintf('%d<',filters(ifilter,1)),'FontWeight','Normal')
        else
            title(sprintf('(%d,%d)',filters(ifilter,1),filters(ifilter,2)),'FontWeight','Normal')
        end
    end
    
    %----------
%     subplot(2,2,3)
%     title('(17,18.5)')
%     subplot(2,2,4)
%     title('18.5<')
    %-------------------
    
    set(gca,'FontName',axis_font,'FontSize',font_size)
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 ImageSizeX_half_page ImageSizeY_half_page])
    print('-dpng', strcat([exporting_folder,'filter_' num2str(iteration) '.png']) , strcat('-r',num2str(dpi)))
%     close(my_fig)
    
%     subplot(2,2,4)
%     surf(X,Y,real(ifftn(ifftshift(test2)))./c_relative,'edgecolor','none')
%     axis([0 xlen 0 ylen 0 1])
%     camlight

    
    
    %%---EXTRA----- Variance
    figure(555555)
    plot(X,trend(1,:),X,trend(2,:),X,trend(3,:),X,trend(4,:))
    xlabel('x');ylabel('Concentration variance, \sigma^2')
    grid on;grid minor;
    %====FIX EVERTY TIME
    warning('fix this every time')
    legend('(r_1,r_2)=(0,3)','(r_1,r_2)=(3,8)','(r_1,r_2)=(8,12)','(r_1,r_2)=(12,\infty)','Location','northwest')
    %====
    set(gca,'FontName',axis_font,'FontSize',font_size)
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 ImageSizeX_half_page ImageSizeY_half_page])
    print('-dpng', strcat([exporting_folder,'variance_' num2str(iteration) '.png']) , strcat('-r',num2str(dpi)))
    
%     0 3
%         3 8
%         8 12
%         12 inf
    %--------------------

    
end

%=======================================================

if potential_on==1
    
    %-----
%     iteration=1;
    %--------
    

    
    X=linspace(dx/2,xlen-dx/2,nex);
    Y=linspace(dy/2,ylen-dy/2,ney);
    
    for i=1:1:6
        
        c=dlmread([c_folder num2str(picked_steps(i))]);
        potential=evaluate_potential(c,nex,ney,ny,n1,n2,Two_chi_n1,...
            weights,grad_T,coef_T);
        subplot(3,2,i)
        
        potential=potential+100000*ones(ney,nex);
        potential=log(potential);
        
        surf(X,Y,potential,'edgecolor','none')
        xlabel('x');ylabel('y');zlabel('\mu_2-\mu_1')
        grid on
        grid minor
        camlight
        time=dlmread([t_folder num2str(picked_steps(i))]);
        title(sprintf('t = %d',time),'FontWeight','Normal')
        
    end
    
end

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
    
    if ~exist([folder 'exporting_figures/FT'],'dir')
        error('FT data has not been generated yet. Please generate them first by setting fourier_transform_on==1')
    end
    if dlmread([folder 'exporting_figures/FT/FT_processing_status'])==0
        for i=1:1:6
            if ~exist([folder 'exporting_figures/FT/real_term/iteration_' num2str(picked_steps(i))],'file')
                error('FT data has not been generated yet. Please generate them first by setting fourier_transform_on==1')
            end
        end
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

if profile_on==1 || phases_on==1 || FT_specific==1
    if exist([exporting_folder 'data_folder/picked_steps.dat'],'file')
        previous=dlmread([exporting_folder 'data_folder/picked_steps.dat']);
        for i=1:1:6
            if previous(i)~=picked_steps(i)
                error('Conflicted with previously picked step')
            end
        end
    end
    dlmwrite([exporting_folder 'data_folder/picked_steps.dat'],picked_steps)
end

if FT_circular_deliberate_specific==1
    
    %--------
    beta=8;
    frequencies=[1 2.1 4.1 6];
    %-------------=--FIX
%     nIterations=10;
%     [allocation,extra,load_less,load_more]=wise_task_splitter_v2(...
%         nIterations,nworkers);
    %-------------------------
    my_fig=figure(ifig);

    twopi=2*pi;
    half_angle=twopi/(2*beta);
    angle_list=linspace(half_angle,twopi-half_angle,beta);
    
    cos_list_raw=cos(angle_list);
    sin_list_raw=sin(angle_list);
    
    temp=zeros(nworkers,load_more,2);
    
    
    hold on
    for ifreq=1:1:4
        
        frequency=frequencies(ifreq);
        
        cos_list=cos_list_raw*frequency;
        sin_list=sin_list_raw*frequency;
    
    
        parfor worker=1:1:nworkers
            if worker<=extra
                temp(worker,:,:)=FT_circular_deliberate_specific_assist(...
                    allocation(worker,:),load_more,load_more,...
                    nx,ny,c_folder,t_folder,ci_ave,cos_list,sin_list,twopi,beta);
            elseif worker>extra
                temp(worker,:,:)=FT_circular_deliberate_specific_assist(...
                    allocation(worker,:),load_less,load_more,...
                    nx,ny,c_folder,t_folder,ci_ave,cos_list,sin_list,twopi,beta);
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
        
        
        
    end
    hold off
    
    %=========FIX EVERY TIME========
    
    warning('warn about legend')
    legend('f=1','f=2.1','f=4.1','f=6','Location','northeast')
    
    %================================
    
    grid on;grid minor
    xlabel('k');ylabel('|| {\fontname{Lucida Calligraphy}F} ||^2');
    set(gca,'FontName',axis_font,'FontSize',font_size)
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 ImageSizeX_half_page ImageSizeY_half_page])
    print('-dpng', strcat(exporting_folder,'FT_circular_specific.png') , strcat('-r',num2str(dpi)))
    close(my_fig)
    
end

if fourier_transform_specific_iteration_on==1
    
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
        export_fourier_transform_data_specific_2d(folder,nworkers,...
            fc,nx,ny,xlen,ylen,ci_ave,sig_dig,f_resolution,picked_steps)
%         dlmwrite([folder 'exporting_figures/FT/FT_processing_status'],1)
    end
    
end

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

if fourier_transform_new_on==1
    
    
    error('Just stop here for now. I need to save xfft and yfft somewhere')
    
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
        export_fourier_transform_data_2d_new(folder,nworkers,nIterations,...
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
    warning('beta being set to high value')
    beta=50;
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

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

function sol=char_f_trend_x_dir(c_relative,nx,ny,ylen)
sol=zeros(1,nx);
for xth=1:1:nx
    data=c_relative(:,xth);
    [fft1d,xfft,~,~]=fftn_new(data,[ylen],[ny]);
    [max_val,index]=max(abs(fft1d));
    
    %-----------------------------
    max_val
    if max_val>10
        sol(xth)=xfft(index);
    end
    %-----------------------------
        
end
sol=abs(sol);
end

function sol=get_phase_specific(real_terms,imaginary_terms,resolution)

sol=zeros(resolution,resolution);

for i=1:1:resolution
    for j=1:1:resolution
        sol(i,j)=atan_sp(real_terms(i,j),imaginary_terms(i,j));
%         sol(i,j)=real_terms(i,j)^2+imaginary_terms(i,j)^2;
%         sol(i,j)=imaginary_terms(i,j);

    end
end


end

function sol=atan_sp(x,y)

if y>0
    if x>0
        sol=atan(y/x);
    elseif x<0
        sol=pi-atan(y/-x);
    elseif x==0
        sol=pi/2;
    end
elseif y<0
    if x>0
        sol=2*pi-atan(-y/x);
    elseif x<0
        sol=pi+atan(y/x);
    elseif x==0
        sol=3*pi/2;
    end
elseif y==0
    if x>0
        sol=0;
    elseif x<0
        sol=pi;
    elseif x==0
        error('The axis has zero length')
    end
end

end

function inverted=invert_FT_2d(resolution,k1_domain,k2_domain,fourier,fc,nworkers)

i_complex=sqrt(-1);

k_domain_new=linspace(-fc,fc,resolution);
[k1_domain_new,k2_domain_new] = meshgrid(k_domain_new);
fourier_new=interp2(k1_domain,k2_domain,fourier',k1_domain_new,k2_domain_new);
fourier_new=fourier_new';

[allocation,extra,load_less,load_more]=wise_task_splitter_v2(...
    resolution,nworkers);

twopi=2*pi;
temp=zeros(nworkers,load_more,resolution);
parfor worker=1:1:nworkers
    
    if worker<=extra
        temp(worker,:,:)=inverse_FT_assist(allocation(worker,:),load_more,load_more,...
            resolution,k_domain_new,fourier_new,i_complex,twopi);
    elseif worker>extra
        temp(worker,:,:)=inverse_FT_assist(allocation(worker,:),load_less,load_more,...
            resolution,k_domain_new,fourier_new,i_complex,twopi);
    end
end

inverted=zeros(resolution,resolution);
index=0;
for worker=1:1:nworkers
    if worker<=extra
        load=load_more;
    elseif worker>extra
        load=load_less;
    end
    for i=1:1:load
        index=index+1;
        for j=1:1:resolution
            inverted(index,j)=temp(worker,i,j);
        end
    end
end

end

function sol=inverse_FT_assist(keys,load,load_more,...
    resolution,k_domain,fourier,i_complex,twopi)
sol=zeros(load_more,resolution);
for i=1:1:load
    for j=1:1:resolution
        sol(i,j)=inverse_FT_2d(fourier,k_domain,keys(i),j,twopi,resolution,i_complex);
    end
end
end

function inverted=inverse_FT_2d(fourier,k_domain,m,n,twopi,resolution,i_complex)
inverted=0;
mM=(m-1)/resolution;
nN=(n-1)/resolution;
for ik1=1:1:resolution
    for ik2=1:1:resolution
        term=twopi*(mM*k_domain(ik1)+nN*k_domain(ik2));
        inverted=inverted+fourier(ik1,ik2)*(cos(term)+i_complex*sin(term));
        
    end
end
inverted=inverted/resolution^2;
end


function [real_terms,imaginary_terms,real_resolution,k1_domain,k2_domain]=....
    load_FT_w_axial_2d(folder,iteration,fc,twopi,nx,ny,ci_ave)

% time=dlmread([t_folder num2str(iteration)]);
[real_terms,imaginary_terms,resolution]=load_FT_2d(folder,iteration);
k1_domain=linspace(-fc,fc,resolution);
k2_domain=linspace(-fc,fc,resolution);
% A=real_terms.^2+imaginary_terms.^2;
if mod(resolution,2)==0
    real_resolution=resolution+1;
    half_resolution=resolution/2;
    iCenter=half_resolution+1;
    real_terms=[real_terms(1:1:half_resolution,1:1:half_resolution) zeros(half_resolution,1) real_terms(1:1:half_resolution,half_resolution+1:1:resolution);
        zeros(1,real_resolution);
        real_terms(half_resolution+1:1:resolution,1:1:half_resolution) zeros(half_resolution,1) real_terms(half_resolution+1:1:resolution,half_resolution+1:1:resolution)];
    imaginary_terms=[imaginary_terms(1:1:half_resolution,1:1:half_resolution) zeros(half_resolution,1) imaginary_terms(1:1:half_resolution,half_resolution+1:1:resolution);
        zeros(1,real_resolution);
        imaginary_terms(half_resolution+1:1:resolution,1:1:half_resolution) zeros(half_resolution,1) imaginary_terms(half_resolution+1:1:resolution,half_resolution+1:1:resolution)];
    k1_domain=[k1_domain(1:1:half_resolution) 0 k1_domain(half_resolution+1:1:resolution)];
    k2_domain=[k2_domain(1:1:half_resolution) 0 k2_domain(half_resolution+1:1:resolution)];
    c=dlmread([folder 'concentration/iteration_' num2str(iteration)]);
    c_nodal=extract_nodal_weights_2D(c,nx,ny);
    c_relative=c_nodal-ci_ave*ones(ny,nx);
    for j=1:1:real_resolution
        [real,imaginary]=FT_real_imaginary_2d(...
            c_relative,k2_domain(j),0,ny,nx,twopi);
        real_terms(iCenter,j)=real;
        imaginary_terms(iCenter,j)=imaginary;
        [real,imaginary]=FT_real_imaginary_2d(...
            c_relative,0,k1_domain(j),ny,nx,twopi);
        real_terms(j,iCenter)=real;
        imaginary_terms(j,iCenter)=imaginary;
    end
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


function export_fourier_transform_data_specific_2d(folder,nworkers,...
    fc,nx,ny,xlen,ylen,ci_ave,sig_dig,f_resolution,picked_steps)

% nIterations=1;

% ONLY GOOD FOR EQUAL MESH

if nx~=ny || xlen~=1 || ylen~=1
    error('This is designed specifically for equal spaced cube of volume = 1 unit squared')
end

[allocation,extra,load_less,load_more]=wise_task_splitter_v2(...
    6,nworkers);

% temp=zeros(nworkers,load_more,f_resolution,f_resolution,f_resolution,2);


f_x=linspace(-fc,fc,f_resolution);
f_y=linspace(-fc,fc,f_resolution);
% f_z=linspace(-fc,fc,f_resolution);

twopi=2*pi;

parfor worker=1:1:nworkers %This is parfor
    if worker<=extra
        FT_data_assist_specific_2d(...
            allocation(worker,:),load_more,load_more,nx,ny,f_x,f_y,...
            folder,ci_ave,twopi,f_resolution,sig_dig,picked_steps);
    elseif worker>extra
        FT_data_assist_specific_2d(...
            allocation(worker,:),load_less,load_more,nx,ny,f_x,f_y,...
            folder,ci_ave,twopi,f_resolution,sig_dig,picked_steps);
    end

end

end

function sol=FT_data_assist_specific_2d(...
    keys,load,load_more,nx,ny,...
    f_x,f_y,folder,ci_ave,twopi,...
    f_resolution,sig_dig,picked_steps)

% THIS ONLY WRITES THE FILES NOW, SO IT IS NOT NECESSARY TO OUTPUT SOL
% FIX LATER

% [real,imaginary]=FT_real_imaginary_3d(c,k1,k2,k3,M,N,O,twopi)
sol=zeros(load_more,f_resolution,f_resolution,2);
real_terms=zeros(f_resolution,f_resolution);
imaginary_terms=zeros(f_resolution,f_resolution);

for i=1:1:load
    
    if exist([folder 'exporting_figures/FT/real_term/iteration_' num2str(picked_steps(keys(i)))],'file')
        return
    end

    c=dlmread([folder 'concentration/iteration_' num2str(picked_steps(keys(i)))]);
    c=extract_nodal_weights_2D(c,nx,ny);
    c=c-ci_ave*ones(ny,nx);
    
    %=======================================
    for xth=1:1:f_resolution
        for yth=1:1:f_resolution
                
            [real,imaginary]=FT_real_imaginary_2d(...
                c,f_y(yth),f_x(xth),ny,nx,twopi);

            sol(i,yth,xth,1)=real;
            sol(i,yth,xth,2)=imaginary;
        end
    end
    

    
%     [data_fft,xfft,yfft,zfft]=fftn_new(data,lens,mesh)
% 
%     for xth=1:1:f_resolution
%         for yth=1:1:f_resolution
%                 
%             [real,imaginary]=FT_real_imaginary_2d(...
%                 c,f_y(yth),f_x(xth),ny,nx,twopi);
% 
%             sol(i,yth,xth,1)=real;
%             sol(i,yth,xth,2)=imaginary;
%         end
%     end
%     
    
    %==================================================
    
    %------------EXPORTING HERE---------------------------
    for xth=1:1:f_resolution
        for yth=1:1:f_resolution
            real_terms(xth,yth)=sol(i,yth,xth,1);
            imaginary_terms(xth,yth)=sol(i,yth,xth,2);
        end
    end
    
    dlmwrite([folder 'exporting_figures/FT/real_term/iteration_' num2str(picked_steps(keys(i)))],real_terms,'precision',sig_dig)
    dlmwrite([folder 'exporting_figures/FT/imaginary_term/iteration_' num2str(picked_steps(keys(i)))],imaginary_terms,'precision',sig_dig)
    %------------------------------------------------------
    
end

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

function export_fourier_transform_data_2d_new(folder,nworkers,nIterations,...
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
        FT_data_assist_2d_new(...
            allocation(worker,:),load_more,load_more,nx,ny,f_x,f_y,...
            folder,ci_ave,twopi,f_resolution,sig_dig,xlen,ylen);
    elseif worker>extra
        FT_data_assist_2d_new(...
            allocation(worker,:),load_less,load_more,nx,ny,f_x,f_y,...
            folder,ci_ave,twopi,f_resolution,sig_dig,xlen,ylen);
    end

end

end

function sol=FT_data_assist_2d_new(...
    keys,load,load_more,nx,ny,...
    f_x,f_y,folder,ci_ave,twopi,...
    f_resolution,sig_dig,xlen,ylen)

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
    
    %=======================================
    
    [data_fft,xfft,yfft,~]=fftn_new(c,[ylen xlen],[ny,nx]);

%     for xth=1:1:f_resolution
%         for yth=1:1:f_resolution
%                 
%             [real,imaginary]=FT_real_imaginary_2d(...
%                 c,f_y(yth),f_x(xth),ny,nx,twopi);
% 
%             sol(i,yth,xth,1)=real;
%             sol(i,yth,xth,2)=imaginary;
%         end
%     end
    
    
    %==================================================
    
    %------------EXPORTING HERE---------------------------
    for xth=1:1:length(xfft)
        for yth=1:1:length(yfft)
            real_terms(xth,yth)=real(data_fft(yth,xth));
            imaginary_terms(xth,yth)=imag(data_fft(yth,xth));
        end
    end
    
    dlmwrite([folder 'exporting_figures/FT/real_term/iteration_' num2str(keys(i))],real_terms,'precision',sig_dig)
    dlmwrite([folder 'exporting_figures/FT/imaginary_term/iteration_' num2str(keys(i))],imaginary_terms,'precision',sig_dig)
    %------------------------------------------------------
    
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
    
    %=======================================
    for xth=1:1:f_resolution
        for yth=1:1:f_resolution
                
            [real,imaginary]=FT_real_imaginary_2d(...
                c,f_y(yth),f_x(xth),ny,nx,twopi);

            sol(i,yth,xth,1)=real;
            sol(i,yth,xth,2)=imaginary;
        end
    end
    

    
%     [data_fft,xfft,yfft,zfft]=fftn_new(data,lens,mesh)
% 
%     for xth=1:1:f_resolution
%         for yth=1:1:f_resolution
%                 
%             [real,imaginary]=FT_real_imaginary_2d(...
%                 c,f_y(yth),f_x(xth),ny,nx,twopi);
% 
%             sol(i,yth,xth,1)=real;
%             sol(i,yth,xth,2)=imaginary;
%         end
%     end
%     
    
    %==================================================
    
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

function sol=FT_circular_deliberate_specific_assist(keys,load,load_more,...
    nx,ny,c_folder,t_folder,ci_ave,cos_list,sin_list,twopi,beta)
sol=zeros(load_more,2);
for i=1:1:load
    c=dlmread([c_folder num2str(keys(i))]);
    c_nodal=extract_nodal_weights_2D(c,nx,ny);
    c_relative=c_nodal-ci_ave*ones(ny,nx);
    summation=0;
    for j=1:1:beta
        summation=summation+get_mag_2d(c_relative,...
            cos_list(j),sin_list(j),ny,nx,twopi);
    end
    sol(i,1)=log(summation/beta);
    sol(i,2)=dlmread([t_folder num2str(keys(i))]);
end
end