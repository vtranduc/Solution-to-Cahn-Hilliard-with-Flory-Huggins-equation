function test109

clear
clc
close all

%----------------------------- Specifications
name='5000 0.5 0.3 1';
picked_steps=[100 10 50 80 100 158]; %Make it so only 

profile_6_figs=1;


if profile_6_figs==1
    ImageSizeX_profile_6_figs=6;
    ImageSizeY_profile_6_figs=8;

end

profile_on=0;
isosurface_on=0;

fourier_transform_on=1;

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


% k_target=get_characteristic_frequency(diff,ci_ave,mean(T),n1,n2,entropy);
if ~exist([folder 'exporting_figures'],'dir')
    mkdir([folder 'exporting_figures'])
end

exporting_folder=[folder 'exporting_figures/'];

%--------------------------------
spec=dlmread([folder 'specification/spec_general']);


nex=spec(1);
ney=spec(2);
nez=spec(3);

% dt
% obs_t
% diff
% include_fluc
ci_fluc=spec(8);
% tol_nr
% tol_jk
% time_out
% max_frame
% entropy
% T_theta
% dt_down
% dt_up
% dt_min
% dt_ideal
% bypass_tol
% n1
% n2

% --- Assumption ---

xlen=1;
ylen=1;
zlen=1;

%---------------------------------

nx=nex+1;
ny=ney+1;
nz=nez+1;

x=linspace(0,xlen,nx);
y=linspace(0,ylen,ny);
z=linspace(0,zlen,nz);

ci_ave=mean(ci);




% time_domain=[];
% mag_sq_log=[];

if profile_on==1
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
            title(sprintf('t=%d',time),'FontWeight','Normal')
            hold on
            plot3([0 0],[0 0],[0 1],'k')
            plot3([0 1],[0 0],[1 1],'k')
            plot3([0 0],[0 1],[1 1],'k')
                hold off
            set(gca,'FontName',axis_font,'FontSize',font_size)
            set(gcf,'PaperUnits','inches','PaperPosition',[0 0 ImageSizeX_profile_6_figs ImageSizeY_profile_6_figs])
            print('-dpng', strcat(exporting_folder,'isosurface.png') , strcat('-r',num2str(dpi)))
        end
        close(my_fig)
    end
end

if fourier_transform_on==1
    
    f_resolution=10;
    sig_dig=17;
    
    if ~exist([folder 'exporting_figures/FT'],'dir')
        mkdir([folder 'exporting_figures/FT'])
        mkdir([folder 'exporting_figures/FT/real_term'])
        mkdir([folder 'exporting_figures/FT/imaginary_term'])
        dlmwrite([folder 'exporting_figures/FT/FT_processing_status'],0)
        
    end
    if dlmread([folder 'exporting_figures/FT/FT_processing_status'])==0
        
        nIterations=determine_nIterations(folder);
        export_fourier_transform_data(folder,nworkers,nIterations,...
            fc,nx,ny,nz,xlen,ylen,zlen,ci_ave,sig_dig,f_resolution)

        dlmwrite([folder 'exporting_figures/FT/FT_processing_status'],1)
        
    end
    
    
    
end



end

%==========================================================================
%==========================================================================
%==========================================================================
%==========================================================================
%==========================================================================
%==========================================================================
%==========================================================================

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
fprintf('Parallel computing fourier transform now\n')
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

% simul=dlmread([folder 'exporting_figures/FT/real_term/iteration_' num2str(1)]);
% 
% [n_resolution,~]=size(simul);
% 
% % simul
% % 
% % return
% 
% real_terms=zeros(n_resolution,n_resolution,n_resolution);
% imaginary_terms=zeros(n_resolution,n_resolution,n_resolution);
% 
% 
% for zth=1:1:n_resolution
%     index_y=(zth-1)*f_resolution;
%     for xth=1:1:n_resolution
%         for yth=1:1:n_resolution
%             real_terms(xth,yth,zth)=simul(xth,index_y+yth);
%             imaginary_terms(xth,yth,zth)=simul(xth,index_y+yth);
%         end
%     end
% end
% 
% A=sqrt(real_terms.^2+imaginary_terms.^2);
% 
% A
% 
% slice(f_x,f_y,f_z,A,0,0,0)
% colorbar('gray')

end

function sol=FT_data_assist(keys,load,load_more,nx,ny,nz,f_x,f_y,f_z,folder,ci_ave,twopi,f_resolution,sig_dig)

% [real,imaginary]=FT_real_imaginary_3d(c,k1,k2,k3,M,N,O,twopi)
sol=zeros(load_more,f_resolution,f_resolution,f_resolution,2);
real_terms=zeros(f_resolution,f_resolution^2);
imaginary_terms=zeros(f_resolution,f_resolution^2);

for i=1:1:load
    
%     c=dlmread([folder 'concentration/iteration_' num2str(keys(i))]);
    
    c=dlmread([folder 'concentration/iteration_' num2str(100)]);
    
    c=extract_abs_result_serendipity(c,nx,ny,nz);
    c=c-ci_ave*ones(ny,nx,nz);
    
    for xth=1:1:f_resolution
        for yth=1:1:f_resolution
            for zth=1:1:f_resolution
                
                [real,imaginary]=FT_real_imaginary_3d(...
                    c,f_x(xth),f_y(yth),f_z(zth),ny,nx,nz,twopi);
                
                sol(i,xth,yth,zth,1)=real;
                sol(i,xth,yth,zth,2)=imaginary;
                
%                 [i yth xth zth]
                
            end
        end
    end
    
    %------------EXPORTING HERE---------------------------
%     index_y=0;
    for zth=1:1:f_resolution
        index_y=(zth-1)*f_resolution;
        for xth=1:1:f_resolution
            for yth=1:1:f_resolution
                real_terms(xth,index_y+yth)=sol(i,xth,yth,zth,1);
                imaginary_terms(xth,index_y+yth)=sol(i,xth,yth,zth,2);
            end
        end
    end
    
    dlmwrite([folder 'exporting_figures/FT/real_term/iteration_' num2str(keys(i))],real_terms,'precision',sig_dig)
    dlmwrite([folder 'exporting_figures/FT/imaginary_term/iteration_' num2str(keys(i))],imaginary_terms,'precision',sig_dig)
    %------------------------------------------------------
    
end

end

function n=determine_nIterations(folder)
n=0;
while 1==1
    n=n+1;
    if ~exist([folder 'concentration/iteration_' num2str(n)],'file')
        break
    end
end
n=n-1;
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