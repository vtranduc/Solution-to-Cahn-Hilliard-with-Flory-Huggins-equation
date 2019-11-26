function main_serendipity_3d_2018

clear
clc
close all


nex=31;
ney=31;
nez=31;

show_figure=2;
export_data=1;

diff=5000;
ci=0.5;
% T=[0.38 0.32];
T=[0.3];
n2=1;
ci_fluc=1.0*10^-6;

%-----------------------------
T_type=4;
T_spec=3;

wetting=1;

thermophoresis=-0;
%-----------------------------

n1=1;
tol_jk=5;
entropy=1;

dt=1.0*10^-10;

nworkers=0;
tol_nr=1.0*10^-6;
communication_reliever=1;
max_frame=100000;
time_out=3600*(60/60)*24;
nr_up=2;
dt_down=0.5;
dt_up=1.1;
dt_min=1.0e-20;
obs_t=1;
bypass_tol=5;
dt_ideal=5.0e-10;
include_fluc=1;
T_theta=1;

% Set up

dto=inf;
nworkers=cpu_initializer(nworkers);
counter_init=tic;
nx=nex+1;ny=ney+1;nz=nez+1;n=nx*ny*nz;nfour=4*n;
nexney=nex*ney;nxny=nx*ny;ne=nexney*nez;
xlen=1;ylen=1;zlen=1;
weights=generateWeights_serendipity_3d(xlen,ylen,zlen,nex,ney,nez);
nnz_=nnz_sj_serendipity(nx,ny,nz);
[gbfs1,gbfs2]=sj_mapper_serendipity(nx,ny,nz,nnz_); %This takes long time
wTerms=generate_wTerms_serendipity_3d(weights);
wwTerms=generate_wwTerms_serendipity_3d(weights);
dxyz=xlen*ylen*zlen/ne;
if length(T)==1
    coef_T=diff*T;
    chi=dimensionless_chi(T,entropy);
    two_chi_n1=2*chi/n1;
    grad_T=0;
    sj_assist=NaN;
elseif length(T)==2
    [coef_T,two_chi_n1]=generate_T_dist_3d(T_type,T_spec,...
        T(1),T(2),ne,nex,ney,nez,xlen,ylen,zlen,entropy,n1);
    grad_T=1;
    sj_assist=2*entropy/n1;
end


%-------Illustrate temperature gradient
% size(coef_T)
% size(two_chi_n1)
% test=zeros(ney,nex,nez);
% test2=zeros(1,nez);test2_spec=[round(nex/2),round(ney/2)];
% e=0;
% for zth=1:1:nez
%     for yth=1:1:ney
%         for xth=1:1:nex
%             e=e+1;
% %             test(yth,xth,zth)=two_chi_n1(e,2,2,2);
%             test(yth,xth,zth)=coef_T(e,2,2,2,3);
%             if xth==test2_spec(1) && yth==test2_spec(2)
%                 test2(zth)=test(yth,xth,zth);
%             end
%         end
%     end
% end
% X=linspace(0,xlen,nex);Y=linspace(0,ylen,ney);Z=linspace(0,zlen,nez);
% subplot(1,2,1)
% slice(X,Y,Z,test,0.5,0.5,0.5)
% subplot(1,2,2)
% plot(Z,test2)
% grid on
% grid minor
% 
% % subplot(1,2,2)
% % test3=zeros(ney,nex);
% % zth=7;
% % e=nex*ney*(zth-1);
% % for yth=1:1:ney
% %     for xth=1:1:nex
% %         e=e+1;
% %         test3(yth,xth)=coef_T(e,2,3,1,1);
% % %         test3(yth,xth)=two_chi_n1(e,3,2,1);
% %     end
% % end
% % size(test3)
% % size(coef_T)
% % surf(X,Y,test3)
% 
% return
%-------------------------------------------------



% [iZeros,iOnes]=brute_force_iZeros_iOnes_serendipity_3d(...
%     gbfs1,gbfs2,nnz_,nworkers,nx,ny,nz);

if wetting==1
    
    wetting_type=1;
    wetting_spec=8;
    hs=[-0.00005 0.00005];
    gs=[-0.00005 0.00005];
    
    [g_x_minus,g_y_minus,g_z_minus,g_x_plus,g_y_plus,g_z_plus,...
        h_x_minus,h_y_minus,h_z_minus,h_x_plus,h_y_plus,h_z_plus]=...
        generate_wetting_pattern(wetting_type,wetting_spec,hs,gs,nx,ny,nz,...
        xlen,ylen,zlen);
    
%     h_z_plus
% %     
%     return
    
%     g_x_minus=0*ones(ny,nz);
%     g_y_minus=0*ones(nx,nz);
%     g_z_minus=0*ones(nx,ny);
%     g_x_plus=0*ones(ny,nz);
%     g_y_plus=0*ones(nx,nz);
%     g_z_plus=-0.0005*ones(nx,ny);
% 
%     h_x_minus=0*ones(ny,nz);
%     h_y_minus=0*ones(nx,nz);
%     h_z_minus=0*ones(nx,ny);
%     h_x_plus=0*ones(ny,nz);
%     h_y_plus=0*ones(nx,nz);
%     h_z_plus=0.00005*ones(nx,ny);
    
    %---Good for 1000 0.3 0.3 1---
%     g_x_minus=0*ones(ny,nz);
%     g_y_minus=0*ones(nx,nz);
%     g_z_minus=0*ones(nx,ny);
%     g_x_plus=0*ones(ny,nz);
%     g_y_plus=0*ones(nx,nz);
%     g_z_plus=-0.0005*ones(nx,ny);
% 
%     h_x_minus=0*ones(ny,nz);
%     h_y_minus=0*ones(nx,nz);
%     h_z_minus=0*ones(nx,ny);
%     h_x_plus=0*ones(ny,nz);
%     h_y_plus=0*ones(nx,nz);
%     h_z_plus=0.00005*ones(nx,ny);
end


if wetting==0
    [iZeros,iOnes]=brute_force_iZeros_iOnes_serendipity_3d(...
        gbfs1,gbfs2,nnz_,nworkers,nx,ny,nz);
elseif wetting==1
    [iZeros,iOnes,g_bc,bc_wetting]=brute_force_iZeros_iOnes_serendipity_wetting_3d(...
        gbfs1,gbfs2,nnz_,nworkers,nx,ny,nz,...
        g_x_minus,g_y_minus,g_z_minus,g_x_plus,g_y_plus,g_z_plus,...
        h_x_minus,h_y_minus,h_z_minus,h_x_plus,h_y_plus,h_z_plus);
end


bypass=0;
if show_figure~=0
    x=linspace(0,xlen,nx);
    y=linspace(0,ylen,ny);
    z=linspace(0,zlen,nz);
end
if show_figure==2
    alpha=100;
    beta=8;
    gamma=8;
    fc=20;
    time_domain=[];
    mag_sq_log=[];
end
if length(ci)==1
    ci_ave=ci;
elseif length(ci)==2
    ci_ave=mean(ci);
end
if export_data==1
    sig_dig=17;
end

if export_data==1
    output_path=folder_setUp_3d(export_data);
    sf=[nex;ney;nez;dt;obs_t;diff;include_fluc;ci_fluc;tol_nr;...
        tol_jk;time_out;max_frame;entropy;T_theta;...
        dt_down;dt_up;dt_min;dt_ideal;bypass_tol;n1;n2];
    dlmwrite([output_path 'specification/spec_general'],sf);
    dlmwrite([output_path 'specification/status'],0);
    dlmwrite([output_path 'specification/spec_T'],T);
    dlmwrite([output_path 'specification/spec_ci'],ci);
end

% Generate initial condition

co=generate_co_serendipity(ci,1,ci_fluc,n,nfour,nx,ny,nz);

if wetting==1
    co=bc_for_c_serendipity_wetting_3d(co,bc_wetting);
end

% Initialize first step

time=0;
c=co;
iframe=1;

% return

if show_figure==1
    frame_pos=[0 0];
    frame_size=[720 720];

    face_colors={'yellow','blue','red'};
    % The order is {Mean Lower Higher}
    view_angle=3;
    transparency=0.9;
    edge_color='none';
    ifig_simulation=1;

    domain_EdgeColor='green';
    domain_FaceColor='white';
    domain_face_transparency=0.1;
    domain_edge_transparency=0.8;

    c_obs=[mean(ci) mean(ci)-ci_fluc mean(ci)+ci_fluc];
    
    X=zeros(1,n);
    Y=zeros(1,n);
    Z=zeros(1,n);
    index=0;
    for iz=1:1:nz
        for iy=1:1:ny
            for ix=1:1:nx
                index=index+1;
                X(index)=x(ix);
                Y(index)=y(iy);
                Z(index)=z(iz);
            end
        end
    end
    fig_domain=[min(X) max(X) min(Y) max(Y) min(Z) max(Z)];

    [x_minus_x,x_minus_y,x_minus_z,x_plus_x,x_plus_y,x_plus_z,...
        y_minus_x,y_minus_y,y_minus_z,y_plus_x,y_plus_y,y_plus_z,...
        z_minus_x,z_minus_y,z_minus_z,z_plus_x,z_plus_y,z_plus_z,...
        x_plot,y_plot,z_plot]=...
        setUp_illustration_curved_structure_3d(nx,ny,nz,X,Y,Z);

    fig=figure(ifig_simulation);
    set(fig,'Position',[frame_pos(1) frame_pos(2) frame_size(1) frame_size(2)])
    clf(fig)
    visual(iframe)=illustrate_and_analyze_serendipity(...
        c,c_obs,fig_domain,face_colors,edge_color,...
        nx,ny,nz,view_angle,transparency,time,...
        x_minus_x,x_minus_y,x_minus_z,x_plus_x,x_plus_y,x_plus_z,...
        y_minus_x,y_minus_y,y_minus_z,y_plus_x,y_plus_y,y_plus_z,...
        z_minus_x,z_minus_y,z_minus_z,z_plus_x,z_plus_y,z_plus_z,...
        x_plot,y_plot,z_plot,domain_EdgeColor,domain_FaceColor,...
        domain_face_transparency,domain_edge_transparency);
elseif show_figure==2
    [time_domain,mag_sq_log]=test105(...
        co,nx,ny,nz,x,y,z,xlen,ylen,zlen,ci_ave,ci_fluc,...
        alpha,beta,gamma,fc,nworkers,time_domain,mag_sq_log,time);
    
end

if export_data==1
    dlmwrite([output_path 'time/iteration_1'],0);
    dlmwrite([output_path 'concentration/iteration_1'],c,'precision',sig_dig);
end

% return

% Enter Time step method

%------------------------------------------ CONTINUE ALGORITHM

% co=dlmread([pwd '/Results_3d/5000 0.3 [0.36 0.35] 1 WET(type=4,spec=NaN,hs=-1,gs=0) d2s[3,0]/concentration/iteration_415']);
% c=dlmread([pwd '/Results_3d/5000 0.3 [0.36 0.35] 1 WET(type=4,spec=NaN,hs=-1,gs=0) d2s[3,0]/concentration/iteration_416']);
% 
% t1=dlmread([pwd '/Results_3d/5000 0.3 [0.36 0.35] 1 WET(type=4,spec=NaN,hs=-1,gs=0) d2s[3,0]/time/iteration_415']);
% t2=dlmread([pwd '/Results_3d/5000 0.3 [0.36 0.35] 1 WET(type=4,spec=NaN,hs=-1,gs=0) d2s[3,0]/time/iteration_416']);
% 
% co=co';
% c=c';
% 
% time=t2;
% iframe=416;
% 
% dto=t2-t1;
% dt=dto;

%-------------------------------------------

while 1
    fprintf('Starting new time step!\n')
    err=inf;
    % Guess the solution
    cp=dt*(c-co)/dto+c;
    [cp,range_check]=adjust_c_serendipity_3d(cp,ci_ave,nx,ny,nz,nfour);
    if range_check==0
        bypass=bypass+1;
        if bypass>bypass_tol
            disp('The concentration is too close to the limit. Simulation will terminate')
            break
        else
            dt=dt*dt_down;
            if dt<dt_min
                disp('Time step is too small!')
                break
            else
                continue
            end
        end
    end
    if wetting==0
        cp=bc_for_c_serendipity_3d(cp,iOnes,gbfs1);
    elseif wetting==1
        cp=bc_for_c_serendipity_wetting_3d(cp,bc_wetting);
    end
    % Preparing for NR
    coo=co;co=c;c=cp;
    jk=0;
    
    % === STARTING NEWTON RAPHSON ===
    while err>tol_nr
        jk=jk+1;
        if jk>tol_jk
            disp('Too many Newton Raphson iterations')
            jk=-1;
            break
        end
        if wetting==1
            fprintf('Force bc within NR...\n')
            c=bc_for_c_serendipity_wetting_3d(c,bc_wetting);
        end
        fprintf('Computing concentrations...\n')
        sj=compute_conc_serendipity_3d(ne,c,weights,nex,nexney,nxny,nx);
        conco_=compute_conco_serendipity_3d(ne,co,weights,nexney,nex,nx,nxny);
        fprintf('Computing terms...\n')
        [sf_terms,sj_terms]=compute_terms_serendipity_3d(...
            ne,sj,two_chi_n1,n1,n2,grad_T,coef_T,...
            conco_,dt,...
            diff,thermophoresis,sj_assist);
        fprintf('Computing residual...\n')
        sf=sf_serendipity_3d(...
            nfour,nx,ny,nz,nex,ney,weights,dxyz,...
            wTerms,sf_terms);
        fprintf('Computing Jacobian...\n')
        sj=sj_serendipity_3d(...
            gbfs1,gbfs2,nnz_,dt,weights,nex,ney,nx,ny,nz,...
            dxyz,wTerms,wwTerms,sj_terms,...
            nworkers,communication_reliever);
        fprintf('Applying bc...\n')
        [sf,sj]=bc_serendipity_3d(sf,sj,iZeros,iOnes,gbfs1);
        if wetting==1
            sj=bc_Jacobian_serendipity_wetting_3d(sj,g_bc);
        end
        
        
        
        fprintf('Sparsing matrix...\n')
        sj=sparse(gbfs1,gbfs2,sj);
        fprintf('Carrying out matrix division...\n')
        sj=sj\-sf';
        fprintf('Evaluating error...\n')
        err=sqrt(sum(sj.^2))
        fprintf('Updating iteration...\n')
        c=c+sj';
        toc(counter_init)
        if isnan(err)
            disp('Likely matrix is too close to or equal to singular.')
%             jk=-2;
            jk=-1;
            break
        elseif ~isreal(sj) || (jk>1 && err>=erro)
            disp('Newton Raphson is diverging or giving imaginary!')
            jk=-1;
            break
        else
            erro=err;
        end
    end
    % === ENDING NEWTON RAPHSON ===
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
    [c,range_check]=adjust_c_serendipity_3d(c,ci_ave,nx,ny,nz,nfour);
    if range_check==0
        disp('Concentration is too close to the limit. Further simulation will be problematic!')
        break
    end
    %At this point, the result for current time step is ready
    time=time+dt;
    iframe=iframe+1;
    if export_data==1
        dlmwrite([output_path 'time/iteration_' num2str(iframe)],time,'precision',sig_dig);
        dlmwrite([output_path 'concentration/iteration_' num2str(iframe)],c,'precision',sig_dig);
    end
    if show_figure==1
        clf(fig)
        visual(iframe)=illustrate_and_analyze_serendipity(...
            c,c_obs,fig_domain,face_colors,edge_color,...
            nx,ny,nz,view_angle,transparency,time,...
            x_minus_x,x_minus_y,x_minus_z,x_plus_x,x_plus_y,x_plus_z,...
            y_minus_x,y_minus_y,y_minus_z,y_plus_x,y_plus_y,y_plus_z,...
            z_minus_x,z_minus_y,z_minus_z,z_plus_x,z_plus_y,z_plus_z,...
            x_plot,y_plot,z_plot,domain_EdgeColor,domain_FaceColor,...
            domain_face_transparency,domain_edge_transparency);
    elseif show_figure==2
        [time_domain,mag_sq_log]=test105(...
            c,nx,ny,nz,x,y,z,xlen,ylen,zlen,ci_ave,ci_fluc,...
            alpha,beta,gamma,fc,nworkers,time_domain,mag_sq_log,time);
        
    end
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
%======CONCLUDE==================

if export_data==1
    steps=[dt;dto];
    dlmwrite([output_path 'specification/time_steps'],steps,'precision',sig_dig);
    dlmwrite([output_path 'specification/elapsed'],toc(counter_init),'precision',sig_dig);
%     dlmwrite([output_path 'specification/termination'],termination_type);
    dlmwrite([output_path 'specification/status'],1);
%     dlmwrite([output_path 'specification/shape'],shape);
end


end