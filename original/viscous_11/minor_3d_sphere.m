function minor_3d_sphere

clear
clc

% === Start the timer =====================================================

counter_init=tic;

% === User defined variables ==============================================

% --- Define mesh ---------------------------------------------------------

nex=5;
ney=4;
nez=3;

% --- Polymer properties --------------------------------------------------

n1=1;
n2=10;
 
diff=-5000;

%======================
ci=0.45;
T=[0.62];
%======================

include_fluc=1;
ci_fluc = 0.0001;

entropy=1;
T_theta=1;

ci_critical=0; % This bypasses ci specified above
T_critical=0; % This bypasses T specified above

% --- Simulation properties -----------------------------------------------

parallel_computing=0; % Non parallel computing uses traditional method

show_figure=1;
export_data=0; % Yes=1, No=0


dt=1.0*10^-7;

obs_t=2.0e-1;

nr_tol=1.0e-6;
nr_max_iteration=20;
%======================
time_out=3600*16;
max_frame=100;
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

% === Validate inputs =====================================================

% VALIDATE THEM HERE

% === Set up variables ====================================================

[nx,ny,nz,n,neight]=setUp_3d(nex,ney,nez);

if ci_critical==1 && T_critical==1
    [ci,T]=identify_critical_point(n1,n2,entropy,T_theta);
elseif ci_critical==1
    [ci,~]=identify_critical_point(n1,n2,entropy,T_theta);
elseif T_critical==1
    [~,T]=identify_critical_point(n1,n2,entropy,T_theta);
end

[gbfs1,gbfs2,iZeros,iOnes]=sj_mapper_3d(nx,ny,nz);
nnz_=nnz_sj_3d(nx,ny,nz);


c_obs=[mean(ci)-ci_fluc mean(ci) mean(ci)+ci_fluc];

%----------------------------------------------------
%---------------------------------------------------- Sphere geometry will
%be applied here instead

% dx=1/nex;
% dy=1/ney;
% dz=1/nez;

% weights=weight_adjuster_3d(generateWeights_3d(),dx,dy,dz);


% dxyz=dx*dy*dz;
% 
% x_coord=linspace(0,1,nx);
% y_coord=linspace(0,1,ny);
% z_coord=linspace(0,1,nz);

[X,Y,Z]=generate_sphere_mesh_3d(nx,ny,nz,n);

if show_figure==1
    [x_minus_x,x_minus_y,x_minus_z,x_plus_x,x_plus_y,x_plus_z,...
    y_minus_x,y_minus_y,y_minus_z,y_plus_x,y_plus_y,y_plus_z,...
    z_minus_x,z_minus_y,z_minus_z,z_plus_x,z_plus_y,z_plus_z]=...
    setUp_illustration_sphere_3d(nx,ny,nz,X,Y,Z);
end




%----------------------------------------------------
%----------------------------------------------------

ne=nex*ney*nez;
nexney=nex*ney;



dto=inf;
if predict_with_acceleration==1
    dtoo=inf;
end

bypass=0;

%------ Test -------------

[n1xth,n1yth,n1zth]=get_n1xyzth_3d(e,nex,nexney);




% % weight_transformer_3d(generateWeights_3d(),ne,X,Y,Z)
% 
% e=ne;
% 
% [n1xth,n1yth,n1zth]=get_n1xyzth_3d(e,nex,nexney);
% 
% nodes=zeros(1,8);
% 
% nodes(1)=get_node_3d(n1xth,n1yth,n1zth,nx,ny);
% nodes(2)=get_node_3d(n1xth+1,n1yth,n1zth,nx,ny);
% nodes(3)=get_node_3d(n1xth,n1yth+1,n1zth,nx,ny);
% nodes(4)=get_node_3d(n1xth+1,n1yth+1,n1zth,nx,ny);
% 
% nodes(5)=get_node_3d(n1xth,n1yth,n1zth+1,nx,ny);
% nodes(6)=get_node_3d(n1xth+1,n1yth,n1zth+1,nx,ny);
% nodes(7)=get_node_3d(n1xth,n1yth+1,n1zth+1,nx,ny);
% nodes(8)=get_node_3d(n1xth+1,n1yth+1,n1zth+1,nx,ny);
% 
% plot3(X,Y,Z,'b.')
% 
% hold on
% for i=1:1:8
%     plot3(X(nodes(i)),Y(nodes(i)),Z(nodes(i)),'ro')
% end
% 
% xlabel('x')
% ylabel('y')
% zlabel('z')
% 
% grid on
% 
% isoparametrically_mapped_gps=isoparametric_mapping_gps_3d(ne,nex,nexney,nx,ny,X,Y,Z,parallel_computing);
% % 
% % % isoparametrically_mapped_gps
% % 
% % isoparametrically_mapped_gps(e,:,:,:,:);
% % 
% % e=13
% 
% for i=1:1:3
%     for j=1:1:3
%         for k=1:1:3
%             plot3(isoparametrically_mapped_gps(e,i,j,k,1),isoparametrically_mapped_gps(e,i,j,k,2),isoparametrically_mapped_gps(e,i,j,k,3),'gx')
%         end
%     end
% end
% 
% hold off

return

% ---------------------------------------------------

if parallel_computing==1 %This can be implemented on traditional algorithm as well
    wTerms=generate_wTerms_3d(weights);
end

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

% === Generate initial condition ==========================================

co=generate_co_3d(ci,include_fluc,ci_fluc,n,neight,nx,ny,nz);

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
    visual(iframe)=illustrate_and_analyze_3d_sphere(c,nx,ny,nz,time,...
        x_minus_x,x_minus_y,x_minus_z,x_plus_x,x_plus_y,x_plus_z,...
        y_minus_x,y_minus_y,y_minus_z,y_plus_x,y_plus_y,y_plus_z,...
        z_minus_x,z_minus_y,z_minus_z,z_plus_x,z_plus_y,z_plus_z);
end

% --- Output the first time step --- %

if export_data==1
    dlmwrite([output_path 'time/iteration_1'],0);
    dlmwrite([output_path 'concentration/iteration_1'],c);
end

% ==== Enter main loop ====================================================

end