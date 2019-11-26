function main_analyze_informal_2d


if nargin==0
    clear
end
clc
close all

if ~exist('name','var')
    name='Run_1';
end
if ~exist('distribution_type','var')
    distribution_type=0;
end
if ~exist('T_distribution_spec','var')
    T_distribution_spec=NaN;
end

picked_steps=[1 111 211 295 1049 2493]; %Make it so only 

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
FT_circular_deliberate=0;
FT_specific=0;
video_create_on=0; % THIS WILL ALSO CREATE LOG PLOT OF K_CIRCULAR, K1, AND K2
zoom_on=0;

%----------------
fft_filter_on=0;

test2_on=1;
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

% energy_requirement=-166.3341;

% energy_homogeneous=compute_homogeneous_energy(...
%     ci_ave,diff,coef_T,n1,n2,Two_chi_n1,xlen,ylen,ne,dxdy,grad_T);
% 
% energy=dlmread([folder 'exporting_figures/data_folder/energy.dat']);
% if energy(length(energy))+energy_homogeneous>energy_requirement
%     error('It will never reach required energy')
% end
% 
% for i=1:1:length(energy)
%     if energy(i)+energy_homogeneous<=energy_requirement
%         break
%     end
% end

nIterations=determine_nIterations(folder)

diff_2nd_quench=10000;
T_2nd_quench=0.35;
n1_2nd_quench=1;
n2_2nd_quench=1;
grad_T_2nd_quench=0;

coef_T_2nd_quench=diff_2nd_quench*T_2nd_quench;
Two_chi_n1_2nd_quench=2*dimensionless_chi(T_2nd_quench,entropy)/n1_2nd_quench;
%     Two_chi_n1=2*chi/n1;

energy_requirement=compute_homogeneous_energy(...
        ci_ave,diff_2nd_quench,coef_T_2nd_quench,...
        n1_2nd_quench,n2_2nd_quench,Two_chi_n1_2nd_quench,...
        xlen,ylen,ne,dxdy,grad_T_2nd_quench);

for iteration=53:1:nIterations
    c=dlmread([c_folder num2str(iteration)]);
%     c_nodal=extract_nodal_weights_2D(c,nx,ny);
    
    energy=evaluate_total_energy_v2(...
        coef_T_2nd_quench,n1_2nd_quench,n2_2nd_quench,...
        ne,ney,Two_chi_n1_2nd_quench,dxdy,...
        ny,c,weights,diff_2nd_quench,...
        grad_T_2nd_quench);
    
    [iteration energy energy_requirement]
    
    if energy<=energy_requirement
        break
    end
    
end

iteration
energy

% c=dlmread([c_folder num2str(keys(i))]);
% sol(i)=evaluate_total_energy_v2(...
%     coef_T,n1,n2,ne,ney,Two_chi_n1,dxdy,...
%     ny,c,weights,diff,...
%     grad_T);

% i
% [energy(i-1)+energy_homogeneous energy(i)+energy_homogeneous]


end