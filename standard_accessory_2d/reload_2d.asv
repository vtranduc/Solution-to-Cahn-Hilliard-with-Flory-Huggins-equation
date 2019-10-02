function [goal_reached,modification,temp_output_path,nex,ney,diff,...
    include_fluc,ci_fluc,nr_tol,nr_max_iteration,entropy,T_theta,...
    dt_down,dt_up,dt_min,dt_ideal,bypass_tol,n1,n2,dt_updated,...
    sim_data,c,co,iteration_last,T,ci,time_last]=reload_2d(...
    simulation_folder,obs_t,dt,max_frame,...
    view_simulation_only)
goal_reached=0;
modification=0;
if ~exist(simulation_folder,'dir')
    error('The folder does not exist!')
end
if ~exist([simulation_folder 'specification/status'],'file')
    error('This simulation does not have data exported!')
end
if dlmread([simulation_folder 'specification/status'])==0
    error('The last run did not complete!')
end
iteration_conc=[simulation_folder 'concentration/iteration_'];
iteration_last=1;
while exist([iteration_conc num2str(iteration_last)],'file')
    iteration_last=iteration_last+1;
end
iteration_last=iteration_last-1;
if iteration_last==0
    error('No iteration found!')
end
ci=dlmread([simulation_folder 'specification/spec_ci']);
T=dlmread([simulation_folder 'specification/spec_T']);
sim_data=dlmread([simulation_folder 'specification/spec_general']);
nex=round(sim_data(1));
ney=round(sim_data(2));
diff=sim_data(5);
include_fluc=round(sim_data(6));
ci_fluc=sim_data(7);
nr_tol=sim_data(8);
nr_max_iteration=sim_data(9);
entropy=sim_data(12);
T_theta=sim_data(13);
dt_down=sim_data(14);
dt_up=sim_data(15);
dt_min=sim_data(16);
dt_ideal=sim_data(17);
bypass_tol=round(sim_data(18));
n1=sim_data(19);
n2=sim_data(20);
if length(n2)==21
    n2(2)=sim_data(21);
end
time_last=dlmread([simulation_folder 'time/iteration_' num2str(iteration_last)]);

if iteration_last==1
    dt_updated=dt;
else
    dt_updated=time_last-dlmread([simulation_folder 'time/iteration_' num2str(iteration_last-1)]);
end

if time_last>=obs_t
    disp('Last iteration has already reached specified observation time. Please increase observation time for further simulation!')
    goal_reached=1;
end
if sim_data(4)<obs_t
    sim_data(4)=obs_t;
    modification=1;
end
if iteration_last>=max_frame
    disp('Maximum iteration specified has been reached! Please specify higher max_frame!')
    if goal_reached==0
        goal_reached=1;
    end
end
if sim_data(11)<max_frame
    sim_data(11)=max_frame;
    if modification==0
        modification=1;
    end
end
temp_output_path=[simulation_folder 'temporary/'];
if exist(temp_output_path,'dir')
    rmdir(temp_output_path,'s')
end
co=dlmread([simulation_folder 'concentration/iteration_' num2str(iteration_last-1)]);
c=dlmread([simulation_folder 'concentration/iteration_' num2str(iteration_last)]);

if view_simulation_only==0
    mkdir(temp_output_path);
    mkdir([temp_output_path 'specification/'])
    mkdir([temp_output_path 'time/'])
    mkdir([temp_output_path 'concentration/'])
    dlmwrite([temp_output_path 'specification/status'],0);
end
end