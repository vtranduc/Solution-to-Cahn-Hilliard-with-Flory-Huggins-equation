function test29

clear
clc

run_number=8

[nex,ney,nez,obs_t,diff,include_fluc,ci_fluc,nr_tol,...
    nr_max_iteration,time_out,max_frame,entropy,T_theta,...
    dt_down,dt_up,dt_min,dt_ideal,bypass_tol,n1,n2,...
    elapsed,ci,T,dts,c,co,coo,time,...
    last_iteration,output_path,cancellation]=...
    import_data_3d(run_number);

time_out/3600

max_frame

obs_t

bypass_tol
dt_down

dt_min

cancellation

end