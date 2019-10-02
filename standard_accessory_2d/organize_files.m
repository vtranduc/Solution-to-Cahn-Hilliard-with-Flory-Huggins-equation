function organize_files(simulation_folder,temporary_path,nframes_last_sim,...
    nframes_new,view_simulation_only)
if view_simulation_only==0
    spec_destination=[simulation_folder 'specification/'];
    spec_new=[temporary_path 'specification/spec_general'];
    if exist(spec_new,'file')
        delete([simulation_folder 'specification/spec_general'])
        movefile(spec_new,spec_destination)
    end
    conc_temp=[temporary_path 'concentration/iteration_'];
    time_temp=[temporary_path 'time/iteration_'];
    conc_destination=[simulation_folder 'concentration/'];
    time_destination=[simulation_folder 'time/'];
    for frame=nframes_last_sim+1:1:nframes_new
        movefile([conc_temp num2str(frame)],conc_destination)
        movefile([time_temp num2str(frame)],time_destination)
    end
    nRuns=2;
    elapsed_file=[simulation_folder 'specification/elapsed_'];
    while exist([elapsed_file num2str(nRuns)],'file')
        nRuns=nRuns+1;
    end
    elapsed_new=[elapsed_file num2str(nRuns)];
    movefile([temporary_path 'specification/elapsed'],elapsed_new)
end
delete([simulation_folder '*.avi'])
movefile([temporary_path '*.avi'],simulation_folder)
rmdir(temporary_path,'s')
end