function [nex,ney,nez,obs_t,diff,include_fluc,ci_fluc,nr_tol,...
    nr_max_iteration,time_out,max_frame,entropy,T_theta,...
    dt_down,dt_up,dt_min,dt_ideal,bypass_tol,n1,n2,...
    elapsed,ci,T,dts,c,co,coo,time,...
    last_iteration,output_path,cancellation]=import_data_3d(run_number)

output_path=[pwd '/Results_3d/Run_' num2str(run_number) '/'];

status=[output_path 'specification/status'];
if ~exist(status,'file')
    error('The data for the run number specified was not found!')
elseif dlmread(status)==0
    error('Simulation specified was not completed!')
end

time_path=[output_path 'time/iteration_'];
con_path=[output_path 'concentration/iteration_'];
spec_path=[output_path 'specification/'];

[nex,ney,nez,~,obs_t,diff,include_fluc,ci_fluc,nr_tol,...
    nr_max_iteration,time_out,max_frame,entropy,T_theta,...
    dt_down,dt_up,dt_min,dt_ideal,bypass_tol,n1,n2]=...
    import_spec_general_3d(spec_path);

elapsed=dlmread([spec_path 'elapsed']);
ci=dlmread([spec_path 'spec_ci']);
T=dlmread([spec_path 'spec_T']);
termination_type=dlmread([spec_path 'termination']);
dts=dlmread([spec_path 'time_steps']);
% dt=dts(1);
% dto=dts(2);
% if predict_with_acceleration==1
%     dtoo=dts(3);
% elseif predict_with_acceleration==0
%     dtoo=0; % This output should be ignored if predict_with_acceleration==0
% end

last_iteration=get_next_iteration_3d(time_path)-1;

time=dlmread([time_path num2str(last_iteration)]);

c=dlmread([con_path num2str(last_iteration)]);
try
    co=dlmread([con_path num2str(last_iteration-1)]);
    try
        coo=dlmread([con_path num2str(last_iteration-2)]);
    catch
        warning('3rd last iteration is not found. coo will be assumed to be equal to co')
        coo=co;
    end
catch
    warning('2nd last iteration is not found. co and coo will be assumed to be equal to c')
    co=c;
    coo=co;
end

% --- Resolve how to go about continuing simulation --- %

% termination_type=4;

 dialog_title='New user specification';
cancellation=0;
 
%--------------------------------------------------------
if termination_type==1

    explanation=sprintf(['Last simulation was terminated because maximum number of frames has been reached.\nIn last simulation, ' num2str(max_frame) ' frames was specified.\n\nPlease enter new maximum number of frames (iterations), including previous ones:']);
    while 1
        answer=inputdlg(explanation,dialog_title);
        if isempty(answer)
            cancellation=1;
            return
        end
        answer=answer{1};
        if ~isempty(answer)
            answer=str2double(answer);
            if isempty(answer)
                err_dialog_3d('Entry is not a number!')
                continue
            end
        else
            err_dialog_3d('Nothing was entered')
            continue
        end
        validation=validate_user_input_3d(answer);
        if validation==0
            err_dialog_3d('Please enter only real scalar number except infinity and NaN. Computation within the box is not allowed.')
            continue
        end
        if answer>max_frame
            max_frame=answer;
            break
        else
            err_dialog_3d(['Please enter frame number higher than ' num2str(max_frame)])
        end
    end
%--------------------------------------------------------
elseif termination_type==2
    explanation=sprintf(['Last simulation was terminated because time limit has been reached.\nIn last simulation, ' num2str(time_out/3600) ' h was specified.\n\nPlease enter new time limit in hours, including previously elapsed time:']);
    while 1
        answer=inputdlg(explanation,dialog_title);
        if isempty(answer)
            cancellation=1;
            return
        end
        answer=answer{1};
        if ~isempty(answer)
            answer=str2double(answer);
            if isempty(answer)
                err_dialog_3d('Entry is not a number!')
                continue
            end
        else
            err_dialog_3d('Nothing was entered')
            continue
        end
        validation=validate_user_input_3d(answer);
        if validation==0
            err_dialog_3d('Please enter only real scalar number except infinity and NaN. Computation within the box is not allowed.')
            continue
        end
        new_time_out=answer*3600;
        if new_time_out>time_out
            time_out=new_time_out;
            break
        else
            err_dialog_3d(['Please enter hours longer than ' num2str(time_out/3600)])
        end
    end
%--------------------------------------------------------
elseif termination_type==3
    
    explanation=sprintf(['Last simulation was terminated because observation time goal has been reached.\nIn last simulation, ' num2str(obs_t) ' dimensionless time unit was specified.\n\nPlease enter new time goal to be reached:']);
    while 1
        answer=inputdlg(explanation,dialog_title);
        if isempty(answer)
            cancellation=1;
            return
        end
        answer=answer{1};
        if ~isempty(answer)
            answer=str2double(answer);
            if isempty(answer)
                err_dialog_3d('Entry is not a number!')
                continue
            end
        else
            err_dialog_3d('Nothing was entered')
            continue
        end
        validation=validate_user_input_3d(answer);
        if validation==0
            err_dialog_3d('Please enter only real scalar number except infinity and NaN. Computation within the box is not allowed.')
            continue
        end
        if answer>obs_t
            obs_t=answer;
            break
        else
            err_dialog_3d(['Please enter time longer than ' num2str(obs_t)])
        end
    end
%--------------------------------------------------------
elseif termination_type==4
    
    explanation=sprintf(['Last simulation was terminated because maximum number of bypasses has been reached.\nIn last simulation, ' num2str(bypass_tol) ' maximum bypasses\nand dt_down of ' num2str(dt_down) ' were specified.\n\nPlease enter new maximum number of bypasses:']);
    explanation={explanation,'Please enter new dt_down:'};
    while 1
        
        answer=inputdlg(explanation,dialog_title);
        if isempty(answer)
            cancellation=1;
            return
        end
        
        answer1=answer{1};
        if ~isempty(answer1)
            answer1=str2double(answer1);
            if isempty(answer1)
                err_dialog_3d('Entry is not a number!')
                continue
            end
        else
            err_dialog_3d('Maximum number of bypasses was not entered')
            continue
        end
        validation=validate_user_input_3d(answer1);
        if validation==0
            err_dialog_3d('Please enter only real scalar number except infinity and NaN. Computation within the box is not allowed.')
            continue
        end
        
        answer2=answer{2};
        if ~isempty(answer2)
            answer2=str2double(answer2);
            if isempty(answer2)
                err_dialog_3d('Entry is not a number!')
                continue
            end
        else
            err_dialog_3d('dt_down was not entered')
            continue
        end
        validation=validate_user_input_3d(answer2);
        if validation==0
            err_dialog_3d('Please enter only real scalar number except infinity and NaN. Computation within the box is not allowed.')
            continue
        end
        if answer1>bypass_tol && answer2<dt_down
            bypass_tol=answer1;
            dt_down=answer2;
            break
        elseif answer1>bypass_tol && answer2==dt_down
            bypass_tol=answer1;
            break
        elseif answer1==bypass_tol && answer2<dt_down
            dt_down=answer2;
            break
        else
            err_dialog_3d(['Please enter new bypass_tol higher than ' num2str(bypass_tol) ' and/or dt_down lower than ' num2str(dt_down)])
        end
    end
%--------------------------------------------------------
elseif termination_type==5

    explanation=sprintf(['Last simulation was terminated because minimum time step has been reached.\nIn last simulation, time step of ' num2str(dt_min) ' dimension time unit was specified.\n\nPlease enter new minimum time step:']);
    while 1
        answer=inputdlg(explanation,dialog_title);
        if isempty(answer)
            cancellation=1;
            return
        end
        answer=answer{1};
        if ~isempty(answer)
            answer=str2double(answer);
            if isempty(answer)
                err_dialog_3d('Entry is not a number!')
                continue
            end
        else
            err_dialog_3d('Nothing was entered')
            continue
        end
        validation=validate_user_input_3d(answer);
        if validation==0
            err_dialog_3d('Please enter only real scalar number except infinity and NaN. Computation within the box is not allowed.')
            continue
        end
        if answer<dt_min
            dt_min=answer;
            break
        else
            err_dialog_3d(['Please enter minimum time step lower than ' num2str(dt_min)])
        end
    end
end

%-----------------------------------------------------------------------------
if time_out<=elapsed
    explanation=sprintf(['time_out is lower than than or equal to elapsed time.\n' num2str(elapsed/3600) ' h has elapsed.\n\nPlease enter new time limit in hours, including previously elapsed time:']);
    while 1
        answer=inputdlg(explanation,dialog_title);
        if isempty(answer)
            cancellation=1;
            return
        end
        answer=answer{1};
        if ~isempty(answer)
            answer=str2double(answer);
            if isempty(answer)
                err_dialog_3d('Entry is not a number!')
                continue
            end
        else
            err_dialog_3d('Nothing was entered')
            continue
        end
        validation=validate_user_input_3d(answer);
        if validation==0
            err_dialog_3d('Please enter only real scalar number except infinity and NaN. Computation within the box is not allowed.')
            continue
        end
        new_time_out=answer*3600;
        if new_time_out>elapsed
            time_out=new_time_out;
            break
        else
            err_dialog_3d(['Please enter hours longer than ' num2str(elapsed/3600) ' h'])
        end
    end
end
if max_frame<=last_iteration
    explanation=sprintf(['max_frame is lower than than or equal to current frame.\n' num2str(last_frame) ' frames have been computed.\n\nPlease enter new max_frame, including previously computed frames:']);
    while 1
        answer=inputdlg(explanation,dialog_title);
        if isempty(answer)
            cancellation=1;
            return
        end
        answer=answer{1};
        if ~isempty(answer)
            answer=str2double(answer);
            if isempty(answer)
                err_dialog_3d('Entry is not a number!')
                continue
            end
        else
            err_dialog_3d('Nothing was entered')
            continue
        end
        validation=validate_user_input_3d(answer);
        if validation==0
            err_dialog_3d('Please enter only real scalar number except infinity and NaN. Computation within the box is not allowed.')
            continue
        end
        new_max_frame=answer;
        if new_max_frame>last_iteration
            max_frame=new_max_frame;
            break
        else
            err_dialog_3d(['Please enter maximum number of frames larger than ' num2str(last_iteration)])
        end
    end
end

end

function next_iteration=get_next_iteration_3d(time_path)

next_iteration=1;
while exist([time_path num2str(next_iteration)],'file')
    next_iteration=next_iteration+1;
end
end

function [nex,ney,nez,dt,obs_t,diff,include_fluc,ci_fluc,nr_tol,...
            nr_max_iteration,time_out,max_frame,entropy,T_theta,...
            dt_down,dt_up,dt_min,dt_ideal,bypass_tol,n1,n2]=...
            import_spec_general_3d(spec_path)

spec=dlmread([spec_path 'spec_general']);

nex=spec(1);
ney=spec(2);
nez=spec(3);
dt=spec(4);
obs_t=spec(5);
diff=spec(6);
include_fluc=spec(7);
ci_fluc=spec(8);
nr_tol=spec(9);
nr_max_iteration=spec(10);
time_out=spec(11);
max_frame=spec(12);
entropy=spec(13);
T_theta=spec(14);
dt_down=spec(15);
dt_up=spec(16);
dt_min=spec(17);
dt_ideal=spec(18);
bypass_tol=spec(19);
n1=spec(20);
n2=spec(21);

end