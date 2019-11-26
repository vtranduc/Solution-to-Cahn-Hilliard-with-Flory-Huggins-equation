function test24

clear
clc

% tic


% myCluster=parcluster('local')
% myCluster.NumWorkers
% 
% return

nworkers=0;
% nworkers=0 if you wanna use maximum number of cores
% nworkers=-1 if you wanna specify it as user inpput
% nworkers=1 if no parallel computing is desired
% nworkers>1 for parallel computing with specified number of workers

myCluster=parcluster('local');
if nworkers>myCluster.NumWorkers
    warning('nworkers specified exceeds number of workers available.\nThe number of workers will be reduced to maximum available, which is %d.',myCluster.NumWorkers)
elseif nworkers==0
    nworkers=myCluster.NumWorkers;
elseif nworkers==-1
    while 1
        prompt=['Please specify number of workers to be used. The maximum available is ' num2str(myCluster.NumWorkers) '\n'];
        nworkers=input(prompt);
        if length(nworkers)==1 && mod(nworkers,1)==0 && nworkers<=myCluster.NumWorkers
            break
        else
            warning('The input is invalid. Please specify a number equal to or lower than %d.',myCluster.NumWorkers)
        end
    end
end
if nworkers~=1
    while 1
        try
            delete(gcp('nocreate'))
            parpool('local',nworkers)
            break
        catch
            nworkers=nworkers-1;
            warning('Though resources are physically available, they are not usable at the moment. nworkers will be reduced to %d.',nworkers')
        end
    end
end
if nworkers==1
    fprintf('%d worker will be used for this run.\n',nworkers)
else
    fprintf('%d workers will be used for this run.\n',nworkers)
end

% return
% 
%     
% ncores=input('How many cores do you wanna use?')
% 
% if ncores>1
%     try
%         delete(gcp('nocreate'))
%         parpool('local',ncores)
%     catch
%         error('Too many!')
%     end
% end
% fprintf('done!\n')
% disp('SUCCESSFUL BABY!')
% 
% toc


end