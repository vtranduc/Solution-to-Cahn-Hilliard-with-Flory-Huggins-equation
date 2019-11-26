function current_nworkers=cpu_initializer(nworkers)
myCluster=parcluster('local');
if nworkers>myCluster.NumWorkers
    warning('nworkers specified exceeds number of workers available.\nThe number of workers will be reduced to maximum available, which is %d.',myCluster.NumWorkers)
    nworkers=myCluster.NumWorkers;
elseif nworkers==0
    nworkers=myCluster.NumWorkers;
elseif nworkers==-1
    prompt=['Please specify number of workers to be used. The maximum available is ' num2str(myCluster.NumWorkers) '\n'];
    while 1
        nworkers=input(prompt);
        if length(nworkers)==1 && mod(nworkers,1)==0 && nworkers<=myCluster.NumWorkers && nworkers>=1
            break
        else
            warning('The input is invalid. Please specify a number equal to or lower than %d.',myCluster.NumWorkers)
        end
    end
end
current=gcp('nocreate');
if isempty(current)
    current_nworkers=1;
else
    current_nworkers=current.NumWorkers;
end

if current_nworkers==nworkers
    return
end
while 1
    try
        delete(gcp('nocreate'))
        parpool('local',nworkers);
        break
    catch
        nworkers=nworkers-1;
        if nworkers==0
            error('No core is available for simulation! Simulation will terminate!')
        end
        warning('Though resources are physically available, they are not usable at the moment. nworkers will be reduced to %d.',nworkers')
    end
end
current=gcp('nocreate');
if isempty(current)
    current_nworkers=1;
else
    current_nworkers=current.NumWorkers;
end
if nworkers==1
    fprintf('%d worker will be used for this run.\n',current_nworkers)
else
    fprintf('%d workers will be used for this run.\n',current_nworkers)
end
end