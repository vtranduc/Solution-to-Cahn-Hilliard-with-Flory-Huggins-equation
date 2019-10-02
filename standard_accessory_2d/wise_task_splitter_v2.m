function [allocation,extra,load_less,load_more]=wise_task_splitter_v2(nTasks,nworkers)

extra=mod(nTasks,nworkers);
load_less=(nTasks-extra)/nworkers;
if extra==0
    load_more=load_less;
else
    load_more=load_less+1;
end

extra=mod(nTasks,nworkers);

if extra==0
    load=nTasks/nworkers;
    allocation=zeros(nworkers,load);
    index=0;
    for worker=1:1:nworkers
        for i=1:1:load
            index=index+1;
            allocation(worker,i)=index;
        end
    end
else
    load=(nTasks-extra)/nworkers+1;
    allocation=zeros(nworkers,load);
    index=0;
    for worker=1:1:extra
        for i=1:1:load
            index=index+1;
            allocation(worker,i)=index;
        end
    end
    for worker=extra+1:1:nworkers
        for i=1:1:load-1
            index=index+1;
            allocation(worker,i)=index;
        end
    end
end

end