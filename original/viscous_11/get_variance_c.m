function variance=get_variance_c(...
    ci_ave,ne,ney,nx,ny,dx,dy,weights,c_folder,...
    allocation,extra,load_less,load_more,nworkers,...
    nIterations)

variance=zeros(1,nIterations);
temp=zeros(nworkers,load_more);
parfor worker=1:1:nworkers
    if worker<=extra
        temp(worker,:)=compute_variance_assist(allocation(worker,:),...
            load_more,load_more,c_folder,...
            ci_ave,ne,ney,nx,ny,dx,dy,weights);
    elseif worker>extra
        temp(worker,:)=compute_variance_assist(allocation(worker,:),...
            load_less,load_more,c_folder,...
            ci_ave,ne,ney,nx,ny,dx,dy,weights); 
    end
end
index=0;
for worker=1:1:nworkers
    if worker<=extra
        load=load_more;
    elseif worker>extra
        load=load_less;
    end
    for i=1:1:load
        index=index+1;
        variance(index)=temp(worker,i);
    end
end
end

function solution=compute_variance_assist(...
    keys,load,load_more,c_folder,...
    ci_ave,ne,ney,nx,ny,dx,dy,weights)
solution=zeros(1,load_more);
for i=1:1:load
    c=dlmread([c_folder num2str(keys(i))]);
    solution(i)=iteration_variance_c(c,ci_ave,ne,ney,nx,ny,dx,dy,weights);
end
end

function variance=iteration_variance_c(c,ci_ave,ne,ney,nx,ny,dx,dy,weights)
summation=0;
w=[5/18 4/9 5/18];
variance_local=elemental_variance(ne,nx,ny,ney,c,ci_ave,weights);
for e=1:1:ne
    for ix=1:1:3
        for iy=1:1:3
            summation=summation+w(ix)*w(iy)...
                *variance_local(e,ix,iy);  
        end
    end
end
variance=summation*dx*dy;
end

function solution=elemental_variance(ne,nx,ny,ney,c,ci_ave,weights)
solution=zeros(ne,3,3);
c_=zeros(1,16);
e=0;
for ix=1:1:nx-1
    for iy=1:1:ney
        gbfs=elemental_gbf(ix,iy,ny);
        for i=1:1:16
            c_(i)=c(gbfs(i));
        end
        e=e+1;
        solution(e,:,:)=(conc(c_,weights,1)-ci_ave*ones(3,3)).^2;
    end
end
end