function [iZeros,iOnes]=brute_force_iZeros_iOnes_serendipity_3d(...
    gbfs1,gbfs2,nnz_,nworkers,nx,ny,nz)

%This one has very ugly algorithm

allocation=wise_task_splitter(nnz_,nworkers);

[~,n]=size(allocation);
extra=mod(nnz_,nworkers);
load_less=(nnz_-extra)/nworkers;
if extra==0
    load_more=load_less;
else
    load_more=load_less+1;
end
temp=zeros(nworkers,n);
parfor worker=1:1:nworkers
    if worker<=extra
        load=load_more;
    else
        load=load_less;
    end
    temp(worker,:)=get_keyVals_brute_force(n,allocation(worker,:),load,gbfs1,gbfs2,nx,ny,nz);
end
iZeros=zeros(1,nnz_);
iOnes=zeros(1,4*nx*ny*nz);
i_iZeros=0;
i_iOnes=0;
index=0;
for worker=1:1:nworkers
    if worker<=extra
        load=load_more;
    else
        load=load_less;
    end
    for i=1:1:load
        index=index+1;
        if temp(worker,i)==1
            i_iZeros=i_iZeros+1;
            iZeros(i_iZeros)=index;
        elseif temp(worker,i)==2
            i_iOnes=i_iOnes+1;
            iOnes(i_iOnes)=index;
        end
    end
end
reducer=nnz(iZeros)+1;
iZeros(reducer:1:nnz_)=[];
reducer=nnz(iOnes)+1;
iOnes(reducer:1:4*nx*ny*nz)=[];
% pre_index=-n;
% vals=[];
% for worker=1:1:nworkers
%     if worker<=extra+1
%         pre_index=pre_index+load_more;
%     else
%         pre_index=pre_index+load_less;
%     end
%     nnz_element=find(temp(worker,:));
%     vals=[vals nnz_element+pre_index*ones(1,length(nnz_element))];
% end

end

function sol=get_keyVals_brute_force(n,keys,load,gbfs1,gbfs2,nx,ny,nz)

sol=zeros(1,n);
for i=1:1:load
    [node1,type1]=analyze_gbs_serendipity(gbfs1(keys(i)));
    [xth1,yth1,zth1]=get_xyzth_3d(node1,nx,ny);
    if type1==2 && (xth1==1 || xth1==nx)

        if gbfs1(keys(i))==gbfs2(keys(i))
            sol(i)=2;
        else
            sol(i)=1;
        end

    
    elseif type1==3 && (yth1==1 || yth1==ny)

        if gbfs1(keys(i))==gbfs2(keys(i))
            sol(i)=2;
        else
            sol(i)=1;
        end
            

    elseif type1==4 && (zth1==1 || zth1==nz)

        if gbfs1(keys(i))==gbfs2(keys(i))
            sol(i)=2;
        else
            sol(i)=1;
        end

    end
end
end