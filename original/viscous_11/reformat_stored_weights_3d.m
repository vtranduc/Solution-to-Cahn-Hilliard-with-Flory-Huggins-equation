function sol=reformat_stored_weights_3d(weights,parallel_computing)
ne=size(weights);
ne=ne(1);
sol=zeros(ne,8,8,7,3,3,3);
if parallel_computing==1
    parfor e=1:1:ne
        index=0;
        eInfo=weights(e,:);
        eWeights=zeros(8,8,7,3,3,3);
        for iorientation=1:1:8
            for itype=1:1:8
                for iorder=1:1:7
                    for iz=1:1:3
                        for iy=1:1:3
                            for ix=1:1:3
                                index=index+1;
                                eWeights(iorientation,itype,iorder,ix,iy,iz)=eInfo(index);
                            end
                        end
                    end
                end
            end
        end
        sol(e,:,:,:,:,:,:)=eWeights;
    end
elseif parallel_computing==0
    for e=1:1:ne
        index=0;
        eInfo=weights(e,:);
        for iorientation=1:1:8
            for itype=1:1:8
                for iorder=1:1:7
                    for iz=1:1:3
                        for iy=1:1:3
                            for ix=1:1:3
                                index=index+1;
                                sol(e,iorientation,itype,iorder,ix,iy,iz)=eInfo(index);
                            end
                        end
                    end
                end
            end
        end
    end
end
end