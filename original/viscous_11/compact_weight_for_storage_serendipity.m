function sol=compact_weight_for_storage_serendipity(weights,parallel_computing)
ne=size(weights);
ne=ne(1);
sol=zeros(ne,6048);
if parallel_computing==1
    parfor e=1:1:ne
        index=0;
        eInfo=zeros(1,6048);
        for iorientation=1:1:8
            for itype=1:1:4
                for iorder=1:1:7
                    for iz=1:1:3
                        for iy=1:1:3
                            for ix=1:1:3
                                index=index+1;
                                eInfo(index)=weights(e,iorientation,itype,iorder,ix,iy,iz);
                            end
                        end
                    end
                end
            end
        end
        sol(e,:)=eInfo;
    end
elseif parallel_computing==0
    for e=1:1:ne
        index=0;
        eInfo=zeros(1,6048);
        for iorientation=1:1:8
            for itype=1:1:4
                for iorder=1:1:7
                    for iz=1:1:3
                        for iy=1:1:3
                            for ix=1:1:3
                                index=index+1;
                                eInfo(index)=weights(e,iorientation,itype,iorder,ix,iy,iz);
                            end
                        end
                    end
                end
            end
        end
        sol(e,:)=eInfo;
    end
end
end