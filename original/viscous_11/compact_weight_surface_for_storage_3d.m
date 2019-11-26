function sol=compact_weight_surface_for_storage_3d(weights_surface,parallel_computing)
ne=size(weights_surface);
ne=ne(1);
sol=zeros(ne,4032);
if parallel_computing==1
    parfor e=1:1:ne
        index=0;
        eInfo=zeros(1,4032);
        for iorientation=1:1:8
            for itype=1:1:8
                for iorder=1:1:7
                    for idim1=1:1:3
                        for idim2=1:1:3
                            index=index+1;
                            eInfo(index)=weights_surface(e,iorientation,itype,iorder,idim1,idim2);
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
        eInfo=zeros(1,4032);
        for iorientation=1:1:8
            for itype=1:1:8
                for iorder=1:1:7
                    for idim1=1:1:3
                        for idim2=1:1:3
                            index=index+1;
                            eInfo(index)=weights_surface(e,iorientation,itype,iorder,idim1,idim2);
                        end
                    end
                end
            end
        end
        sol(e,:)=eInfo;
    end
end
end