function solution=compute_elemental_weights_3d(e,c,weights,nexney,nex,nx,ny,gbfs)
solution=zeros(7,3,3,3);
if nargin==7
    gbfs=get_gbfs_of_element(e,nexney,nex,nx,ny);
end
for iorder=1:1:7
    igbf=0;
    for iorientation=1:1:8
        for itype=1:1:8
            igbf=igbf+1;
            for igpz=1:1:3
                for igpy=1:1:3
                    for igpx=1:1:3
                        solution(iorder,igpx,igpy,igpz)=...
                            solution(iorder,igpx,igpy,igpz)...
                            +weights(iorientation,itype,iorder,igpx,igpy,igpz)...
                            *c(gbfs(igbf));
                    end
                end
            end
        end
    end
end
end