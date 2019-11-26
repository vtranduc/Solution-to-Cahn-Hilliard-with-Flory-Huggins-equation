function conc_=compute_conc_serendipity(ne,c,weights,nex,nexney,nxny,nx)

conc_=zeros(ne,7,3,3,3);

parfor e=1:1:ne
    conc_(e,:,:,:,:)=compute_elemental_elemental_conc_serendipity(e,c,weights,nex,nexney,nxny,nx);
end

end

function elemental_conc_=compute_elemental_elemental_conc_serendipity(e,c,weights,nex,nexney,nxny,nx)

elemental_conc_=zeros(7,3,3,3);

gbfs=get_gbfs_of_element_serendipity(e,nex,nexney,nxny,nx);

for order=1:1:7
    for ix=1:1:3
        for iy=1:1:3
            for iz=1:1:3
                index=0;
                for orientation=1:1:8
                    for type=1:1:4
                        index=index+1;
                        elemental_conc_(order,ix,iy,iz)=elemental_conc_(order,ix,iy,iz)...
                            +c(gbfs(index))...
                            *weights(e,orientation,type,order,ix,iy,iz);
                    end
                end
            end
        end
    end
end

end