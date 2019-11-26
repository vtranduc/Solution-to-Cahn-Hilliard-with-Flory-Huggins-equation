function conc_=compute_conco_serendipity_3d(ne,c,weights,nexney,nex,nx,nxny)
conc_=zeros(ne,3,3,3);
parfor e=1:1:ne %Apparently, parfor makes it worse
    conc_(e,:,:,:)=conc_compute_elemental_conc_order_zero_serendipity(e,c,weights,nexney,nex,nx,nxny);
end
end

function elemental_conc_=conc_compute_elemental_conc_order_zero_serendipity(e,c,weights,nexney,nex,nx,nxny)
elemental_conc_=zeros(3,3,3);
gbfs=get_gbfs_of_element_serendipity(e,nex,nexney,nxny,nx);
for ix=1:1:3
    for iy=1:1:3
        for iz=1:1:3
            index=0;
            for orientation=1:1:8
                for type=1:1:4
                    index=index+1;
                    elemental_conc_(ix,iy,iz)=elemental_conc_(ix,iy,iz)...
                        +c(gbfs(index))*weights(orientation,type,1,ix,iy,iz);
                end
            end
        end
    end
end
end