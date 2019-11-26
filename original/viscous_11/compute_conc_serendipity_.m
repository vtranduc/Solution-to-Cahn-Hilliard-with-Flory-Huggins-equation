function conc_=compute_conc_serendipity_(nodal_vals,weights)

conc_=zeros(7,3,3,3);

for order=1:1:7
    for ix=1:1:3
        for iy=1:1:3
            for iz=1:1:3
                index=0;
                for orientation=1:1:8
                    for type=1:1:4
                        index=index+1;
                        conc_(order,ix,iy,iz)=conc_(order,ix,iy,iz)...
                            +nodal_vals(index)...
                            *weights(orientation,type,order,ix,iy,iz);
                    end
                end
            end
        end
    end
end

end