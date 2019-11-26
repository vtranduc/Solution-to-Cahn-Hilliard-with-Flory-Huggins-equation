function wTerms=generate_wTerms_curved_3d(ne,weights)
wTerms=zeros(ne,8,8,3,3,3);
for e=1:1:ne
    for orientation=1:1:8
        for type=1:1:8
            for iz=1:1:3
                for iy=1:1:3
                    for ix=1:1:3
                        wTerms(e,orientation,type,ix,iy,iz)=...
                            weights(e,orientation,type,5,ix,iy,iz)+...
                            weights(e,orientation,type,6,ix,iy,iz)+...
                            weights(e,orientation,type,7,ix,iy,iz);
                    end
                end
            end
        end
    end
end
end