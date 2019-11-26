function wTerms=generate_wTerms_3d(weights)
wTerms=zeros(8,8,3,3,3);
for orientation=1:1:8
    for type=1:1:8
        for iz=1:1:3
            for iy=1:1:3
                for ix=1:1:3
                    wTerms(orientation,type,ix,iy,iz)=...
                        weights(orientation,type,5,ix,iy,iz)+...
                        weights(orientation,type,6,ix,iy,iz)+...
                        weights(orientation,type,7,ix,iy,iz);
                end
            end
        end
    end
end
end