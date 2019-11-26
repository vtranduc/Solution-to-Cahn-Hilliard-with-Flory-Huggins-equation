function wwTerms=generate_wwTerms_serendipity_3d(weights)
wwTerms=zeros(8,4,8,4,3,3,3);
for orientation1=1:1:8
    for type1=1:1:4
        
        for orientation2=1:1:8
            for type2=1:1:4
        
                for iz=1:1:3
                    for iy=1:1:3
                        for ix=1:1:3
                            wwTerms(orientation1,type1,orientation2,type2,ix,iy,iz)=...
                                weights(orientation1,type1,2,ix,iy,iz)...
                                *weights(orientation2,type2,2,ix,iy,iz)...
                                +weights(orientation1,type1,3,ix,iy,iz)...
                                *weights(orientation2,type2,3,ix,iy,iz)...
                                +weights(orientation1,type1,4,ix,iy,iz)...
                                *weights(orientation2,type2,4,ix,iy,iz);
                        end
                    end
                end
                
            end
        end

    end
end
end