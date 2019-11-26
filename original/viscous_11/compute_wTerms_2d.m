function wTerms=compute_wTerms_2d(weights)
wTerms=zeros(3,3,16);
for ix=1:1:3
    for iy=1:1:3
        for ibasis=1:1:16
            wTerms(ix,iy,ibasis)=...
                weights(ix,iy,ibasis,4)...
                +weights(ix,iy,ibasis,5);
        end
    end
end
end