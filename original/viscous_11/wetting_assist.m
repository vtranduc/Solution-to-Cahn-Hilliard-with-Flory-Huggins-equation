function sol=wetting_assist(weights)
sol=zeros(3,3,16,16);
for ilocal=1:1:16
    for jlocal=1:1:16
        for ix=1:1:3
            for iy=1:1:3
                sol(ix,iy,ilocal,jlocal)=...
                    weights(ix,iy,ilocal,2)*weights(ix,iy,jlocal,2)...
                    +weights(ix,iy,ilocal,3)*weights(ix,iy,jlocal,3);
            end
        end
    end
end
end