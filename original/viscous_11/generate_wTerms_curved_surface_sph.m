function wTerms_surface=generate_wTerms_curved_surface_sph(n_eSurf,weights_surface)

wTerms_surface=zeros(n_eSurf(6),8,8,3,3);

% weights_surface=zeros(n_eSurf(6),8,8,7,3,3);

for e=1:1:n_eSurf(6)
    for iorientation=1:1:8
        for itype=1:1:8
            for dim1=1:1:3
                for dim2=1:1:3
                    wTerms_surface(e,iorientation,itype,dim1,dim2)=...
                        weights_surface(e,iorientation,itype,5,dim1,dim2)...
                        +weights_surface(e,iorientation,itype,6,dim1,dim2)...
                        +weights_surface(e,iorientation,itype,7,dim1,dim2);
                end
            end
        end
    end
end

end