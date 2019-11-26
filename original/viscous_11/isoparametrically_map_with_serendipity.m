function weights=isoparametrically_map_with_serendipity()
weights=zeros(ne,8,4,7,3,3,3);
phis=get_phis();

if parallel_computing==1
    for e=1
        
    end
elseif parallel_computing==0
    for e=1
        
    end
end


end

function sol=get_phis()
sol=zeros(8,4,10,3,3,3);
gps=[(-sqrt(3/5)+1)/2 0.5 (sqrt(3/5)+1)/2];
coeffs=get_Hermitian_pol_coeffs_3d();
for orientation=1:1:8
    for type=1:1:4
        for order=1:1:10
            for ix=1:1:3
                for iy=1:1:3
                    for iz=1:1:3
                        sol(orientation,type,order,ix,iy,iz)=...
                            compute_weight_specific(...
                            gps(ix),gps(iy),gps(iz),orientation,type,order,coeffs);
                    end
                end
            end
        end
    end
end
end