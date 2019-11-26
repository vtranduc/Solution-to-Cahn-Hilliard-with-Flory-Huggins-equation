function sol=get_phis()
% This function is not to be called repetitively, because it calls the
% function get_Hermitian_pol_coeffs_3d().
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