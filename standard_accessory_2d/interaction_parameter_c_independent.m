function sol=interaction_parameter_c_independent(coef_T,chi_a,chi_b,type)

[ne,~,~,~]=size(coef_T);
sol=zeros(ne,3,3,2);

if type==1
    for e=1:1:ne
        for ix=1:1:3
            for iy=1:1:3
                sol(e,ix,iy,1)=chi_a+chi_b/coef_T(e,ix,iy,1);
                sol(e,ix,iy,2)=-chi_b/(coef_T(e,ix,iy,1)^2);
            end
        end
    end
end

end