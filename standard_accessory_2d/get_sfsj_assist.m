function sfsj_assist=get_sfsj_assist(ne,weights,wTerms,...
    coef_M, coef_T,viscosity,thermophoresis)

sfsj_assist=zeros(ne,3,3,16,3);

for e=1:1:ne
    for ix=1:1:3
        for iy=1:1:3
            for ilocal=1:1:16

                sfsj_assist(e,ix,iy,ilocal,1)=...
                    weights(ix,iy,ilocal,1)-viscosity*wTerms(ix,iy,ilocal);

                sfsj_assist(e,ix,iy,ilocal,2)=...
                    coef_M(e,ix,iy,2)*weights(ix,iy,ilocal,2)...
                    +coef_M(e,ix,iy,3)*weights(ix,iy,ilocal,3)...
                    +coef_M(e,ix,iy,1)*wTerms(ix,iy,ilocal);

                sfsj_assist(e,ix,iy,ilocal,3)=thermophoresis...
                    *(coef_T(e,ix,iy,2)*weights(ix,iy,ilocal,2)...
                    +coef_T(e,ix,iy,3)*weights(ix,iy,ilocal,3));

            end
        end
    end
end

end