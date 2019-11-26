function [terms_sf,terms_sj]=terms_2d(c,ne,n1,n2,Two_chi_n1,...
    grad_T,coef_T,diff,terms_assist,thermophoresis,wetting)

if wetting==0
    if grad_T==0
        %=====================MUST FIX THIS==============================
    %     terms=zeros(ne,3,3,6);

        terms_sf=zeros(ne,3,3,2);
        terms_sj=zeros(ne,3,3,6);

        for e=1:1:ne
    %         [inx,~]=inxiny_elemental(e,ney);
            for ix=1:1:3
                for iy=1:1:3

                    sub=c(ix,iy,2,e)^2+c(ix,iy,3,e)^2;
                    c2=1-c(ix,iy,1,e);

                    terms_sf(e,ix,iy,2)=c(ix,iy,4,e)+c(ix,iy,5,e);

                    terms_sj(e,ix,iy,1)=(c(ix,iy,1,e)*n1)^-1+(c2*n2)^-1-Two_chi_n1;
                    terms_sj(e,ix,iy,2)=-(c(ix,iy,1,e)^2*n1)^-1+(c2^2*n2)^-1;
                    terms_sj(e,ix,iy,3)=terms_sf(e,ix,iy,2);
                    terms_sj(e,ix,iy,4)=...
                        2*((c(ix,iy,1,e)^3*n1)^-1+(c2^3*n2)^-1)*sub;
                    terms_sj(e,ix,iy,5)=c(ix,iy,2,e);
                    terms_sj(e,ix,iy,6)=c(ix,iy,3,e);

                    terms_sf(e,ix,iy,1)=coef_T*(terms_sj(e,ix,iy,1)*terms_sj(e,ix,iy,3)...
                        +terms_sj(e,ix,iy,2)*sub);
                end
            end
        end
        %==================================================================
    elseif grad_T==1

        if thermophoresis==0

            terms_sf=zeros(ne,3,3,4);
            terms_sj=zeros(ne,3,3,8);

            for e=1:1:ne

                for ix=1:1:3
                    for iy=1:1:3

                        sub1=c(ix,iy,2,e)^2+c(ix,iy,3,e)^2;
                        c2=1-c(ix,iy,1,e);

                        sub2=diff*(...
                            (log(c(ix,iy,1,e))+1)/n1...
                            -(log(c2)+1)/n2+Two_chi_n1(e,ix,iy)...
                            *(0.5-c(ix,iy,1,e))...
                            +terms_assist(e,ix,iy)*(c(ix,iy,1,e)-0.5)...
                            );
                        terms_sf(e,ix,iy,2)=sub2*coef_T(e,ix,iy,2);
                        terms_sf(e,ix,iy,3)=sub2*coef_T(e,ix,iy,3);
                        terms_sf(e,ix,iy,4)=c(ix,iy,4,e)+c(ix,iy,5,e);

                        terms_sj(e,ix,iy,1)=(c(ix,iy,1,e)*n1)^-1+(c2*n2)^-1-Two_chi_n1(e,ix,iy);
                        terms_sj(e,ix,iy,2)=-(c(ix,iy,1,e)^2*n1)^-1+(c2^2*n2)^-1;
                        terms_sj(e,ix,iy,3)=(coef_T(e,ix,iy,2)*c(ix,iy,2,e)...
                            +coef_T(e,ix,iy,3)*c(ix,iy,2,3))...
                            +coef_T(e,ix,iy,1)*terms_sf(e,ix,iy,4);

                        terms_sj(e,ix,iy,4)=diff...
                            *2*((c(ix,iy,1,e)^3*n1)^-1+(c2^3*n2)^-1)...
                            *coef_T(e,ix,iy,1)*sub1;

                        sub3=diff*(terms_sj(e,ix,iy,1)+terms_assist(e,ix,iy));

                        terms_sj(e,ix,iy,5)=sub3*coef_T(e,ix,iy,2);
                        terms_sj(e,ix,iy,6)=sub3*coef_T(e,ix,iy,3);
                        sub4=2*coef_T(e,ix,iy,1);
                        terms_sj(e,ix,iy,7)=sub4*c(ix,iy,2,e);
                        terms_sj(e,ix,iy,8)=sub4*c(ix,iy,3,e);

                        dot_Tc=coef_T(e,ix,iy,2)*c(ix,iy,2,e)...
                            +coef_T(e,ix,iy,3)*c(ix,iy,2,3);

                        terms_sf(e,ix,iy,1)=diff*(terms_sj(e,ix,iy,1)*(...
                            dot_Tc+coef_T(e,ix,iy,1)*terms_sf(e,ix,iy,4))...
                            ...
                            +terms_sj(e,ix,iy,2)*coef_T(e,ix,iy,1)*sub1...
                            +terms_assist(e,ix,iy)*dot_Tc);

                        terms_sj(e,ix,iy,1)=diff*terms_sj(e,ix,iy,1)*coef_T(e,ix,iy,1);
                        terms_sj(e,ix,iy,2)=diff*terms_sj(e,ix,iy,2); 


                    end
                end
            end

        elseif thermophoresis~=0
            terms_sf=zeros(ne,3,3,4);
            terms_sj=zeros(ne,3,3,10);
            for e=1:1:ne
                for ix=1:1:3
                    for iy=1:1:3

                        sub1=c(ix,iy,2,e)^2+c(ix,iy,3,e)^2;
                        c2=1-c(ix,iy,1,e);

                        sub2=diff*(...
                            (log(c(ix,iy,1,e))+1)/n1...
                            -(log(c2)+1)/n2+Two_chi_n1(e,ix,iy)...
                            *(0.5-c(ix,iy,1,e))...
                            +terms_assist(e,ix,iy)*(c(ix,iy,1,e)-0.5)...
                            ...
                            ...%==============THERMOPHORESIS============
                            )+thermophoresis*c(ix,iy,1,e)*c2...
                            ...%=========================================
                            ...
                            ;


                        terms_sf(e,ix,iy,2)=sub2*coef_T(e,ix,iy,2);
                        terms_sf(e,ix,iy,3)=sub2*coef_T(e,ix,iy,3);
                        terms_sf(e,ix,iy,4)=c(ix,iy,4,e)+c(ix,iy,5,e);

                        terms_sj(e,ix,iy,1)=(c(ix,iy,1,e)*n1)^-1+(c2*n2)^-1-Two_chi_n1(e,ix,iy);
                        terms_sj(e,ix,iy,2)=-(c(ix,iy,1,e)^2*n1)^-1+(c2^2*n2)^-1;
                        terms_sj(e,ix,iy,3)=(coef_T(e,ix,iy,2)*c(ix,iy,2,e)...
                            +coef_T(e,ix,iy,3)*c(ix,iy,2,3))...
                            +coef_T(e,ix,iy,1)*terms_sf(e,ix,iy,4);

                        terms_sj(e,ix,iy,4)=diff...
                            *2*((c(ix,iy,1,e)^3*n1)^-1+(c2^3*n2)^-1)...
                            *coef_T(e,ix,iy,1)*sub1;

                        sub3=diff*(terms_sj(e,ix,iy,1)+terms_assist(e,ix,iy));

                        terms_sj(e,ix,iy,5)=sub3*coef_T(e,ix,iy,2);
                        terms_sj(e,ix,iy,6)=sub3*coef_T(e,ix,iy,3);
                        sub4=2*coef_T(e,ix,iy,1);
                        terms_sj(e,ix,iy,7)=sub4*c(ix,iy,2,e);
                        terms_sj(e,ix,iy,8)=sub4*c(ix,iy,3,e);

                        dot_Tc=coef_T(e,ix,iy,2)*c(ix,iy,2,e)...
                            +coef_T(e,ix,iy,3)*c(ix,iy,2,3);

                        terms_sf(e,ix,iy,1)=diff*(terms_sj(e,ix,iy,1)*(...
                            dot_Tc+coef_T(e,ix,iy,1)*terms_sf(e,ix,iy,4))...
                            ...
                            +terms_sj(e,ix,iy,2)*coef_T(e,ix,iy,1)*sub1...
                            +terms_assist(e,ix,iy)*dot_Tc);

                        terms_sj(e,ix,iy,1)=diff*terms_sj(e,ix,iy,1)*coef_T(e,ix,iy,1);
                        terms_sj(e,ix,iy,2)=diff*terms_sj(e,ix,iy,2);
                        %=============THERMOPHORESIS================
                        sub5=sub3+thermophoresis*(1-2*c(ix,iy,1,e));
                        terms_sj(e,ix,iy,9)=sub5*coef_T(e,ix,iy,2);
                        terms_sj(e,ix,iy,10)=sub5*coef_T(e,ix,iy,3);
                        %===========================================
                    end
                end
            end
        end
    end
    
elseif wetting==1
    terms_sf=zeros(ne,3,3,3);
    terms_sj=zeros(ne,3,3,3);
    for e=1:1:ne
        for ix=1:1:3
            for iy=1:1:3
                c2=1-c(ix,iy,1,e);
                terms_sj(e,ix,iy,3)=coef_T*((c(ix,iy,1,e)*n1)^-1+(c2*n2)^-1-Two_chi_n1);
                terms_sf(e,ix,iy,1)=terms_sj(e,ix,iy,3)*c(ix,iy,2,e);
                terms_sf(e,ix,iy,2)=terms_sj(e,ix,iy,3)*c(ix,iy,3,e);
                terms_sf(e,ix,iy,3)=c(ix,iy,4,e)+c(ix,iy,5,e);
                sub=coef_T*(-(c(ix,iy,1,e)^2*n1)^-1+(c2^2*n2)^-1);
                terms_sj(e,ix,iy,1)=sub*c(ix,iy,2,e);
                terms_sj(e,ix,iy,2)=sub*c(ix,iy,3,e);
            end
        end
    end
end


end