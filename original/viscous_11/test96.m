function [terms_sf,terms_sj]=test96(c,ne,n1,n2,Two_chi_n1,...
    grad_T,term_persistent,coef_T,diff)

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
                
%                 sub1=-(c(ix,iy,1,e)^2*n1)^-1+((1-c(ix,iy,1,e))^2*n2)^-1;
%                 sub2=c(ix,iy,2,e)^2+c(ix,iy,3,e)^2;
%                 sub3=(c(ix,iy,1,e)*n1)^-1+((1-c(ix,iy,1,e))*n2)^-1-Two_chi_n1;
%                 sub4=c(ix,iy,4,e)+c(ix,iy,5,e);
%                 terms(e,ix,iy,1)=sub1*sub2;
%                 terms(e,ix,iy,2)=sub3*sub4;
%                 terms(e,ix,iy,3)=sub4;
%                 terms(e,ix,iy,4)=...
%                     2*((c(ix,iy,1,e)^3*n1)^-1+((1-c(ix,iy,1,e))^3*n2)^-1)...
%                     *sub2;
%                 terms(e,ix,iy,5)=sub1;
%                 terms(e,ix,iy,6)=sub3;
            end
        end
    end
    %==================================================================
elseif grad_T==1
%     terms=zeros(ne,3,3,6);

    terms_sf=zeros(ne,3,3,2);
    terms_sj=zeros(ne,3,3,6);
    
    for e=1:1:ne
%         [inx,~]=inxiny_elemental(e,ney);
        for ix=1:1:3
            for iy=1:1:3

                sub=c(ix,iy,2,e)^2+c(ix,iy,3,e)^2;
                c2=1-c(ix,iy,1,e);
                
                dot_Tc=coef_T(e,ix,iy,2)*c(ix,iy,2,e)...
                    +coef_T(e,ix,iy,3)*c(ix,iy,2,3);
                
                terms_sf(e,ix,iy,2)=c(ix,iy,4,e)+c(ix,iy,5,e);
                
                terms_sj(e,ix,iy,1)=(c(ix,iy,1,e)*n1)^-1+(c2*n2)^-1-Two_chi_n1(e,ix,iy);
                terms_sj(e,ix,iy,2)=-(c(ix,iy,1,e)^2*n1)^-1+(c2^2*n2)^-1;
                terms_sj(e,ix,iy,3)=2*dot_Tc...
                    +coef_T(e,ix,iy,1)*terms_sf(e,ix,iy,2);
                terms_sj(e,ix,iy,4)=...
                    2*((c(ix,iy,1,e)^3*n1)^-1+(c2^3*n2)^-1)...
                    *coef_T(e,ix,iy,1)*sub;
                
                terms_sf(e,ix,iy,1)=diff*(((log(c(ix,iy,1,e))+1)/n1...
                    -(log(c2)+1)/n2+Two_chi_n1(e,ix,iy)...
                    *(0.5-c(ix,iy,1,e)))*coef_T(e,ix,iy,4)...
                    +terms_sj(e,ix,iy,1)...
                    *(2*dot_Tc+coef_T(e,ix,iy,1)*terms_sf(e,ix,iy,2))...
                    +terms_sj(e,ix,iy,2)*coef_T(e,ix,iy,1)*sub...
                    +term_persistent(e,ix,iy)*(...
                    (0.5-c(ix,iy,1,e))*coef_T(e,ix,iy,4)-2*dot_Tc));
                
                terms_sj(e,ix,iy,5)=c(ix,iy,2,e);
                terms_sj(e,ix,iy,6)=c(ix,iy,3,e);
                
                
                %----------------------------
                
% %                 terms_sf(e,ix,iy,1)=terms_sf(e,ix,iy,1)*coef_T(e,ix,iy,1);
%                 
%                 
%                 if e==51 && ix==2 && iy==3
%                     warning('correct one')
% %                     term_persistent(e,ix,iy)*(...
% %                         (0.5-c(ix,iy,1,e))*coef_T(e,ix,iy,4)-2*dot_Tc)
% %                     diff
%                     terms_sf(e,ix,iy,1)
%                 end
                %--
                
%                 terms(e,ix,iy,3)=((c(ix,iy,1,e)^3*n1)^-1+(c2^3*n2)^-1)*sub;
                
                
                
                
%                 terms(e,ix,iy,4)=-(c(ix,iy,1,e)^2*n1)^-1+(c2^2*n2)^-1;
%                 terms(e,ix,iy,5)=coef_T(inx,ix)*terms(e,ix,iy,2)...
%                     +two_slope*c(ix,iy,2,e);
                
%                 terms(e,ix,iy,6)=(c(ix,iy,1,e)*n1)^-1+(c2*n2)^-1-Two_chi_n1(inx,ix);
                
%                 terms(e,ix,iy,1)=coef_T(inx,ix)*(terms(e,ix,iy,4)*sub...
%                     -2*term_persistent(inx,ix)*c(ix,iy,2,e))...
%                     +terms(e,ix,iy,6)*(coef_T(inx,ix)*terms(e,ix,iy,2)...
%                     +two_slope*c(ix,iy,2,e));
                
                
                
            end
        end
    end
end
end