function terms=compute_terms_2d(c,ne,n1,n2,Two_chi_n1,ney,...
    grad_T,term_persistent,two_slope,coef_T)

if grad_T==0
    terms=zeros(ne,3,3,6);
    for e=1:1:ne
%         [inx,~]=inxiny_elemental(e,ney);
        for ix=1:1:3
            for iy=1:1:3
                sub1=-(c(ix,iy,1,e)^2*n1)^-1+((1-c(ix,iy,1,e))^2*n2)^-1;
                sub2=c(ix,iy,2,e)^2+c(ix,iy,3,e)^2;
                sub3=(c(ix,iy,1,e)*n1)^-1+((1-c(ix,iy,1,e))*n2)^-1-Two_chi_n1;
                sub4=c(ix,iy,4,e)+c(ix,iy,5,e);
                terms(e,ix,iy,1)=sub1*sub2;
                terms(e,ix,iy,2)=sub3*sub4;
                terms(e,ix,iy,3)=sub4;
                terms(e,ix,iy,4)=...
                    2*((c(ix,iy,1,e)^3*n1)^-1+((1-c(ix,iy,1,e))^3*n2)^-1)...
                    *sub2;
                terms(e,ix,iy,5)=sub1;
                terms(e,ix,iy,6)=sub3;
            end
        end
    end
elseif grad_T==1
    terms=zeros(ne,3,3,6);
    for e=1:1:ne
        [inx,~]=inxiny_elemental(e,ney);
        for ix=1:1:3
            for iy=1:1:3
                
                sub=c(ix,iy,2,e)^2+c(ix,iy,3,e)^2;
                c2=1-c(ix,iy,1,e);
                
                terms(e,ix,iy,2)=c(ix,iy,4,e)+c(ix,iy,5,e);
                terms(e,ix,iy,3)=((c(ix,iy,1,e)^3*n1)^-1+(c2^3*n2)^-1)*sub;
                
                terms(e,ix,iy,4)=-(c(ix,iy,1,e)^2*n1)^-1+(c2^2*n2)^-1;
                terms(e,ix,iy,5)=coef_T(inx,ix)*terms(e,ix,iy,2)...
                    +two_slope*c(ix,iy,2,e);
                terms(e,ix,iy,6)=(c(ix,iy,1,e)*n1)^-1+(c2*n2)^-1-Two_chi_n1(inx,ix);
                
                terms(e,ix,iy,1)=coef_T(inx,ix)*(terms(e,ix,iy,4)*sub...
                    -2*term_persistent(inx,ix)*c(ix,iy,2,e))...
                    +terms(e,ix,iy,6)*(coef_T(inx,ix)*terms(e,ix,iy,2)...
                    +two_slope*c(ix,iy,2,e));
                
                
                
            end
        end
    end
end
end