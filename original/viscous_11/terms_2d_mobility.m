function [terms_sf,terms_sj]=terms_2d_mobility(c,ne,n1,chi,...
    coef_T,coef_n2,diff,thermophoresis,...
    nworkers,mobility_type,mobility_spec,ci_ave)

terms_sf=zeros(ne,3,3,3);
terms_sj=zeros(ne,3,3,3);

for e=1:1:ne
    for ix=1:1:3
        for iy=1:1:3
            
            c2=1-c(ix,iy,1,e);
            c1c2=c(ix,iy,1,e)*c2;
            
            diffT=diff*coef_T(e,ix,iy,1);
            
            sub1=1-2*c(ix,iy,1,e);
            
            df_dc=(log(c(ix,iy,1,e))+1)/n1-(log(c2)+1)/coef_n2(e,ix,iy,1)+chi(e,ix,iy,1)*sub1/n1;
            
            d2f_dc2=(c(ix,iy,1,e)*n1)^-1+(c2*coef_n2(e,ix,iy,1))^-1-2*chi(e,ix,iy,1)/n1;
            
            
            d3f_dc3=-(c(ix,iy,1,e)^2*n1)^-1+(c2^2*coef_n2(e,ix,iy,1))^-1;
            
            d2f_dcdT=sub1*chi(e,ix,iy,2)/n1;
            
            d3f_dc2dT=-2*chi(e,ix,iy,2)/n1;
            
            d2f_dcdn2=(log(c2)+1)/coef_n2(e,ix,iy,1)^2;
            
            d3f_dc2dn2=-(c2*coef_n2(e,ix,iy,1)^2)^-1;
            
%             d2f_dcdn2=0;
%             d3f_dc2dn2=0;
            
            
            sub2=diff*(df_dc+coef_T(e,ix,iy,1)*d2f_dcdT)+thermophoresis*c1c2;
            sub3=diffT*d2f_dc2;
            sub4=diffT*d2f_dcdn2;
            
            terms_sf(e,ix,iy,1)=sub2*coef_T(e,ix,iy,2)+sub3*c(ix,iy,2,e)+sub4*coef_n2(e,ix,iy,2);
            terms_sf(e,ix,iy,2)=sub2*coef_T(e,ix,iy,3)+sub3*c(ix,iy,3,e)+sub4*coef_n2(e,ix,iy,3);
            terms_sf(e,ix,iy,3)=c(ix,iy,4,e)+c(ix,iy,5,e);
            
            
            sub5=diff*(d2f_dc2+coef_T(e,ix,iy,1)*d3f_dc2dT)+thermophoresis*sub1;
            sub6=diffT*d3f_dc3;
            sub7=diffT*d3f_dc2dn2;

            terms_sj(e,ix,iy,1)=sub5*coef_T(e,ix,iy,2)+sub6*c(ix,iy,2,e)+sub7*coef_n2(e,ix,iy,2);
            terms_sj(e,ix,iy,2)=sub5*coef_T(e,ix,iy,3)+sub6*c(ix,iy,3,e)+sub7*coef_n2(e,ix,iy,3);
            
            terms_sj(e,ix,iy,3)=sub3;
            
            %----------------------April 2019-----------------------
            
            if mobility_type==1
                M=1;
                dM_dc=0;
            elseif mobility_type==2
                M=2-exp(mobility_spec*(c(ix,iy,1,e)-ci_ave)^2);
                dM_dc=-2*mobility_spec*(c(ix,iy,1,e)-ci_ave)*exp(mobility_spec*(c(ix,iy,1,e)-ci_ave)^2);
            else
                error('No such mobility_type is available')
            end
            
            terms_sf(e,ix,iy,1)=M*terms_sf(e,ix,iy,1);
            terms_sf(e,ix,iy,2)=M*terms_sf(e,ix,iy,2);
            
            
            
            
        end
    end
end
    
end
