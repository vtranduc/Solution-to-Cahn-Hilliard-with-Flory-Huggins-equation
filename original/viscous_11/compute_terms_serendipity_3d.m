function [sf_terms,sj_terms]=compute_terms_serendipity_3d(...
    ne,conc_,two_chi_n1,n1,n2,grad_T,coef_T,...
    conco_,dt,...
    diff,thermophoresis,sj_assist)

% === Compute terms =======================================================

% terms=zeros(ne,6,3,3,3);
sf_terms=zeros(ne,5,3,3,3);
sj_terms=zeros(ne,4,3,3,3);
if grad_T==0
%     size(conc_)
    parfor e=1:1:ne %Parfor must be feasible!!!
%         sf_terms(e,:,:,:,:)=compute_sf_sj_terms_serendipity_3d(e,conc_,two_chi_n1,n1,n2,coef_T);
        
        [sf,sj]=compute_sf_sj_terms_serendipity_3d(...
            e,conc_,two_chi_n1,n1,n2,coef_T,conco_,dt);
        sf_terms(e,:,:,:,:)=sf;
        sj_terms(e,:,:,:,:)=sj;
        
    end
elseif grad_T==1
    parfor e=1:1:ne
        
        [sf,sj]=compute_sf_sj_terms_serendipity_grad_T_3d(...
            e,conc_,two_chi_n1,n1,n2,coef_T,conco_,dt,...
            diff,thermophoresis,sj_assist);
        sf_terms(e,:,:,:,:)=sf;
        sj_terms(e,:,:,:,:)=sj;
        
    end
end

end

function [sf,sj]=compute_sf_sj_terms_serendipity_grad_T_3d(...
    e,conc_,two_chi_n1,n1,n2,coef_T,conco_,dt,...
    diff,thermophoresis,sj_assist)

sf=zeros(5,3,3,3);
sj=zeros(4,3,3,3);

for iz=1:1:3
    for iy=1:1:3
        for ix=1:1:3
            
%             c1=conc_(e,1,ix,iy,iz);
%             
%             cx=conc_(e,2,ix,iy,iz);
%             cy=conc_(e,3,ix,iy,iz);
%             cz=conc_(e,4,ix,iy,iz);
%             
%             c2=1-c1;
%             
%             chi=two_chi_n1(e,ix,iy,iz);
%             T=coef_T(e,ix,iy,iz,1);
%             Tx=coef_T(e,ix,iy,iz,2);
%             Ty=coef_T(e,ix,iy,iz,3);
%             Tz=coef_T(e,ix,iy,iz,4);
%             
%             f1=log(c1)/n1-log(c2)/n2+1/n1-1/n2+chi*(1-2*c1)/n1;
%             f2=1/(c1*n1)+1/(c2*n2)-2*chi/n1;
%             f3=-1/((c1^2)*n1)+1/((c2^2)*n2);
%             
%             diffT=diff*T;
            
            
            
            %-------------------
            c2=1-conc_(e,1,ix,iy,iz);
            
            sf(5,ix,iy,iz)=(conc_(e,1,ix,iy,iz)-conco_(e,ix,iy,iz))/dt;
            
            z1=1-2*conc_(e,1,ix,iy,iz);
            
            z2=diff*(log(conc_(e,1,ix,iy,iz))/n1-log(c2)/n2+1/n1-1/n2...
                +two_chi_n1(e,ix,iy,iz)*z1/2+...
                -z1/(n1*coef_T(e,ix,iy,iz,1)))...
                +thermophoresis*conc_(e,1,ix,iy,iz)*c2;
            
%             sj_assist/coef_T(e,ix,iy,iz,1)

            diffT=diff*coef_T(e,ix,iy,iz,1);
            
            z3=1/(conc_(e,1,ix,iy,iz)*n1)+1/(c2*n2)-two_chi_n1(e,ix,iy,iz);
            
            z4=diffT*z3;
            
            sf(1,ix,iy,iz)=coef_T(e,ix,iy,iz,2)*z2+conc_(e,2,ix,iy,iz)*z4;
            sf(2,ix,iy,iz)=coef_T(e,ix,iy,iz,3)*z2+conc_(e,3,ix,iy,iz)*z4;
            sf(3,ix,iy,iz)=coef_T(e,ix,iy,iz,4)*z2+conc_(e,4,ix,iy,iz)*z4;
            
            sf(4,ix,iy,iz)=...
                conc_(e,5,ix,iy,iz)...
                +conc_(e,6,ix,iy,iz)...
                +conc_(e,7,ix,iy,iz);
            
            z5=diffT*(-1/((conc_(e,1,ix,iy,iz)^2)*n1)+1/((c2^2)*n2));
            
            z6=diff*(z3+sj_assist/coef_T(e,ix,iy,iz,1))+thermophoresis*z1;
            
            
            sj(1,ix,iy,iz)=conc_(e,2,ix,iy,iz)*z5+coef_T(e,ix,iy,iz,2)*z6;
            sj(2,ix,iy,iz)=conc_(e,3,ix,iy,iz)*z5+coef_T(e,ix,iy,iz,3)*z6;
            sj(3,ix,iy,iz)=conc_(e,4,ix,iy,iz)*z5+coef_T(e,ix,iy,iz,4)*z6;
            
            sj(4,ix,iy,iz)=z4;
            
%             c2=1-conc_(e,1,ix,iy,iz);
%             
%             sj(4,ix,iy,iz)=coef_T*(1/(n1*conc_(e,1,ix,iy,iz))+1/(n2*c2)-two_chi_n1);
%             
%             sub=coef_T*(-1/(n1*conc_(e,1,ix,iy,iz)^2)+1/(n2*c2^2));
%             
%             sj(1,ix,iy,iz)=sub*conc_(e,2,ix,iy,iz);
%             sj(2,ix,iy,iz)=sub*conc_(e,3,ix,iy,iz);
%             sj(3,ix,iy,iz)=sub*conc_(e,4,ix,iy,iz);
%             
%             sf(1,ix,iy,iz)=sj(4,ix,iy,iz)*conc_(e,2,ix,iy,iz);
%             sf(2,ix,iy,iz)=sj(4,ix,iy,iz)*conc_(e,3,ix,iy,iz);
%             sf(3,ix,iy,iz)=sj(4,ix,iy,iz)*conc_(e,4,ix,iy,iz);
%             
%             sf(4,ix,iy,iz)=...
%                 conc_(e,5,ix,iy,iz)...
%                 +conc_(e,6,ix,iy,iz)...
%                 +conc_(e,7,ix,iy,iz);
%             
%             sf(5,ix,iy,iz)=(conc_(e,1,ix,iy,iz)-conco_(e,ix,iy,iz))/dt;
            
        end
    end
end
end

function [sf,sj]=compute_sf_sj_terms_serendipity_3d(...
    e,conc_,two_chi_n1,n1,n2,coef_T,conco_,dt)
sf=zeros(5,3,3,3);
sj=zeros(4,3,3,3);
for iz=1:1:3
    for iy=1:1:3
        for ix=1:1:3
            
            c2=1-conc_(e,1,ix,iy,iz);
            
            sj(4,ix,iy,iz)=coef_T*(1/(n1*conc_(e,1,ix,iy,iz))+1/(n2*c2)-two_chi_n1);
            
            sub=coef_T*(-1/(n1*conc_(e,1,ix,iy,iz)^2)+1/(n2*c2^2));
            
            sj(1,ix,iy,iz)=sub*conc_(e,2,ix,iy,iz);
            sj(2,ix,iy,iz)=sub*conc_(e,3,ix,iy,iz);
            sj(3,ix,iy,iz)=sub*conc_(e,4,ix,iy,iz);
            
            sf(1,ix,iy,iz)=sj(4,ix,iy,iz)*conc_(e,2,ix,iy,iz);
            sf(2,ix,iy,iz)=sj(4,ix,iy,iz)*conc_(e,3,ix,iy,iz);
            sf(3,ix,iy,iz)=sj(4,ix,iy,iz)*conc_(e,4,ix,iy,iz);
            
            sf(4,ix,iy,iz)=...
                conc_(e,5,ix,iy,iz)...
                +conc_(e,6,ix,iy,iz)...
                +conc_(e,7,ix,iy,iz);
            
            sf(5,ix,iy,iz)=(conc_(e,1,ix,iy,iz)-conco_(e,ix,iy,iz))/dt;
            
        end
    end
end
end