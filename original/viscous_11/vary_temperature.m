function [temperature,Two_chi_n1,terms_assist]=vary_temperature(T,type,spec,time,...
    ne,nex,ney,xlen,ylen,entropy,n1,...
    distribution_type,T_distribution_spec)

gp=[(-sqrt(3/5)+1)/2 0.5 (sqrt(3/5)+1)/2];
dx=xlen/nex;
gp=gp*dx;
temperature=zeros(ne,3,3,3);

if type==1
    
    %----------------------
%     varying_end=rate_coef*time+T(1);
    if time<spec
        varying_end=((T(2)-T(1))/spec)*time+T(1);
    elseif time>=spec
        varying_end=T(2);
    end
    
    temperature=T_characterization(distribution_type,ne,[T(1) varying_end],...
        T_distribution_spec,nex,ney,xlen,ylen);
    
%     slope=(varying_end-T(1))/xlen;
%     %----------------------
%     
%     
%     for e=1:1:ne
%         [xth,~]=inxiny_elemental(e,ney);
%         ifirst=(xth-1)*dx;
%         for ix=1:1:3
%             pos=ifirst+gp(ix);
%             temperature(e,ix,1,1)=slope*pos+T(1);
%             temperature(e,ix,1,2)=slope;
%             for iy=2:1:3
%                 temperature(e,ix,iy,1)=temperature(e,ix,1,1);
%                 temperature(e,ix,iy,2)=slope;
%             end
%         end
%     end
%     
end

Two_chi_n1=get_two_chi_n1(ne,temperature,entropy,n1);

terms_assist=assist_global(...
    1,ne,temperature,entropy,n1);

end