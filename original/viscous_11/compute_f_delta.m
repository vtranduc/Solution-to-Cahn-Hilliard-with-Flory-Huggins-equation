function [c_domain,f_delta]=compute_f_delta(T,n1,n2,entropy,tol,...
    ci_ave,nDomain)

c_domain=linspace(tol,1-tol,nDomain);
chi=compute_chi(T,entropy);
c2=1-ci_ave;
f_co=ci_ave*log(ci_ave)/n1+c2*log(c2)/n2+ci_ave*c2*chi/n1;
df_dc_co=log(ci_ave)/n1-log(c2)/n2+1/n1-1/n2+chi*(1-2*ci_ave)/n1;
f_delta=zeros(1,nDomain);
for i=1:1:nDomain
    c2=1-c_domain(i);
    f_delta(i)=...
        c_domain(i)*log(c_domain(i))/n1+c2*log(c2)/n2+c_domain(i)*c2*chi/n1...
        -(c_domain(i)-ci_ave)*df_dc_co-f_co;
end

% spinodal=zeros(1,2);
% 
% chi=compute_chi(T,entropy);
% alpha=2*chi*n1*n2;
% beta=n1-n2-alpha;
% radicant=beta^2-4*alpha*n2;
% if abs(radicant)<=tol && radicant~=0
%     radicant=0;
% elseif radicant<0
%     error('Error detected in computation of spinodal curve\nTemperature=%d does not yield inflection points\nCritical temperature is %d',T,Tcutoff1)
% end
% spinodal(1)=(-beta-sqrt(radicant))/(2*alpha);
% spinodal(2)=(-beta+sqrt(radicant))/(2*alpha);

end