function spinodal=get_specific_spinodal_binodal(T,n1,n2,entropy,tol)

spinodal=zeros(1,2);

chi=get_chi(T,entropy,1);
alpha=2*chi*n1*n2;
beta=n1-n2-alpha;
radicant=beta^2-4*alpha*n2;
if abs(radicant)<=tol && radicant~=0
    radicant=0;
elseif radicant<0
    error('Error detected in computation of spinodal curve\nTemperature=%d does not yield inflection points\nCritical temperature is %d',T,Tcutoff1)
end
spinodal(1)=(-beta-sqrt(radicant))/(2*alpha);
spinodal(2)=(-beta+sqrt(radicant))/(2*alpha);

end