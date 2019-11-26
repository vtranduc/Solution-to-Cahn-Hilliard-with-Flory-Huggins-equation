function k=get_characteristic_frequency(diff,ci,T,n1,n2,entropy)

% diff=10000;
% T=0.5;
% ci=0.5;
% n1=1;
% n2=10;
% entropy=1;

chi=0.5-entropy*(1-1/T);
d2f_dc2=(ci*n1)^-1+((1-ci)*n2)^-1-2*chi/n1;
if d2f_dc2>0
    error('This is not within unstable region!')
end
k=sqrt(-0.5*diff*T*d2f_dc2)/(2*pi);
end