function chi=dimensionless_chi(T,entropy)
chi=zeros(1,length(T));
for i=1:1:length(T)
    chi(i)=0.5-entropy*(1-1/T(i));
end
end