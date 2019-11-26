function solution=structure_factor_1D_multiple(c_nodal,...
    freq_domain,isf1d)

n=length(isf1d);

solution=zeros(n,length(freq_domain));

for i=1:1:n
    solution(i,:)=structure_factor_1D(c_nodal(:,isf1d(i)),freq_domain);
end

end