function solution=structure_factor_1D(data,freq_domain)

N=length(data);

data_=zeros(1,N);

ci_ave=mean(data);

for i=1:1:N
    data_(i)=data(i)-ci_ave;
end

solution=zeros(1,length(freq_domain));

twopiN=2*pi/N;

isol=0;
for k=freq_domain
    csum=0;
    ssum=0;
    kN=twopiN*k;
    for i=1:1:N
        csum=csum+data_(i)*cos((i-1)*kN);
        ssum=ssum-data_(i)*sin((i-1)*kN);
    end
    isol=isol+1;
    solution(isol)=(csum^2+ssum^2)/N;
end

end