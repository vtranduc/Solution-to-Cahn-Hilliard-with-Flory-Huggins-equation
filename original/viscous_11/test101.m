function [sol,test]=test101(data,fc,nfrequencies)

sol=zeros(nfrequencies,nfrequencies);
frequency_domain=linspace(-fc,fc,nfrequencies);
[M,N]=size(data);
MN=M*N;

test=zeros(2,nfrequencies^2);
itest=0;

for ik1=1:1:nfrequencies
    for ik2=1:1:nfrequencies
        csum=0;
        ssum=0;
        term1=2*pi*frequency_domain(ik1)/M;
        term2=2*pi*frequency_domain(ik2)/N;
        for m=0:1:M-1
            for n=0:1:N-1
                term=term1*m+term2*n;
                csum=csum+data(m+1,n+1)*cos(term);
                ssum=ssum-data(m+1,n+1)*sin(term);
            end
        end
        sol(ik1,ik2)=(csum^2+ssum^2)/MN;
        
        itest=itest+1;
        test(1,itest)=sqrt(frequency_domain(ik1)^2+frequency_domain(ik2)^2);
        test(2,itest)=sol(ik1,ik2);
        
    end
end

end