function solution=structure_factor_2D(alpha,beta,fc,c,co)

solution=zeros(1,alpha);

solution(1)=fft2lensq(c,co,0,0);

angular=2*pi/beta;
radial=fc/(alpha-1.);

radius=0;
for iradius=2:1:alpha
    radius=radius+radial;
    angle=-angular;
    fft2lensq_=0;
    for iangle=1:1:beta
        angle=angle+angular;
        k1=radius*cos(angle);
        k2=radius*sin(angle);
        fft2lensq_=fft2lensq_+fft2lensq(c,co,k1,k2);
    end
    solution(iradius)=fft2lensq_/beta;
end

end

function solution=fft2lensq(c,co,k1,k2)

[M,N]=size(c);

c_=c-co*ones(M,N);

csum=0;
ssum=0;

twopi=2*pi;
k1M=k1/M;
k2N=k2/N;

for m=1:1:M
    for n=1:1:N
        csum=csum+c_(m,n)*cos(twopi*(k1M*(m-1)+k2N*(n-1)));
        ssum=ssum-c_(m,n)*sin(twopi*(k1M*(m-1)+k2N*(n-1)));
    end
end

solution=(csum^2+ssum^2)/(M*N);

end