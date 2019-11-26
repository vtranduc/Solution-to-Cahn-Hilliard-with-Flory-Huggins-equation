function [real,imaginary]=FT_real_imaginary_3d(c,k1,k2,k3,M,N,O,twopi)

real=0;
imaginary=0;
k1M=twopi*k1/M;
k2N=twopi*k2/N;
k3O=twopi*k3/O;

for m=1:1:M
    for n=1:1:N
        for o=1:1:O
            real=real+c(m,n,o)*cos(k1M*(m-1)+k2N*(n-1)+k3O*(o-1));
            imaginary=imaginary-c(m,n,o)*sin(k1M*(m-1)+k2N*(n-1)+k3O*(o-1));
        end
    end
end
divider=sqrt(M*N*O);
real=real/divider;
imaginary=imaginary/divider;
end