function [mag_spec,phase_spec]=test86(alpha,beta,fc,c,co)
mag_spec=zeros(1,alpha);
phase_spec=zeros(1,alpha);
[mag_spec(1),phase_spec(1)]=mag_phase(c,co,0,0);
angular=2*pi/beta;
radial=fc/(alpha-1.);
radius=0;
for iradius=2:1:alpha
    radius=radius+radial;
    angle=-angular;
    mag_phase_spec=[0 0];
    for iangle=1:1:beta
        angle=angle+angular;
        k1=radius*cos(angle);
        k2=radius*sin(angle);
        mag_phase_spec=mag_phase_spec+mag_phase(c,co,k1,k2);
    end
    mag_phase_spec=mag_phase_spec/beta;
    mag_spec(iradius)=mag_phase_spec(1);
    phase_spec(iradius)=mag_phase_spec(2);
end
end

function [mag,phase]=mag_phase(c,co,k1,k2)
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
mag=(csum^2+ssum^2)/(M*N);
phase=atan(ssum/csum);
end