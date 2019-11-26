function [mag_spec,phase_spec]=fourier_analysis_2d(alpha,beta,fc,c,co,get_phase,nworkers)
mag_spec=zeros(1,alpha);
if get_phase==1
    phase_spec=zeros(1,alpha);
    [mag_spec(1),phase_spec(1)]=mag_phase(c,co,0,0);
elseif get_phase==0
    phase_spec=NaN;
    
    [M,N]=size(c);
    c_=c-co*ones(M,N);
    twopi=2*pi;
%     k1M=k1/M;
%     k2N=k2/N;

    mag_spec(1)=mag_only(c_,0,0,M,N,twopi);
end
if nworkers==1 || get_phase==1 %parallel computing is not available for get_phase==1 !!!!!!
    angular=2*pi/beta;
    radial=fc/(alpha-1.);
    radius=0;
end
if get_phase==1
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
elseif get_phase==0
    if nworkers==1
        for iradius=2:1:alpha
            radius=radius+radial;
            angle=-angular;
            magnitude=0;
            for iangle=1:1:beta
                angle=angle+angular;
                k1=radius*cos(angle);
                k2=radius*sin(angle);
                magnitude=magnitude+mag_only(c_,k1,k2,M,N,twopi);
            end
            mag_spec(iradius)=magnitude/beta;
        end
    elseif nworkers>1

        frequency_domain=linspace(0,fc,alpha);
        
        split=split_integer(alpha-1,nworkers);
        
        index_array=zeros(2,nworkers);
        index_array(1,1)=2;
        index_array(2,nworkers)=alpha;
        for i=2:1:nworkers
            index_array(2,i-1)=index_array(1,i-1)+split(i-1)-1;
            index_array(1,i)=index_array(2,i-1)+1;
        end
        angular_domain=linspace(0,twopi-twopi/beta,beta);
        cos_domain=cos(angular_domain);
        sin_domain=sin(angular_domain);
        mag_array=ones(nworkers,split(1));
        nLoad=split(1);
        parfor core=1:1:nworkers
            mag_array(core,:)=mag_only_parallel(...
                core,c_,frequency_domain,index_array,...
                beta,cos_domain,sin_domain,M,N,twopi,nLoad);
        end
        index=1;
        for core=1:1:nworkers
            for i=1:1:split(core)
                index=index+1;
                mag_spec(index)=mag_array(core,i);
            end
        end
    end
end
end

function sol=mag_only_parallel(core,c_,frequency_domain,index_array,...
    beta,cos_domain,sin_domain,M,N,twopi,nLoad)
sol=zeros(1,nLoad);
index=0;
for ik=index_array(1,core):1:index_array(2,core)
    index=index+1;
    magnitude=0;
    for ibeta=1:1:beta
        k1=frequency_domain(ik)*cos_domain(ibeta);
        k2=frequency_domain(ik)*sin_domain(ibeta);
        magnitude=magnitude+mag_only(c_,k1,k2,M,N,twopi);
    end
    sol(index)=magnitude/beta;
end

end

function mag=mag_only(c_,k1,k2,M,N,twopi)
% [M,N]=size(c);
% c_=c-co*ones(M,N);
csum=0;
ssum=0;
% twopi=2*pi;
k1M=k1/M;
k2N=k2/N;
for m=1:1:M
    for n=1:1:N
        csum=csum+c_(m,n)*cos(twopi*(k1M*(m-1)+k2N*(n-1)));
        ssum=ssum-c_(m,n)*sin(twopi*(k1M*(m-1)+k2N*(n-1)));
    end
end
mag=(csum^2+ssum^2)/(M*N);
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