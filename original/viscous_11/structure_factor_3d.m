function solution=structure_factor_3d(alpha_plus,beta,gamma,fc,c,co)

solution=zeros(1,alpha_plus);

solution(1)=fft2lensq(c,co,0,0,0);

twopi=2*pi; %Combine

angular=twopi/beta;
radial=fc/(alpha_plus-1);
elevator=pi/(gamma+1);

radius=0;

total_pts=beta*gamma+2;

% figure(3415153)
% hold on

for iradius=2:1:alpha_plus
    radius=radius+radial;
    
    elevation=0;
    fft2lensq_=0;
    
    for ielevation=1:1:gamma
        elevation=elevation+elevator;
        
        angle=-angular;
        
        for iangle=1:1:beta
            angle=angle+angular;
            k1=radius*sin(elevation)*cos(angle);
            k2=radius*sin(elevation)*sin(angle);
            k3=radius*cos(elevation);
            
%             if iradius==3
%                 plot3(k1,k2,k3,'o')
%             end
            
            fft2lensq_=fft2lensq_+fft2lensq(c,co,k1,k2,k3);
            
        end
    end
    
%     if iradius==3
%         plot3(0,0,radius,'x')
%         plot3(0,0,-radius,'x')
%     end
%     
%     if iradius==3
%         error('adfad')
%     end

    fft2lensq_=fft2lensq_+fft2lensq(c,co,0,0,radius);
    fft2lensq_=fft2lensq_+fft2lensq(c,co,0,0,-radius);
    
    solution(iradius)=fft2lensq_/total_pts;
    
end

end

function solution=fft2lensq(c,co,k1,k2,k3)

[N,M,O]=size(c); %Bring this up

c_=c-co*ones(N,M,O); %Bring this up

csum=0;
ssum=0;

twopi=2*pi; % Bring this up

k2N=k2/N;
k1M=k1/M;
k3O=k3/O;

for n=1:1:N
    for m=1:1:M
        for o=1:1:O
            
            term=twopi*(k2N*(n-1)+k1M*(m-1)+k3O*(o-1));
            
            csum=csum+c_(n,m,o)*cos(term);
            ssum=ssum-c_(n,m,o)*sin(term);
        end
    end
end

solution=(csum^2+ssum^2)/(M*N*O);

end