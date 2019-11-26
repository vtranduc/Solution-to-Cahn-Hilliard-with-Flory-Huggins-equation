function main_1d_test

clear
clc

diff=3400;

load_co=1;
save_c=0;
T=0.6;
frames=2000;

% load_co=2;
% save_c=1;
% T=0.5;
% frames=26;

force_adjustment=1;

%-------------------------------

ne=200;

ci=0.4;


dt=1.0*10^-6;



n1=1;n2=10;

n=ne+1;
ntwo=2*n;

fluc=0.005;

chi=get_chi(T,1,1);

w=[5/18 4/9 5/18];

dx=1/ne;

weights=generate_weightt_1d(dx);

tol=1.0*10^-6;tol_jk=20;

frequency_limit=10;sampling_pts=200;

%-----------------------------

nT=100;
nInitialGuessPts=100;
guess_accuracy=10^-3;

[xspinodal,yspinodal,xbinodal,ybinodal]=spinodal_binodal(...
    n1,n2,1,1,nT,tol,0,nInitialGuessPts,...
    guess_accuracy);

sampling_domain=linspace(0,frequency_limit,sampling_pts);


diffT=diff*T;

%--------------------------

if load_co==0
    co=generate_co_1d(ci,n,ntwo,fluc);
elseif load_co==1
    co=dlmread('rinchan');
elseif load_co==2
    co=zeros(1,ntwo);
    index=-1;
    adder=ci-fluc/n;
    for ijk=1:1:n
        index=index+2;
        co(index)=adder;
    end
    if mod(n,2)==1
        mid=(n-1)/2+1;
    else
        mid=n/2;
    end
    mid=2*(mid-1)+1;
    co(mid)=co(mid)+fluc;
elseif load_co==3
    co=zeros(1,ntwo);
    index=-1;
    for ijk=1:1:n
        index=index+2;
        co(index)=ci;
    end
end

figure(1)
set(gcf, 'Position', [0, 0, 1280, 720])
c_abs=extract_abs_c(co,n);

structure_factor=compute_structure_factor_1d(c_abs,n,ci,sampling_domain,sampling_pts);

time=0;
sampled_time=0;

fourier_trasform_max_log=log(max(structure_factor));
plotter(c_abs,structure_factor,frequency_limit,fourier_trasform_max_log,sampled_time,n,xspinodal,yspinodal,xbinodal,ybinodal,ci,T)
drawnow
% caps=getframe;

%---------------------------------

c=co;
for iTime=1:1:frames
    
    coo=co;co=c;
    
    err=inf;
    erro=inf;
    jk=0;
    
    while err>tol
        
        jk=jk+1;
        if jk==tol_jk
            error('Too many iterations')
        end
        
        sf=zeros(1,ntwo);
        sj=zeros(ntwo,ntwo);
        for e=1:1:ne
            gbfs=2*e-1:1:2*e+2;
            cgbfs=c(gbfs);
            conc_=compute_e_conc(cgbfs,weights);
            terms=compute_terms_1d(conc_,n1,n2,chi);
            cogbfs=co(2*e-1:1:2*e+2);
            conco_=compute_e_conco(cogbfs,weights);
            for ix=1:1:3
                for igbf=1:1:4
                    sf(gbfs(igbf))=sf(gbfs(igbf))+dx*w(ix)*(...
                        ...
                        weights(igbf,1,ix)*(conc_(1,ix)-conco_(ix))/dt...
                        -diffT*weights(igbf,1,ix)...
                        *(terms(1,ix)+terms(2,ix))...
                        +terms(3,ix)*weights(igbf,3,ix)...
                        ...
                        );

                    for jgbf=1:1:4
                        sj(gbfs(igbf),gbfs(jgbf))=sj(gbfs(igbf),gbfs(jgbf))+dx*w(ix)*(...
                            ...
                            weights(igbf,1,ix)*weights(jgbf,1,ix)/dt...
                            -diffT*weights(igbf,1,ix)*(weights(jgbf,1,ix)*terms(4,ix)+terms(5,ix)...
                            *(2*conc_(2,ix)*weights(jgbf,2,ix)+weights(jgbf,1,ix)*terms(3,ix))...
                            +terms(6,ix)*weights(jgbf,3,ix))...
                            +weights(igbf,3,ix)*weights(jgbf,3,ix)...
                            ...
                            ...
                            );
                    end
                end
            end
        end
        
        %-------------------------------------------
        
        sj(2,1)=0;sj(2,2)=1;sj(2,3)=0;sj(2,4)=0;sf(2)=0;
        
        sj(ntwo,ntwo-3)=0;sj(ntwo,ntwo-2)=0;sj(ntwo,ntwo-1)=0;sj(ntwo,ntwo)=1;sf(ntwo)=0;
        
        %--------------------------------------------
        
        c_=sj\-sf';
        
        
        
        %---------------------- Apply well bc
        c_(2)=0;
        c_(ntwo)=0;
        %-----------------
        
        err=sqrt(sum(c_.^2))
        
        if err>erro || isnan(err)
            jk=-1;
            warning('Diverging')
            break
        else
            erro=err;
        end
        c=c+c_';
    end
    
    if jk==-1
        dt=dt/2;
        c=co;co=coo;
        continue
    end
    
    time=time+dt;
    sampled_time(iTime+1)=time;
    
    if jk<=3
        warning('Stepping up')
        dt=dt*1.5;
    end

    c_abs=extract_abs_c(c,n);
    
    %===============
    
    if force_adjustment==1
        shifted_ci=mean(c_abs);
        adjustment=shifted_ci-ci;
        index=-1;
        for i=1:1:n
            index=index+2;
            c(index)=c(index)-adjustment;
        end
        c_abs=extract_abs_c(c,n);
        if max(c_abs)>=1 || min(c_abs)<=0
            error('c has fallen outside the domain due to adjustment')
        end
    end
    
    %===============
    
    if max(c_abs)>=1
        error('max reached')
    end
    
    structure_factor=compute_structure_factor_1d(c_abs,n,ci,sampling_domain,sampling_pts);
    
    fourier_trasform_max_log(iTime+1)=log(max(structure_factor));
    
    plotter(c_abs,structure_factor,frequency_limit,fourier_trasform_max_log,sampled_time,n,xspinodal,yspinodal,xbinodal,ybinodal,ci,T)
%     plotter(c_abs,abs(fft(c_abs)),n,xspinodal,yspinodal,xbinodal,ybinodal,ci,T)
    drawnow
%     caps(iTime+1)=getframe;
    
end

% c_abs=extract_abs_c(c,n);
% 
% structure_factor=compute_structure_factor_1d(c_abs,n,ci,sampling_domain,sampling_pts);
% 
% fourier_trasform_max_log=log(max(structure_factor));
% plotter(c_abs,structure_factor,frequency_limit,fourier_trasform_max_log,n,xspinodal,yspinodal,xbinodal,ybinodal,ci,T)
% drawnow

if save_c==1
    dlmwrite('rinchan',c);
end


% movie(caps)

end

function sol=compute_structure_factor_1d(c_abs,n,ci,sampling_domain,sampling_pts)
sol=zeros(1,sampling_pts);
cn=c_abs-ci*ones(1,n);
for i=1:1:sampling_pts
%     summation=0;
    csum=0;
    ssum=0;
    term1=2*pi*sampling_domain(i)/n;
    for j=0:1:n-1
        term2=term1*j;
%         summation=summation+cn(j+1)...
%             *(cos(term2)^2+sin(term2)^2);
        csum=csum+cn(j+1)*cos(term2);
        ssum=ssum-cn(j+1)*sin(term2);
    end
    sol(i)=csum^2+ssum^2;
end
sol=sol/n;
end

function co=generate_co_1d(ci,n,ntwo,fluc)
co=zeros(1,ntwo);
z1=2*fluc;
z2=0;
for i=1:2:ntwo-1
    fluctuation=(z1*rand(1,1)-1)*fluc;
    co(i)=ci+fluctuation;
    z2=z2+fluctuation;
end
z1=z2/n;
z2=0;
for i=1:2:ntwo-1
    co(i)=co(i)-z1;
    z2=z2+co(i);
end
end

function terms=compute_terms_1d(conc_,n1,n2,chi)
terms=zeros(6,3);
for ix=1:1:3
    sub_term1=-1/(conc_(1,ix)^2*n1)+1/((1-conc_(1,ix))^2*n2);
    sub_term2=conc_(2,ix)^2;
    sub_term3=1/(conc_(1,ix)*n1)+1/((1-conc_(1,ix))*n2)-2*chi;
    sub_term4=conc_(3,ix);
    terms(1,ix)=sub_term1*sub_term2;
    terms(2,ix)=sub_term3*sub_term4;
    terms(3,ix)=sub_term4;
    terms(4,ix)=2*(1/(conc_(1,ix)^3*n1)+1/((1-conc_(1,ix))^3*n2))*sub_term2;
    terms(5,ix)=sub_term1;
    terms(6,ix)=sub_term3;
end

end

function conco_=compute_e_conco(cogbfs,weights)
conco_=zeros(1,3);
for ix=1:1:3
    for ibf=1:1:4
        conco_(ix)=conco_(ix)+...
            cogbfs(ibf)*weights(ibf,1,ix);
    end
end
end

function conc_=compute_e_conc(cgbfs,weights)
conc_=zeros(3,3);
for iorder=1:1:3
    for ix=1:1:3
        for ibf=1:1:4
            conc_(iorder,ix)=conc_(iorder,ix)+...
                cgbfs(ibf)*weights(ibf,iorder,ix);
        end
    end
end
end

function solution=generate_weightt_1d(dx)
solution=zeros(4,3,3);
gps=[(-sqrt(3/5)+1)/2 0.5 (sqrt(3/5)+1)/2];
for igps=1:1:3
    solution(1,1,igps)=basis(gps(igps),0,0,0);
    solution(2,1,igps)=basis(gps(igps),0,1,0)*dx;
    solution(3,1,igps)=basis(gps(igps),1,0,0);
    solution(4,1,igps)=basis(gps(igps),1,1,0)*dx;
    solution(1,2,igps)=basis(gps(igps),0,0,1)/dx;
    solution(2,2,igps)=basis(gps(igps),0,1,1);
    solution(3,2,igps)=basis(gps(igps),1,0,1)/dx;
    solution(4,2,igps)=basis(gps(igps),1,1,1);
    solution(1,3,igps)=basis(gps(igps),0,0,2)/(dx^2);
    solution(2,3,igps)=basis(gps(igps),0,1,2)/dx;
    solution(3,3,igps)=basis(gps(igps),1,0,2)/(dx^2);
    solution(4,3,igps)=basis(gps(igps),1,1,2)/dx;
end
end

function plotter(c_abs,structure_factor,frequency_limit,fourier_trasform_max_log,sampled_time,n,xspinodal,yspinodal,xbinodal,ybinodal,ci,T)
subplot(2,2,1)
plot(xspinodal,yspinodal,xbinodal,ybinodal,ci,T,'o')
grid on
grid minor
% take=zeros(1,n);
% for i=1:1:n
%     take(i)=c((i-1)*2+1);
% end
subplot(2,2,2)
plot(linspace(0,1,n),c_abs)
axis([0 1 0. 1])
grid on
grid minor
subplot(2,2,3)
plot(linspace(0,frequency_limit,length(structure_factor)),structure_factor)
grid on
grid minor
subplot(2,2,4)
plot(sampled_time,fourier_trasform_max_log,'*')
grid on
grid minor
end

function sol=extract_abs_c(c,n)
sol=zeros(1,n);
for i=1:1:n
    sol(i)=c((i-1)*2+1);
end
end