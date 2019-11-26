function test14

clear
clc
subplot(2,1,1)
%==========================================================================
[xspinodal,yspinodal]=compute_spinodal_binodal();
plot(xspinodal,yspinodal,'.')
grid on
%==========================================================================
n1=1; %IT MUST BE 1!
n2=10;
entropy=1; %IT MUST BE 1!
T_theta=1; %IT MUST BE 1!
nT=100;
T_min=0.001; %Must be larger than zero. Should be as close to 0 as possible
nc_eachSide=100;
c_min=0.001; %Must be larger than zero. Should be as close to 0 as possible
c_max=0.9999999999999999999;
%Function-specific literals
tol=1.0e-6;
T_min_default_factor=0.0001; %Has to be less than 1
%==========================================================================


%sol=vpasolve([dfreeE_dc1(c1,n1,n2,chi_symbolic)==dfreeE_dc1(c2,n1,n2,chi_symbolic),m==dfreeE_dc1(c1,n1,n2,chi_symbolic),m*c1+b==freeE(c1,n1,n2,chi_symbolic),m*c2+b==freeE(c2,n1,n2,chi_symbolic)],[T c1 m b],guess);
hold on
plot(xspinodal(nT),yspinodal(nT),'ro')
hold off

%==========================================================================

T=yspinodal(nT-50)
syms c1

testT=[0.2 0.5 0.85];
itest=99; %=================
hold on
plot(xspinodal(itest),yspinodal(itest),'go')
hold off
testT=yspinodal(itest)
% xspinodal;
subplot(2,1,2)
for i=1:1:length(testT)
    if i==2
        hold on
    end
    T=testT(i);
    chi=get_chi(T,entropy,T_theta);
    % sol=vpasolve(dfreeE_dc1(c1,n1,n2,chi)==0,c1,[0.5,xspinodal(nT-1)])
    plot(xspinodal,freeE(xspinodal,n1,n2,chi),'.')
end
if i~=1
    hold off
end
grid on

jtest=nT+(nT-itest);

c1_=xspinodal(itest);
c2_=xspinodal(jtest);
z3=(freeE(c2_,n1,n2,chi)-freeE(c1_,n1,n2,chi))/...
    (c2_-c1_)
z1=c1_/2;
z2=(c2_+1)/2;
% z2=0.99999999999
z4=0;

for i=1:1:10
    [residual,jacobian]=...
        newton_raphson_frame(z1,z2,z3,z4,n1,n2,chi);
    adjust=jacobian\-residual;
    error=sqrt(sum(adjust.^2))
    z1=z1+adjust(1);
    z2=z2+adjust(2);
    z3=z3+adjust(3);
    z4=z4+adjust(4);
end

n=100;
domain=linspace(0,1,n);
range=zeros(1,n);
for i=1:1:n
    range(i)=z3*domain(i)+z4;
end
hold on
plot(domain,range)
hold off

% [residual,jacobian]=newton_raphson_frame(c1,c2,m,b,n1,n2,chi)

% sol=vpasolve(dfreeE_dc1(c1,n1,n2,chi)==0,c1,[0,xspinodal(nT-1)])
% hold on
% plot(sol,T,'x')
% hold off



% test15(yspinodal(nT-1))

% xbinodal=zeros(1,nT*2-1);
% ybinodal=yspinodal;
% 
% xbinodal(nT)=xspinodal(nT);
% syms c1 c2 m b
% T=yspinodal(nT-1)
% guess=[xspinodal(nT-1)/2 (xspinodal(nT+1)+1)/2 0.1 0]
% chi=get_chi(T,entropy,T_theta)
% disp('solving')
% sol=vpasolve([dfreeE_dc1(c1,n1,n2,chi)==dfreeE_dc1(c2,n1,n2,chi),m==dfreeE_dc1(c1,n1,n2,chi),m*c1+b==freeE(c1,n1,n2,chi),m*c2+b==freeE(c2,n1,n2,chi)],[c1 c2 m b],guess);
% disp('end solving')
% sol.c1
% sol.c2
% double([sol.c1 xbinodal(nT) sol.c2])

% disp('Initializing binodal')
% 
% xbinodal=zeros(1,nT*2-1);
% ybinodal=yspinodal;
% xbinodal(nT)=xspinodal(nT);
% syms c1 c2 m b
% guess=zeros(1,4);
% for k=1:1:nT-1
%     i=nT-k;
%     j=nT+k;
%     T=yspinodal(i);
%     chi=get_chi(T,entropy,T_theta);
%     guess(1)=xspinodal(i)/2;
%     guess(2)=(xspinodal(j)+1)/2;
%     guess(3)=(freeE(xspinodal(j),n1,n2,chi)-freeE(xspinodal(i),n1,n2,chi))/...
%         (xspinodal(j)-xspinodal(i));
% %     guess(3)=0;
%     disp('solving')
%     disp(k)
% 	sol=vpasolve([dfreeE_dc1(c1,n1,n2,chi)==dfreeE_dc1(c2,n1,n2,chi),...
%         m==dfreeE_dc1(c1,n1,n2,chi),m*c1+b==freeE(c1,n1,n2,chi),...
%         m*c2+b==freeE(c2,n1,n2,chi)],[c1 c2 m b],guess);
% 	disp('miku chan')
%     %Debating which method is better... Try statement or if statement?
%     %Perhaps if statement cuz you wanna check to make sure we are
%     %proceeding forward and check for imaginary number as well
% 	if length(sol.c1)~=1 || length(sol.c2)~=1 || length(sol.m)~=1 ||...
%             length(sol.b)~=1 || sol.c1>xbinodal(i+1) ||...
%             sol.c2<xbinodal(j-1) || sol.m<0 || isreal(sol.c1)==0 ||...
%             isreal(sol.c2)==0 || isreal(sol.m)==0 || isreal(sol.b)==0
%         warning('Temperature method is no longer working')
%         break
%     else
%         xbinodal(i)=sol.c1;
%         xbinodal(j)=sol.c2;
% 	end
% end
% 
% sol.c1
% sol.c2
%==========================================================================


end

function [residual,jacobian]=newton_raphson_frame(c1,c2,m,b,n1,n2,chi)
jacobian=zeros(4,4);
jacobian(1,1)=d2freeE_dc12(c1,n1,n2,chi);
jacobian(1,3)=-1;
jacobian(2,2)=d2freeE_dc12(c2,n1,n2,chi);
jacobian(2,3)=-1;
jacobian(3,1)=dfreeE_dc1(c1,n1,n2,chi)-m;
jacobian(3,3)=-c1;
jacobian(3,4)=-1;
jacobian(4,2)=dfreeE_dc1(c2,n1,n2,chi)-m;
jacobian(4,3)=-c2;
jacobian(4,4)=-1;

residual=zeros(4,1);
residual(1)=dfreeE_dc1(c1,n1,n2,chi)-m;
residual(2)=dfreeE_dc1(c2,n1,n2,chi)-m;
residual(3)=freeE(c1,n1,n2,chi)-m*c1-b;
residual(4)=freeE(c2,n1,n2,chi)-m*c2-b;
end

function solution=freeE(c1_,n1,n2,chi)
n=length(c1_);
solution=zeros(1,n);
for i=1:1:n
    c1=c1_(i);
    c2=1-c1;
    solution(i)=c1*log(c1)/n1+c2*log(c2)/n2+c1*c2*chi;
end
end

function solution=dfreeE_dc1(c1,n1,n2,chi)
solution=(1/n1-1/n2)+log(c1)/n1-log(1-c1)/n2+(1-2*c1)*chi;
end

function solution=d2freeE_dc12(c1,n1,n2,chi)
solution=1/(n1*c1)+1/(n2*(1-c1))-2*chi;
end