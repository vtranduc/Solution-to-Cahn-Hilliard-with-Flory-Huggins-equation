function test11 

clear
clc

n1=1;
n2=10;

R=8.3144598;
T=0.5;
frac1=0.4;

frac2=1-frac1;

chi=1;

Gm=R*T*(n1*log(frac1)+n2*log(frac2)+n1*frac2*chi)

frac1_=linspace(0,1,100);

n=length(frac1_);

Gm_=zeros(1,n);
for i=1:1:n
    frac1=frac1_(i);
    frac2=1-frac1;
    Gm_(i)=R*T*(n1*log(frac1)+n2*log(frac2)+n1*frac2*chi);
end

plot(Gm_)

grid on


%==================================

n=100;

nex=100;
ney=100;

c_=linspace(0,1,n);

entropy=1;
T_theta=1;

chi=0.5-entropy*(1-T_theta/T);

kb=1.38064852d-23;

v=1/(nex*ney);

const=kb*T/v;

f=zeros(1,n);

const=1;

for i=1:1:n
    c=c_(i);
    f(i)=(c/n1)*log(c)+((1-c)/n2)*log(1-c)+chi*c*(1-c);
    f(i)=f(i)*const;
end

figure(3)

plot(c_,f)
grid on

%=================================================================

clear all

n=100;

n1=1
n2=3

R=8.3144598;
entropy=1;
T_theta=1;

T=0.5

chi=0.5-entropy*(1-T_theta/T)

% chi=0

c_=linspace(0.1,0.9,n);

f=zeros(1,n);



% chi=2.2
% 
% T=entropy*T_theta/(entropy+chi-0.5)

% n1=1;
% n2=1;

for i=1:1:n
    c1=c_(i);
    c2=1-c1;
    f(i)=c1/n1*log(c1)+c2/n2*log(1-c1)+c1*(1-c1)*chi;
%     f(i)=(1/n1-1/n2)+(1/n1)*log(c1)-(1/n2)*log(1-c1)+chi*(1-2*c1);
end

figure(2)

plot(c_,f)
grid on

% clear c
% 
% syms c
% 
% 
% vpasolve((1/n1-1/n2)+(1/n1)*log(c)-(1/n2)*log(1-c)+chi*(1-2*c)==0,c,'random',true)
% 
% vpasolve(c^2-1==0,c)

%=========================================================================

% clear all
% 
% myfun=@(c,n1,n2,chi) (1/n1-1/n2)+(1/n1)*log(c)-(1/n2)*log(1-c)+chi*(1-2*c);
% 
% n1=1;
% n2=1;
% 
% entropy=1;
% T_theta=1;
% 
% T=0.3704
% chi=0.5-entropy*(1-T_theta/T)
% 
% fun=@(c) myfun(c,n1,n2,chi)
% 
% fzero(fun,0.8)

%========================================================================

test_phase(1,20);

end