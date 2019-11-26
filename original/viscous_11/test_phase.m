function test_phase(n1,n2)

entropy=1;
Ttheta=1;
tol=1.0e-6;

nT=1000;
T_min=0.001;
T_min_default_factor=0.0001; %Has to be less than 1

%=====Validate inputs=================

if n1~=1 && n2~=1
    error('n1 or n2 must be 1!')
end

%=====Compute cut-off temp============

gamma=2*n1*n2;
A=gamma^2;
B=2*gamma*(-n2-n1);
C=(n2-n1)^2;

Xcut1=(-B+sqrt(B^2-4*A*C))/(2*A);
Xcut2=(-B-sqrt(B^2-4*A*C))/(2*A);

Tcut1=(1/Ttheta)*((Xcut1-0.5)/entropy+1)^-1
Tcut2=(1/Ttheta)*((Xcut2-0.5)/entropy+1)^-1

if Tcut1>1
    error('Cut-off temperature is larger than 1!')
end

%======Find out temperature range to be evaluated

if T_min<Tcut1 && T_min>0
    T_range=linspace(T_min,Tcut1,nT);
else
    warning('T_min specified is outside acceptable range. Other value is used')
    T_range=linspace(Tcut1*T_min_default_factor,Tcut1,nT);
end

%======Find inflection points===========

spinodal=zeros(nT,2);

for i=1:1:nT-1
    T=T_range(i);
    X=0.5-entropy*(1-Ttheta/T);
    alpha=2*X*n1*n2;
    beta=n1-n2-alpha;
    radical=beta^2-4*alpha*n2;
    if abs(radical)<tol
        %ARE YOU SURE????? MUST RECALL WHEN PHASE DIAGRAM "FLIPS"!!!!!!
%         error('Temperature is out of range!')
        radical=0;
    elseif radical<0
        error('Unreasonable temperature has been detected!')
    end
    radical=sqrt(radical);
    c1=(-beta+radical)/(2*alpha);
    c2=(-beta-radical)/(2*alpha);
    spinodal(i,:)=[c1 c2];
end

spinodal(nT,2)=-beta/(2*alpha);

xspinodal=[spinodal(:,2); spinodal(nT-1:-1:1,1)];
yspinodal=[T_range T_range(nT-1:-1:1)];

%=======Find binodal line========================

syms x1 x2 m b

binodal=zeros(nT-1,2);

for i=nT-1:-1:1
    T=T_range(i);
    guess1=spinodal(i,2)/2;
    guess2=(spinodal(i,1)+1)/2;
    guess3=(freeE(spinodal(i,1),n1,n2,X)-freeE(spinodal(i,2),n1,n2,X))...
        /(spinodal(i,1)-spinodal(i,2));
    
    X=0.5-entropy*(1-Ttheta/T); %COULD HAVE BEEN SAVED IN MEMORY
    
    try
        sol=vpasolve([(1/n1-1/n2)+(log(x1)/n1-log(1-x1)/n2)+(1-2*x1)*X==(1/n1-1/n2)+(log(x2)/n1-log(1-x2)/n2)+(1-2*x2)*X,m*x1+b==x1*log(x1)/n1+(1-x1)*log(1-x1)/n2+x1*(1-x1)*X,m*x2+b==x2*log(x2)/n1+(1-x2)*log(1-x2)/n2+x2*(1-x2)*X,m==(1/n1-1/n2)+(log(x1)/n1-log(1-x1)/n2)+(1-2*x1)*X],[x1, x2, m, b],[guess1,guess2,guess3,0]);
        binodal(i,:)=[sol.x1 sol.x2];
    catch
        break
    end
end

ibinodal=i+1;

xbinodal=[binodal(ibinodal:1:nT-1,1);spinodal(nT,2);binodal(nT-1:-1:ibinodal,2)];
ybinodal=[T_range(ibinodal:1:nT) T_range(nT-1:-1:ibinodal)];

figure(2014)

plot(xspinodal,yspinodal,xbinodal,ybinodal)
grid on
xlim([0 1])
xlabel('Concentration, c')
ylabel('Temperature, T')


end

function solution=freeE(c1,n1,n2,X)

c2=1-c1;
solution=c1*log(c1)/n1+c2*log(c2)/n2+c1*c2*X;

end