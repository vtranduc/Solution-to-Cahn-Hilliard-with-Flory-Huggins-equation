function [xspinodal,yspinodal,xbinodal,ybinodal]=...
    phase_diagram(n1,n2,entropy,Ttheta,nT,T_min)

tol=1.0e-6;

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

Xcut1=(-B+sqrt(B^2-4*A*C))/(2*A);=zeros(nT-1,2);

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

% figure(2014)
% 
% plot(xspinodal,yspinodal,xbinodal,ybinodal)
% grid on
% xlim([0 1])
% xlabel('Concentration, c')
% ylabel('Temperature, T')


end

function solution=freeE(c1,n1,n2,X)

c2=1-c1;
solution=c1*log(c1)/n1+c2*log(c2)/n2+c1*c2*X;

end