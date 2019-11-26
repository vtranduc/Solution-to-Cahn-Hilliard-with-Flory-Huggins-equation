function compute_spinodal_binodal
%function [xspinodal,yspinodal]=compute_spinodal_binodal()

% function [xspinodal,yspinodal,xbinodal,ybinodal]=...
%     compute_spinodal_binodal(n1,n2,entropy,Ttheta,nT,T_min)

clear
clc

format long

n1=1; %IT MUST BE 1!
%=========================================
n2=100;
%=========================================
entropy=1; %IT MUST BE 1!
T_theta=1; %IT MUST BE 1!
nT=100;
T_min=0.001; %Must be larger than zero. Should be as close to 0 as possible
nc_eachSide=100;
c_min=0.001; %Must be larger than zero. Should be as close to 0 as possible

%Function-specific literals
tol=1.0e-6;
T_min_default_factor=0.0001; %Has to be less than 1
max_one_approacher_iteration=10;
one_approacher_iteration=-1;

%Compute cut-off or critical temperature
gamma=2*n1*n2;
A=gamma^2;
B=-2*gamma*(n1+n2);
C=(n2-n1)^2;
Xcutoff1=(-B+sqrt(B^2-4*A*C))/(2*A);
% Xcutoff2=(-B-sqrt(B^2-4*A*C))/(2*A);
Tcutoff1=(((Xcutoff1-0.5)/entropy+1)^-1)/T_theta;
% Tcutoff2=(((Xcutoff2-0.5)/entropy+1)^-1)/T_theta

%Define temperature range to be evaluated based on cut-off temperature
if T_min>=Tcutoff1
    warning('T_min specified is larger than or equal to critical temperature. Default value is used instead')
    T_=linspace(T_min_default_factor*Tcutoff1,Tcutoff1,nT);
else
    T_=linspace(T_min,Tcutoff1,nT);
end
   
%Compute spinodal curve
xspinodal=zeros(1,nT*2-1);
yspinodal=[T_ T_(nT-1:-1:1)];
for i=1:1:nT-1
    T=T_(i);
    chi=get_chi(T,entropy,T_theta);
    alpha=2*chi*n1*n2;
    beta=n1-n2-alpha;
    radicant=beta^2-4*alpha*n2;
    if abs(radicant)<=tol && radicant~=0
        radicant=0;
    elseif radicant<0
        error('Error detected in computation of spinodal curve\nTemperature=%d does not yield inflection points\nCritical temperature is %d',T,Tcutoff1)
    end
    xspinodal(i)=(-beta-sqrt(radicant))/(2*alpha);
    xspinodal(nT*2-i)=(-beta+sqrt(radicant))/(2*alpha);
end
chi=get_chi(Tcutoff1,entropy,T_theta);
alpha=2*chi*n1*n2;
beta=n1-n2-alpha;
xspinodal(nT)=-beta/(2*alpha);

% subplot(1,2,1)
plot(xspinodal,yspinodal,'.')
% grid on

%Compute binodal curve
xbinodal=zeros(1,2*nT-1);
ybinodal=yspinodal;
xbinodal(nT)=xspinodal(nT);
%c2 is guessed as midtpoint between inflection point and higher bound
for k=1:1:nT-1
    i=nT-k;
    j=nT+k;
    chi=get_chi(ybinodal(i),entropy,T_theta)
    g1=xspinodal(i)/2
    g2=(xspinodal(j)+1)/2
%     g2=0.99999999
    g3=(freeE(xspinodal(j),n1,n2,chi)-freeE(xspinodal(i),n1,n2,chi))/...
        (xspinodal(j)-xspinodal(i))
%     g4=0;
%     
%     subplot(1,2,2)
%     g2=0.99999999
%     g1=xspinodal(i)
%     g2=xspinodal(j)
%     freePlot(n1,n2,get_chi(ybinodal(i),entropy,T_theta))
%     yspinodal(nT)
%     hold on
%     plot(g1,freeE(g1,n1,n2,chi),'ro')
%     plot(g2,freeE(g2,n1,n2,chi),'ro')
%     plot(xspinodal(i),freeE(xspinodal(i),n1,n2,chi),'rx')
%     plot(xspinodal(j),freeE(xspinodal(j),n1,n2,chi),'rx')
    g4=-g3
%     hold off
%     grid on
%     g2=0.99
    
    
    g1=xspinodal(i);
    g2=xspinodal(j);
    g3=(freeE(xspinodal(j),n1,n2,chi)-freeE(xspinodal(i),n1,n2,chi))/...
        (xspinodal(j)-xspinodal(i));
    g4=freeE(g1,n1,n2,chi)-g3*g1;
    
    nnn=1000;
    xhira=linspace(0.0001,0.9999,nnn);
    yhira=zeros(1,nnn);
    yhiro=yhira;
    for iii=1:1:nnn
        yhira(iii)=g3*xhira(iii)+g4;
        yhiro(iii)=freeE(xhira(iii),n1,n2,chi);
    end
%     hold on
    plot(xhira,yhira,'.',xhira,yhiro,'.')
    hold on
    plot(xspinodal(i),freeE(xspinodal(i),n1,n2,chi),'o')
    hold off
%     hold off
    grid on
    axis([xspinodal(i)-0.01 xspinodal(j)+0.01 freeE(xspinodal(i),n1,n2,chi)-0.01 freeE(xspinodal(j),n1,n2,chi)])
%     error('Just stop nowetfa')
%     error('jdasfsdafsadg')
%     g2=0.99999999
%     g1=0.6
%     error('fasdfas')
    
    error_=inf;
    disp('ENETER')
    while error_>tol
        [residual,jacobian]=...
            newton_raphson_frame(g1,g2,g3,g4,n1,n2,chi);
        adjust=jacobian\-residual;
        if any(imag(adjust)) || any(isnan(adjust)) || any(isinf(adjust))
            warning('Switch method')
            error_=-1;
            break
        end
        compare=sqrt(sum(adjust.^2));
        if compare>=error_
            warning('Not converging')
            error_=-1;
            break
        end
        error_=compare
        g1=g1+adjust(1);
        if g1<0
            error_=-1;
            k
            g1
            error('g1')
            break
        end
        g2=g2+adjust(2);
        if g2>1
            error_=-1;
            g2
            k
            error('g2')
            break
        end
        g3=g3+adjust(3);
        g4=g4+adjust(4);
    end
    if error_~=-1
        xbinodal(i)=g1;
        xbinodal(j)=g2;
    else
        break
    end
end

if k==nT-1
    return
end

k
error('Stop before making c2 approach 1')

%Try to guess c2 value to be closer to 1
approacher=one_approacher((xspinodal(j)+1)/2)
for l=k:1:nT-1
    i=nT-l;
    j=nT+l;
    chi=get_chi(ybinodal(i),entropy,T_theta)
    g1=xspinodal(i)/2
%     g2=(xspinodal(j)+1)/2
%     g2=0.99999999
    g3=(freeE(xspinodal(j),n1,n2,chi)-freeE(xspinodal(i),n1,n2,chi))/...
        (xspinodal(j)-xspinodal(i))
    g4=-g3;
    
    while 1
        g2=approacher;
        error_=inf;
        while error_>tol
            [residual,jacobian]=...
                newton_raphson_frame(g1,g2,g3,g4,n1,n2,chi);
            adjust=jacobian\-residual;
            if any(imag(adjust)) || any(isnan(adjust)) || any(isinf(adjust))
                warning('Switch method')
                error_=-1;
                break
            end
            compare=sqrt(sum(adjust.^2));
            if compare>=error_
                warning('Not converging')
                error_=-1;
                break
            end
            error_=compare
            g1=g1+adjust(1);
            g2=g2+adjust(2);
            g3=g3+adjust(3);
            g4=g4+adjust(4);
        end
        if error_~=-1
            xbinodal(i)=g1;
            xbinodal(j)=g2;
            break
        else
            approacher=one_approacher(approacher)
            if approacher==1
                warning('switch method further')
                error_=-2;
                break
            end
        end
    end
    if error_==-2
        break
    end
end
l
if l==nT-1
    return
end
%Bruteforcibly assume c2=1
%YOU SHOULD GUESS G1 AS THE THE LOCAL MINUMUM 
for p=l:1:nT-1
    i=nT-p;
    j=nT+p;
    chi=get_chi(ybinodal(i),entropy,T_theta);
    %[residual,jacobian]=newton_raphson_bruteForce(c1,m,b,n1,n2,chi)
    g1=xspinodal(i);
    g1=fminbnd(@(x)freeE(x,n1,n2,chi),c_min,g1)
    %===========================
%     g1=0.05
%     freePlot(n1,n2,chi)
%     hold on
%     plot(g1,freeE(g1,n1,n2,chi),'x')
%     
% %     plot(xspinodal(i),freeE(xspinodal(i),n1,n2,chi),'ro')
%     
%     hold off
%     
%     test1=fminbnd(@(xxx)freeE(xxx,n1,n2,chi),0.01,0.99)
%     
%     hold on
%     plot(test1,freeE(test1,n1,n2,chi),'bo')
%     hold off
%     
%     error('Jsfdaas')
%     a=n1-n2;
%     b=2*n2;
%     c=-n2;
%     g1=(-b+sqrt(b^2-4*a*c))/(2*a)
%     hold on
%     plot(g1,freeE(g1,n1,n2,chi),'yo')
%     hold off
%     domain=linspace(0.01,0.9,100);
%     range=zeros(1,100);
%     for iii=1:1:100
%         range(iii)=d2freeE_dc12(domain(iii),n1,n2,chi);
%     end
%     hold on
%     plot(domain,range,'gx')
%     plot(xspinodal(i),freeE(xspinodal(i),n1,n2,chi),'ro')
%     axis([0 1 -0.5 0.5])
%     hold off
% %     den=d2freeE_dc12(g1,n1,n2,chi)
%     if d2freeE_dc12(g1,n1,n2,chi)==0 %Ideally, this is imppossible
%         error('Inflection point is too close to guessed value')
%     end
%     error_=inf;
% %     while error_>tol
%     for i=1:1:5
%         disp('loop')
%         adjust=dfreeE_dc1(g1,n1,n2,chi)/d2freeE_dc12(g1,n1,n2,chi);
%         error_=abs(adjust)
%         g1=g1-adjust;
% %         break
%     end
%     
%     hold on
%     plot(g1,freeE(g1,n1,n2,chi),'x')
%     hold off
%     grid on
% %     p
%     error('Stop at Newton Raphson')
    %===========================
%     g2=(xspinodal(j)+1)/2
%     g2=0.99999999
    g2=(freeE(xspinodal(j),n1,n2,chi)-freeE(g1,n1,n2,chi))/...
        (xspinodal(j)-g1)
    
%     g2=0.01
%     g1=xspinodal(i+1)
%     freeE(xspinodal(j),n1,n2,chi)
%     freeE(xspinodal(i),n1,n2,chi)
%     g2=0.3
%     g1=0.05
%     freePlot(n1,n2,chi)
%     grid on
%     hold on
%     plot(g1,freeE(g1,n1,n2,chi),'x')
%     hold off
%     plot(xspinodal(i),freeE(xspinodal(i),n1,n2,chi),'o',...
%         xspinodal(j),freeE(xspinodal(j),n1,n2,chi),'o')
%     hold off
%     return
    g3=-g2
    
    error_=inf;
%     error('just stop')
    while error_>tol
        [residual,jacobian]=newton_raphson_bruteForce(g1,g2,g3,n1,n2,chi);
%         g1
%         jacobian
        adjust=jacobian\-residual;
        if any(imag(adjust)) || any(isnan(adjust)) || any(isinf(adjust))
            p
            adjust
            g1
            g2
            g3
            warning('Just stop here man')
            error_=-1;
            break
        end
        compare=sqrt(sum(adjust.^2));
        if compare>error_
            warning('Not converging')
            error_=-1;
            break
        else
            error_=compare;
        end
        
        g1=g1+adjust(1);
        g2=g2+adjust(2);
        g3=g3+adjust(3);

    end
    if error_==-1
        break
    end
end

% hold on
% domain=linspace(0.0001,0.999,1000);
% range=zeros(1,1000);
% for i=1:1:1000
%     range(i)=freeE(domain(i),n1,n2,chi);
% end
% plot(domain,range)
% hold off

hold on
plot(xbinodal,ybinodal,'x')
hold off
return

% % figure(102)
% % plot(xspinodal,yspinodal,'*')
% % grid on
% 
% %Recompute binodal curve
% 
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
%     T=T_(i);
%     chi=get_chi(T,entropy,T_theta);
%     guess(1)=xspinodal(i)/2;
%     guess(2)=(xspinodal(j)+1)/2;
%     guess(3)=(freeE(xspinodal(j),n1,n2,chi)-freeE(xspinodal(i),n1,n2,chi))/...
%         (xspinodal(j)-xspinodal(i));
% %     guess(3)=0;
% 
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
% disp('Done temperature method')
% if k==1
%     warning('Temperaure method was unsuccessful')
%     xbinodal(1:1:nT)=linspace(c_min,xspinodal(nT),nT);
%     guess(1)=ybinodal(nT);
%     guess(2)=(xbinodal(nT)+1)/2;
% %     spacer=xbinodal(2)-xbinodal(1);
% %     guess(2)=xbinodal(nT)+spacer
%     guess(3)=0;
%     c1=xbinodal(nT-1);
%     syms T
%     chi_symbolic=get_chi(T,entropy,T_theta);
%     sol=vpasolve([dfreeE_dc1(c1,n1,n2,chi_symbolic)==dfreeE_dc1(...
%         c2,n1,n2,chi_symbolic),m==dfreeE_dc1(c1,n1,n2,chi_symbolic),...
%         m*c1+b==freeE(c1,n1,n2,chi_symbolic),m*c2+b==freeE(...
%         c2,n1,n2,chi_symbolic)],[T c2 m b],guess);
%     while length(sol.T)~=1 || length(sol.c2)~=1 || length(sol.m)~=1 ||...
%             length(sol.b)~=1 || sol.T>ybinodal(i+1) ||...
%             sol.c2<xbinodal(j-1) || sol.m<0 || isreal(sol.T)==0 ||...
%             isreal(sol.c2)==0 || isreal(sol.m)==0 || isreal(sol.b)==0
%         
%         one_approacher_iteration=one_approacher_iteration+1
%         if one_approacher_iteration>max_one_approacher_iteration
%             error('Phase diagram is too skewed. vpasolve is unable to solve binodal curve. Increase nT and/or decrease n2?')
%         end
%         
%         guess(2)=one_approacher(guess(2));
%         sol=vpasolve([dfreeE_dc1(c1,n1,n2,chi_symbolic)==dfreeE_dc1(...
%         c2,n1,n2,chi_symbolic),m==dfreeE_dc1(c1,n1,n2,chi_symbolic),...
%         m*c1+b==freeE(c1,n1,n2,chi_symbolic),m*c2+b==freeE(...
%         c2,n1,n2,chi_symbolic)],[T c2 m b],guess);
%         
%         warning('Unable to compute binodal curve. Increase nT and/or decrease n2?')
%     end
% 	xbinodal(nT+1)=sol.c2;
% 	ybinodal(nT-1)=sol.T;
% 	ybinodal(nT+1)=sol.T;
% 	k=2;
% %     end
% elseif k==nT-1
%     return
% else
%     c1_remaining=nT-k;
%     xbinodal(1:1:c1_remaining+1)=...
%         linspace(c_min,xbinodal(c1_remaining+1),c1_remaining+1);
%     syms T
%     chi_symbolic=get_chi(T,entropy,T_theta);
% end
% % error('HOLD')
% disp('Initializing c1 method')
% spacer=xbinodal(2)-xbinodal(1);
% 
% for l=k:1:nT-1
%     i=nT-l;
%     j=nT+l;
%     if xbinodal(j-1)==1
%         error('c2 is getting too close to 1. vpasolve can no longer handle operation. n2 is too high')
%     end
%     c1=xbinodal(i);
%     guess(1)=-spacer...
%         *(ybinodal(i+2)-ybinodal(i+1))/(xbinodal(i+2)-xbinodal(i+1))...
%         +ybinodal(i+1);
%     
% %=========== Should I go trough all the computations for the slope? =======
% %     chi=get_chi(guess(1),entropy,T_theta); %Based on guessed value of T
% %     alpha=2*chi*n1*n2;
% %     beta=n1-n2-alpha;
% %     radicant=beta^2-4*alpha*n2;
% %     disp('checking 123')
% %     [l radicant chi guess(1) ybinodal(i+1) ybinodal(i+2) xbinodal(i+1) xbinodal(i+2)]
% %     disp('end checking')
% %     if abs(radicant)<=tol
% %         guess(3)=0;
% %     elseif radicant<0
% %         radicant
% % %         guess(3)=0;
% %         error('Error detected in computation of binodal curve\nGuessed temperature=%d does not yield inflection points\nCritical temperature is %d',guess(1),Tcutoff1)
% %     else
% %         xinflection1=(-beta-sqrt(radicant))/(2*alpha);
% %         xinflection2=(-beta+sqrt(radicant))/(2*alpha);
% %         guess(3)=(freeE(xinflection2,n1,n2,chi)-freeE(xinflection1,n1,n2,chi))...
% %              /(xinflection2-xinflection1);
% %     end
%     
%     guess(3)=0; %Bypass all
% %============== Fin ====================================================
% %     xbinodal(j-1)
%     guess(2)=(xbinodal(j-1)+1)/2;
% %     guess(2)=0.9999999999999
% 
% 
%     sol=vpasolve([dfreeE_dc1(c1,n1,n2,chi_symbolic)==dfreeE_dc1(c2,n1,n2,chi_symbolic),m==dfreeE_dc1(c1,n1,n2,chi_symbolic),m*c1+b==freeE(c1,n1,n2,chi_symbolic),m*c2+b==freeE(c2,n1,n2,chi_symbolic)],[T c2 m b],guess);
% 	[i l sol.T sol.c2 sol.m sol.b nT-1]
%     while length(sol.T)~=1 || length(sol.c2)~=1 || length(sol.m)~=1 ||...
%             length(sol.b)~=1 || sol.T>ybinodal(i+1) ||...
%             sol.c2<xbinodal(j-1) || sol.m<0 || isreal(sol.T)==0 ||...
%             isreal(sol.c2)==0 || isreal(sol.m)==0 || isreal(sol.b)==0
%         
%         sol.T
%         sol.c2
%         
%         break
%         disp('Recomputing...')
%         guess(2)=one_approacher(guess(2))
%         disp('Solving...')
%         sol=vpasolve([dfreeE_dc1(c1,n1,n2,chi_symbolic)==dfreeE_dc1(c2,n1,n2,chi_symbolic),m==dfreeE_dc1(c1,n1,n2,chi_symbolic),m*c1+b==freeE(c1,n1,n2,chi_symbolic),m*c2+b==freeE(c2,n1,n2,chi_symbolic)],[T c2 m b],guess);
%         disp('solved')
%         [i l sol.T sol.c2 sol.m sol.b nT-1]
% 
%         
%     end
%     
% %     try
% %         sol=vpasolve([dfreeE_dc1(c1,n1,n2,chi_symbolic)==dfreeE_dc1(c2,n1,n2,chi_symbolic),m==dfreeE_dc1(c1,n1,n2,chi_symbolic),m*c1+b==freeE(c1,n1,n2,chi_symbolic),m*c2+b==freeE(c2,n1,n2,chi_symbolic)],[T c2 m b],guess);
% %         [i l sol.T sol.c2 sol.m sol.b nT-1]
% %         %c2 is expected to be never below c1
%     try
% 	ybinodal(i)=sol.T;
% 	ybinodal(j)=sol.T;
% 	xbinodal(j)=sol.c2;
%     catch
%         disp('bad')
%         break
%     end
% %     catch
% %         warning('sup')
% %         sol.T
% %         sol.c2
% %         sol.m
% %         sol.b
% %         break
% %     end
% end
% 
% 
% 
% plot(xspinodal,yspinodal,xbinodal,ybinodal,'r*',xspinodal(nT),yspinodal(nT),'bo',xbinodal(nT-1),ybinodal(nT-1),'go',xbinodal(nT+1),ybinodal(nT+1),'go')
% grid on
% % hold on
% % plot(xbinodal,ybinodal,'*')
% % hold off
% % test1=xbinodal(nT-1:1:nT+20)
% % test2=xbinodal(nT+1:-1:nT-20)
% % 
% % [test1' test2']

end

function solution=freeE(c1,n1,n2,chi)
c2=1-c1;
solution=c1*log(c1)/n1+c2*log(c2)/n2+c1*c2*chi;
end

function solution=dfreeE_dc1(c1,n1,n2,chi)
solution=(1/n1-1/n2)+log(c1)/n1-log(1-c1)/n2+(1-2*c1)*chi;
end

function solution=d2freeE_dc12(c1,n1,n2,chi)
solution=1/(n1*c1)+1/(n2*(1-c1))-2*chi;
end

function solution=one_approacher(closing_val)
%closing_val must be between and not equal to 0 and 1
solution=0;
comparison=1-closing_val;
decimal_place=1;
while 1-solution>=comparison
    decimal_place=decimal_place*0.1;
    solution=solution+decimal_place*9;
    if solution==1
        break
    end
end
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

function [residual,jacobian]=newton_raphson_bruteForce(c1,m,b,n1,n2,chi)
jacobian=zeros(3,3);
jacobian(1,2)=1;
jacobian(1,3)=1;
jacobian(2,1)=d2freeE_dc12(c1,n1,n2,chi);
jacobian(2,2)=-1;
jacobian(3,1)=dfreeE_dc1(c1,n1,n2,chi)-m;
jacobian(3,2)=-c1;
jacobian(3,3)=-1;

residual=zeros(3,1);
residual(1)=m+b;
residual(2)=dfreeE_dc1(c1,n1,n2,chi)-m;
residual(3)=freeE(c1,n1,n2,chi)-m*c1-b;
end

function solution=freeE_amplify(amplifier,c1,n1,n2,chi)
solution=amplifier*freeE(c1,n1,n2,chi);
end

function freePlot(n1,n2,chi)
n=1000;
domain=linspace(0.0001,0.9999,n);
range=zeros(1,n);
for i=1:1:n
    range(i)=freeE(domain(i),n1,n2,chi);
end
plot(domain,range)
end