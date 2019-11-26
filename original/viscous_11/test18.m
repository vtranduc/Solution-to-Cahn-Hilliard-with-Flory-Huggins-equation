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
guess_step=0.0001; %MAYBE THERE'S A BETTER WAY TO ADJUST THIS
guess_step_adjust=2;
guess_steps_min=10000; %Must be pretty large number

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
% plot(xspinodal,yspinodal,'.')
% grid on

%THE FOLLOWING WORKS ONLY IF N2>N1!!!!!!
%Compute binodal curve
xbinodal=zeros(1,2*nT-1);
ybinodal=yspinodal;
xbinodal(nT)=xspinodal(nT);

% chi=get_chi(ybinodal(nT-1),entropy,T_theta)
% freePlot(n1,n2,chi)
% min_bound=fminbnd(@(x)freeE(x,n1,n2,chi),0,1)
% g1=(xspinodal(nT-1)+0)/2
% 
% g1=xspinodal(nT-1)
% 
% g2=(xspinodal(nT+1)+1)/2;
% g3=(freeE(xspinodal(nT+1),n1,n2,chi)-freeE(xspinodal(nT-1),n1,n2,chi))/...
% 	(xspinodal(nT+1)-xspinodal(nT-1))
% % g4=freeE(g2,n1,n2,chi)-g2*g3;
% % g3=0
% g4=-g3
% hold on
% plot(g1,freeE(g1,n1,n2,chi),'o')
% plot(g2,freeE(g2,n1,n2,chi),'o')
% hold off
% grid on
% newton_raphson_4_variates(g1,g2,g3,g4,n1,n2,chi);
% 
% g1=xspinodal(nT-1);
% g2=xspinodal(nT+1);
% g3=(freeE(g2,n1,n2,chi)-freeE(g1,n1,n2,chi))/(g2-g1);
% disp('enter func')
[lower,higher]=tangent_finder(xspinodal(nT-1),xspinodal(nT+1),n1,n2,chi)

% for k=1:1:1
%     i=nT-k;
%     j=nT+k;
%     
%     chi=get_chi(ybinodal(i),entropy,T_theta);
%     g1=xspinodal(i);
%     g2=xspinodal(j);
%     
%     min_bound=fminbnd(@(x)freeE(x,n1,n2,chi),0,1);
%     compare=(g1-min_bound)/guess_steps_min;
%     while guess_step>compare
% %         disp('hello peon')
%         guess_step=guess_step/guess_step_adjust;
%     end
%     compare=(1-g2)/guess_steps_min;
%     while guess_step>compare
%         guess_step=guess_step/guess_step_adjust;
%     end
%         
%     for kkk=1:1:5
%         %Work on the left-hand-side
% 
%         %First iteration
%         c2freeE=freeE(g2,n1,n2,chi);
%         [m,b]=local_tangent(g1,n1,n2,chi);
%         error_=abs(m*g2+b-c2freeE);
% 
%         disp('check1')
%         while 1 %Maybe adjust so that the first few steps must be getting smaller?
%             g1=g1-guess_step;
%             guess_step
%             if g1<min_bound
%                 
%                 g1=g1+guess_step;
%                 guess_step=guess_step/guess_step_adjust;
%                 continue
%             end
%             [m,b]=local_tangent(g1,n1,n2,chi);
%             compare=abs(m*g2+b-c2freeE);
% %             if compare==error_
% %                 g1=g1+guess_step;
% %                 guess_step=guess_step*guess_step_adjust;
%             if compare>=error_
%     %             disp('bigger')
%                 g1=g1+guess_step;
%                 guess_step=guess_step/guess_step_adjust;
%             else
%                 break
%             end
%         end
%         %Continue until step is larger
%         disp('check2')
%         while 1
%             g1=g1-guess_step;
%             [m,b]=local_tangent(g1,n1,n2,chi);
%             compare=abs(m*g2+b-c2freeE);
%             if compare>=error_
%                 g1=g1+guess_step;
%                 break
%             end
%         end
%         if g1<min_bound
%             error('g1 is being guessed as minimum bound!')
%         end
% 
%         %Move on to adjust the right-hand-side
% 
%         c1freeE=freeE(g1,n1,n2,chi);
%         [m,b]=local_tangent(g2,n1,n2,chi);
%         error_=abs(m*g1+b-c1freeE);
%         while 1 %Maybe adjust so that the first few steps must be getting smaller?
%             g2=g2+guess_step;
%             if g2>1
%                 g2=g2-guess_step;
%                 guess_step=guess_step/guess_step_adjust;
%                 continue
%             end
%             [m,b]=local_tangent(g2,n1,n2,chi);
%             compare=abs(m*g1+b-c1freeE);
% %             if compare==error_
% %                 g2=g2-guess_step;
% %                 guess_step=guess_step*guess_step_adjust;
%             if compare>=error_
%     %             disp('bigger')
%                 g2=g2-guess_step;
%                 guess_step=guess_step/guess_step_adjust;
%             else
%                 break
%             end
%         end
%         while 1
%             g2=g2+guess_step;
%             if g2>=1
%                 error('g2 thrown out of bound')
%             end
%             [m,b]=local_tangent(g2,n1,n2,chi);
%             compare=abs(m*g1+b-c1freeE);
%             if compare>=error_
%                 g2=g2+guess_step;
%                 break
%             end
%         end
%         [m1,b1]=local_tangent(g1,n1,n2,chi);
%         [m2,b2]=local_tangent(g2,n1,n2,chi);
%         error__=sqrt((m1-m2)^2+(b1-b2)^2)
%     end
%     
%     g3=m2;
%     g4=b2;
%     
%     disp('np method')
% %     for kk=1:1:5
% %         [residual,jacobian]=newton_raphson_frame(g1,g2,g3,g4,n1,n2,chi);
% % 
% %         adjust=jacobian\residual;
% % 
% %         g1=g1+adjust(1);
% %         g2=g2+adjust(2);
% %         if g2>=1
% %             g2=0.999999;
% %         end
% %         g3=g3+adjust(3);
% %         g4=g4+adjust(4);
% % 
% %         error__=sqrt(sum(adjust.^2))
% %     end
%     
% %     while 1
% %         g1=g1-guess_step;
% %         [m,b]=local_tangent(g1,n1,n2,chi);
% %         compare=abs(m*g2+b-c2freeE);
% %         if compare<error_
% %             error_=compare;
% %         else
% %             g1
% % %             error('just stop here')
% %             break
% %         end
% %     end
% end
    
    
%     g3=(freeE(g2,n1,n2,chi)-freeE(g1,n1,n2,chi))/...
%         (g2-g1);
%     g4=freeE(g2,n1,n2,chi)-g3*g2;
%     
%     freePlot(n1,n2,chi)
%     
%     hold on
%     for kk=1:1:1
%         [residual,jacobian]=inflection_method_frame(g1,g2,g3,g4,n1,n2,chi);
% 
%         adjust=jacobian\residual;
%         g1=g1+adjust(1);
%         g3=g3+adjust(2);
%         g4=g4+adjust(3);
%         error_=sqrt(sum(adjust.^2))
%         
%         [m,b]=local_tangent(g1,n1,n2,chi);
%         freeTangentPlot(m,b);
%     end
%     hold off
    


% hold on
% [m,b]=local_tangent(0.2,n1,n2,chi);
% freeTangentPlot(m,b);
% hold off

% hold on
% plot(xbinodal,ybinodal,'o')
% hold off

end

function newton_raphson_4_variates(c1,c2,m,b,n1,n2,chi)
guess=[c1,c2,m,b];
for i=1:1:5
[residual,jacobian]=newton_raphson_frame(...
    guess(1),guess(2),guess(3),guess(4),n1,n2,chi);
adjust=jacobian\residual;
guess=guess+adjust';
error_=sqrt(sum(adjust.^2))
end
end

function [lower_pt,higher_pt]=tangent_finder(...
    infl_left,infl_right,n1,n2,chi)
% guess=[c1,c2,m,b]
tol=1.0d-6;
n=1000;
min_bound_left=fminbnd(@(x)freeE(x,n1,n2,chi),0,infl_left);
trial_pts_left=linspace(min_bound_left,infl_left,n);
min_bound_right=fminbnd(@(x)freeE(x,n1,n2,chi),infl_right,1);
trial_pts_right=linspace(min_bound_right,1-tol,n);
error_=inf;
for i=1:1:n
    for j=1:1:n
        [m1,b1]=local_tangent(trial_pts_left(i),n1,n2,chi);
        [m2,b2]=local_tangent(trial_pts_right(j),n1,n2,chi);
        spacing=area_twoLines(m1,b1,m2,b2);
        if spacing<error_
            ileft=i;
            iright=j;
            if spacing==0
                lower_pt=trial_pts_left(ileft);
                higher_pt=trial_pts_right(iright);
                return
            else
                error_=spacing;
            end
        end
    end
end
%Refine the result========================
accuracy_left=trial_pts_left(2)-trial_pts_left(1);
accuracy_right=trial_pts_right(2)-trial_pts_right(1);
if accuracy_left<=tol && accuracy_right<=tol
    lower_pt=trial_pts_left(ileft);
    higher_pt=trial_pts_right(iright);
    return
end
if accuracy_left>tol
    nleft=round(2*accuracy_left/tol);
else
    nleft=2;
end
if accuracy_right>tol
    nright=round(2*accuracy_right/tol);
else
    nright=2;
end
%=====================================
lower_bound=trial_pts_left(ileft)-accuracy_left;
higher_bound=trial_pts_left(ileft)+accuracy_left;
if lower_bound<min_bound_left
    lower_bound=min_bound_left;
end
if higher_bound>infl_left
    higher_bound=infl_left;
end
trial_pts_left=linspace(lower_bound,higher_bound,nleft);
%======================================
lower_bound=trial_pts_right(iright)-accuracy_right;
higher_bound=trial_pts_right(iright)+accuracy_right;
if lower_bound<infl_right
    lower_bound=infl_right;
end
if higher_bound>=1
    higher_bound=1-tol;
end
trial_pts_right=linspace(lower_bound,higher_bound,nright);
%========================================
error_=inf;
for i=1:1:nleft
    for j=1:1:nright
        [m1,b1]=local_tangent(trial_pts_left(i),n1,n2,chi);
        [m2,b2]=local_tangent(trial_pts_right(j),n1,n2,chi);
        spacing=area_twoLines(m1,b1,m2,b2);
        if spacing<error_
            ileft=i;
            iright=j;
            if spacing==0
                lower_pt=trial_pts_left(ileft);
                higher_pt=trial_pts_right(iright);
                return
            else
                error_=spacing;
            end
        end
    end
end
lower_pt=trial_pts_left(ileft);
higher_pt=trial_pts_right(iright);
end

function solution=area_twoLines(m1,b1,m2,b2)
%Area bounded between x=0 and x=1
if m1==m2
    if b1==b2
        solution=0;
        return
    else
        solution=abs(b1-b2);
        return
    end
end
intersection=(b2-b1)/(m1-m2);
if intersection>=1 || intersection<=0
    a1=0.5*m1+b1;
    a2=0.5*m2+b2;
    solution=abs(a1-a2);
else
    a1=0.5*m1*intersection^2+b1*intersection;
    a2=0.5*m2*intersection^2+b2*intersection;
    a3=0.5*m1+b1-(0.5*m1*intersection^2+b1*intersection);
    a4=0.5*m2+b2-(0.5*m2*intersection^2+b2*intersection);
    solution=abs((a1-a2)+(a4-a3));
end
end

function freeTangentPlot(m,b);
n=1000;
domain=linspace(0,1,n);
range=zeros(1,n);
for i=1:1:n
    range(i)=domain(i)*m+b;
end
plot(domain,range)
end

function [m,b]=local_tangent(c1,n1,n2,chi)
m=dfreeE_dc1(c1,n1,n2,chi);
b=freeE(c1,n1,n2,chi)-m*c1;
end

function [residual,jacobian]=inflection_method_frame(c1,c2,m,b,n1,n2,chi)
jacobian=zeros(3,3);
jacobian(1,1)=d2freeE_dc12(c1,n1,n2,chi);
jacobian(1,2)=-1;
jacobian(2,1)=dfreeE_dc1(c1,n1,n2,chi)-m;
jacobian(2,2)=-c1;
jacobian(2,3)=-1;
jacobian(3,2)=-c2;
jacobian(3,3)=-1;

residual=zeros(3,1);
residual(1)=dfreeE_dc1(c1,n1,n2,chi)-m;
residual(2)=freeE(c1,n1,n2,chi)-m*c1-b;
residual(3)=freeE(c2,n1,n2,chi)-m*c2-b;
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