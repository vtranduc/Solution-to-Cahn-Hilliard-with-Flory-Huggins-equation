function test19
%function [xspinodal,yspinodal]=compute_spinodal_binodal()

% function [xspinodal,yspinodal,xbinodal,ybinodal]=...
%     compute_spinodal_binodal(n1,n2,entropy,Ttheta,nT,T_min)

clear
clc
tic
format long

n1=1; %IT MUST BE 1!
%=========================================
n2=10000; %IT MUST BE 1 OR LARGER!!!!
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
parallel_computing=1; %0=no, 1=yes
nInitialGuessPts=1000; %Number of points on each side at which binodal points are guessed

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

%THE FOLLOWING WORKS ONLY IF N2>N1!!!!!!
%Compute binodal curve
xbinodal=zeros(1,2*nT-1);
ybinodal=yspinodal;
xbinodal(nT)=xspinodal(nT);

if parallel_computing==0
    for k=1:1:nT-1
        chi=get_chi(ybinodal(nT-k),entropy,T_theta);
        [xbinodal(nT-k),xbinodal(nT+k)]=tangent_finder(...
            xspinodal(nT-k),xspinodal(nT+k),n1,n2,chi,tol,nInitialGuessPts);
    end
elseif parallel_computing==1
    lows=zeros(1,nT-1);
    highs=zeros(1,nT-1);
    for k=1:1:nT-1
        lows(k)=xspinodal(nT-k);
        highs(k)=xspinodal(nT+k);
    end
    parfor k=1:1:nT-1
        disp('heaf')
        k
        chi=get_chi(yspinodal(nT-k),entropy,T_theta);
        [lows(k),highs(k)]=tangent_finder(...
            lows(k),highs(k),n1,n2,chi,tol,nInitialGuessPts);
    end
    for k=1:1:nT-1
        xbinodal(nT-k)=lows(k);
        xbinodal(nT+k)=highs(k);
    end
else
    error('parallel_computing must be 0=no or 1=yes!')
end

plot(xspinodal,yspinodal,'.',xbinodal,ybinodal,'.')
grid on

% error('Stop before refinement')

%Refinement
disp('=============Refine=================================')
for k=1:1:nT-1
    chi=get_chi(ybinodal(nT-k),entropy,T_theta);
    c1=xbinodal(nT-k);
    c2=xbinodal(nT+k);
    [c1Refined,c2Refined]=refine_binodal(c1,c2,n1,n2,chi,tol);
    xbinodal(nT-k)=c1Refined;
    xbinodal(nT+k)=c2Refined;
end

hold on
plot(xbinodal,ybinodal,'go')
hold off

toc
end

function [c1Refined,c2Refined]=refine_binodal(c1,c2,n1,n2,chi,tol)
%c2=1, it will simply assume c2=1. Otherwise, it will attempt to refine c2
if c2==1
    c1Refined=refine_c1_binodal(c1,n1,n2,chi,tol);
    c2Refined=1;
    return
end
[m,b]=linePassingTwoPoints(c1,freeE(c1,n1,n2,chi),c2,freeE(c2,n1,n2,chi));
guess=[c1,c2,m,b];
%Start iteration
error_=inf;
while error_>tol
    [residual,jacobian]=newton_raphson_frame(...
        guess(1),guess(2),guess(3),guess(4),n1,n2,chi);
    adjust=jacobian\-residual;
    guess=guess+adjust';
    if guess(2)>=1 || guess(1)<=0
        c1Refined=refine_c1_binodal(c1,n1,n2,chi,tol);
        c2Refined=1;
        return
%         error('bad')
    else
        error__=sqrt(sum(adjust.^2));
        if error__>=error_
            c1Refined=refine_c1_binodal(c1,n1,n2,chi,tol);
            c2Refined=1;
            return
%             error('Not converging')
        else
            error_=error__;
        end
    end
end
c1Refined=guess(1);
c2Refined=guess(2);
end

function solution=refine_c1_binodal(c1,n1,n2,chi,tol)
%Assume c2=0
guess1=c1;
guess2=freeE(c1,n1,n2,chi)/(c1-1);
error_=inf;
while error_>tol
    [residual,jacobian]=NR_c2_equal_1_assumed(guess1,guess2,n1,n2,chi);
    adjust=jacobian\-residual;
    guess1=guess1+adjust(1);
    if guess1<=0
        solution=0;
        return
    end
    guess2=guess2+adjust(2);
    error__=sqrt(sum(adjust.^2));
    if error__>=error_
        solution=0;
        return
    else
        error_=error__;
    end
end
disp('SUCCESS!')
solution=guess1;
end

function [residual,jacobian]=NR_c2_equal_1_assumed(c1,m,n1,n2,chi)
jacobian=zeros(2,2);
jacobian(1,1)=d2freeE_dc12(c1,n1,n2,chi);
jacobian(1,2)=-1;
jacobian(2,1)=dfreeE_dc1(c1,n1,n2,chi)-m;
jacobian(2,2)=1-c1;
residual=zeros(2,1);
residual(1)=dfreeE_dc1(c1,n1,n2,chi)-m;
residual(2)=freeE(c1,n1,n2,chi)-m*(c1-1);
end

function [m,b]=linePassingTwoPoints(x1,y1,x2,y2)
m=(y2-y1)/(x2-x1);
b=y1-m*x1;
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
    infl_left,infl_right,n1,n2,chi,tol,n)
%====================================
if dfreeE_dc1(1-tol,n1,n2,chi)<=0 %Local minima is too close to zero
%     warning('what')
%     higher_pt=0;
%     lower_pt=0;
%     return
    higher_pt=1;
    lower_pt=tangent_finder_c1_only(infl_left,higher_pt,n1,n2,chi,tol,n)
    return
end
%========================================
% n is the number of points on each side to be considered initially
min_bound_left=fminbnd(@(x)freeE(x,n1,n2,chi),0,infl_left);
%Note that this will never be 0, because freeE is undefined there
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
if iright==n %The best tangent is very likely higher than 1-tol
    higher_pt=1-tol;
    lower_pt=tangent_finder_c1_only(infl_left,higher_pt,n1,n2,chi,tol,n);
    return
elseif iright==1 || ileft==1 || iright==n
    error('Perhaps the spacing is too huge. Try to increase n')
end
%Qualitatively evaluate the result========
% % sup=realmin('single')
% % tester=1-sup
% % if tester==1
% %     error('no way')
% % end
% % error('Just stop for testing')
% error_
% if iright==n
% %     trial_pts_right(iright)
% %     trial_pts_right'
% %     1-tol
% %     if trial_pts_right(iright)==1-tol
% %         error('numerical error')
% %     end
%     if ileft==1 || ileft==n
%         ileft
%         [m1,b1]=local_tangent(trial_pts_left(ileft),n1,n2,chi)
%         [m2,b2]=local_tangent(trial_pts_right(iright),n1,n2,chi)
%         spacing=area_twoLines(m1,b1,m2,b2)
%         error_
%         infl_right
%         if infl_right>1-tol
%             error('error detected')
%         end
%         %==============================
%         tester=1-exp(n2*((1/n1-1/n2)-chi))
%         if tester==1
%             error('Pretty bad')
%         end
%         %==============================
%         error('bad again==========')
%     end
% end
% infl_right
% tester=1-exp(n2*((1/n1-1/n2)-chi))
% if tester==1
%     error('Pretty good')
% end
% if ileft==1
%     iright
%     error('CASE IMPOSSIBLE!')
% end
%Refine the result========================
accuracy_left=trial_pts_left(2)-trial_pts_left(1);
accuracy_right=trial_pts_right(2)-trial_pts_right(1);
if accuracy_left<=tol && accuracy_right<=tol
    lower_pt=trial_pts_left(ileft);
    higher_pt=trial_pts_right(iright);
    return
end
if accuracy_left>tol
    nleft=ceil(2*accuracy_left/tol);
else
    nleft=2;
end
if accuracy_right>tol
    nright=ceil(2*accuracy_right/tol);
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

function solution=tangent_finder_c1_only(infl_left,c2,n1,n2,chi,tol,n)
x2=c2;
if x2==1
    y2=0;
else
    y2=freeE(x2,n1,n2,chi);
end
min_bound_left=fminbnd(@(x)freeE(x,n1,n2,chi),0,infl_left);
trial_pts_left=linspace(min_bound_left,infl_left,n);
error_=inf;
for i=1:1:n
    x1=trial_pts_left(i);
    y1=freeE(x1,n1,n2,chi);
    [m1,b1]=linePassingTwoPoints(x1,y1,x2,y2);
    [m2,b2]=local_tangent(x1,n1,n2,chi);
    spacing=area_twoLines(m1,b1,m2,b2);
    if spacing<error_
        if spacing==0
            solution=x1;
            return
        else
            ileft=i;
            error_=spacing;
        end
    end
end
solution=trial_pts_left(ileft);
%Refine further==================
accuracy=trial_pts_left(2)-trial_pts_left(1);
if accuracy<=tol
    return
else
    nleft=ceil(2*accuracy/tol);
    lower_bound=solution-accuracy;
    if lower_bound<min_bound_left
        lower_bound=min_bound_left;
    end
    higher_bound=solution+accuracy;
    if higher_bound>infl_left
        higher_bound=infl_left;
    end
    trial_pts_left=linspace(lower_bound,higher_bound,nleft);
end
%----------------------------------------------
error_=inf;
for i=1:1:nleft
    x1=trial_pts_left(i);
    y1=freeE(x1,n1,n2,chi);
    [m1,b1]=linePassingTwoPoints(x1,y1,x2,y2);
    [m2,b2]=local_tangent(x1,n1,n2,chi);
    spacing=area_twoLines(m1,b1,m2,b2);
    if spacing<error_
        if spacing==0
            solution=x1;
            return
        else
            ileft=i;
            error_=spacing;
        end
    end
end
%----------------------------------------------
solution=trial_pts_left(ileft);

% min_bound_left=fminbnd(@(x)freeE(x,n1,n2,chi),0,infl_left);
% trial_pts_left=linspace(min_bound_left,infl_left,n);
% [m2,b2]=local_tangent(c2,n1,n2,chi);
% error_=inf;
% for i=1:1:n
%     [m1,b1]=local_tangent(trial_pts_left(i),n1,n2,chi);
%     spacing=area_twoLines(m1,b1,m2,b2);
%     if spacing<error_
%         ileft=i;
%         if spacing==0
%             solution=trial_pts_left(ileft);
%             return
%         else
%             error_=spacing;
%         end
%     end
% end
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