function test90

clear
all; clc;
%=======================================
===================================
nx=500;
ny=500;
dx=0.03; %e x
direction
dy=0.03 ; %grid size y
direction
dt=0.0003 ; %time step
tau=0.0003;
epsilonbar= 0.01;
mu=1.0;
k=4 ;
delta=0.02;
anisotropy= 4.0 ;
alpha= 0.9;
gamma= 10.0;
teq=1.0 ;
x= 1:nx;
y= 1:ny;
t=zeros(nx,ny) ; epsilon
= zeros(nx,ny);
epsilon_derivative =
zeros(nx,ny); phi1= zeros(nx,ny);
grad_phi_x=zeros(nx,ny);
grad_phi_y=zeros(nx,ny);
lap_phi=zeros(nx,ny);
lap_t=zeros(nx,ny);
phinew=zeros(nx,ny);
tnew= zeros(nx,ny);
angl= zeros(nx,ny);
times=0;
%====================setting the
nuclei==================================
==
for i=1:nx
for j=1:ny
phi1(i,j)=0;
if ((i-nx/2)*(i-nx/2)+(jny/2)*(j-ny/2)<20);
phi1(i,j) =1;
end
end
end
phi=phi1;
%=====================gradient&Laplace==
===================================
for times=1:2000000
for i=1:nx
for j=1:ny
    
    %periodic boundary
condition jp = j+1;
jm = j-1;
ip = i+1;
im = i-1;
if (im==0)
im=nx;
elseif (ip==nx+1)
ip=1;
end
if (jm==0)
jm=ny;
elseif (jp==ny+1)
jp=1;
end
%gradient
grad_phi_x(i,j) = (phi(ip,j)
- phi(im,j))/dx;
grad_phi_y(i,j) = (phi(i,jp)
- phi(i,jm))/dy;
%laplacian
lap_phi(i,j) =
(2.0*(phi(ip,j)+phi(im,j)+phi(i,jp)+
phi(i,jm)) +
phi(ip,jp)+phi(im,jm)+phi(im,jp)+
phi(ip,jm) - 12.0*phi(i,j))/(3.0*dx*dx);
lap_t(i,j) =
(2.0*(t(ip,j)+t(im,j)+t(i,jp)+t(i,jm)) +
t(ip,jp)+t(im,jm)+t(im,jp)+t(ip,jm) -
12.0*t(i,j))/(3.0*dx*dx);
if (grad_phi_x(i,j)==0) if
(grad_phi_y(i,j)<0 )
angl(i,j) =-0.5*pi;
elseif (grad_phi_y(i,j)>0 )
angl(i,j) = 0.5*pi;
end
end
if (grad_phi_x(i,j)>0) if
(grad_phi_y(i,j)<0)
angl(i,j)= 2.0*pi +
atan(grad_phi_y(i,j)/grad_phi_x(i,j));
elseif (grad_phi_y(i,j)>0)
angl(i,j)=
atan(grad_phi_y(i,j)/grad_phi_x(i,j));
end
end
if (grad_phi_x(i,j)<0)
angl(i,j) = pi +
atan(grad_phi_y(i,j)/grad_phi_x(i,j));
end
epsilon(i,j) = epsilonbar*(1.0 +
delta*cos(anisotropy*(angl(i,j))));
epsilon_derivative(i,j) = -
epsilonbar*anisotropy*delta*sin(anisotro
py*(angl(i,j)));
grad_epsilon2_x = (epsilon(ip,j)^2
- epsilon(im,j)^2)/dx;
grad_epsilon2_y = (epsilon(i,jp)^2
- epsilon(i,jm)^2)/dy;

end
end
%===================evolution
========================================
=====
for i=1:nx
for j=1:ny
% periodic boundary
condition jp = j+1;
jm = j-1;
ip = i+1;
im = i-1;
if (im==0)
im=nx;
elseif (ip==nx+1)
ip=1;
end
if (jm==0)
jm=ny;
elseif (jp==ny+1)
jp=1;
end
term1=
(epsilon(i,jp)*epsilon_derivative(i,jp)*
grad_phi_x(i,jp) - epsilon(i,jm)*epsilon
_derivative(i,jm)*g rad_phi_x(i,jm))/dy;
term2= -
(epsilon(ip,j)*epsilon_derivative(ip,j)*
grad_phi_y(ip,j) - epsilon(im,j)*epsilon
_derivative(im,j)*g rad_phi_y(im,j))/dx;
term3= grad_epsilon2_x*grad_phi_x(i,j)+
grad_epsilon2_y*grad_phi_y(i,j);

phiold = phi(i,j);
m = alpha/pi * atan(gamma*(teq-t(i,j)));
%time evolution
phinew(i,j) = phi(i,j) + (term1 + term2
+ epsilon(i,j)*epsilon(i,j)*lap_phi(i,j)
+ term3 + phiold*(1.0-
phiold)*(phiold-0.5+m))*dt/tau;
tnew(i,j) =t(i,j) + lap _t(i,j)*dt
+ k*(phi(i,j) - phiold);
phi(i,j)=phinew(i,j);
t(i,j)=tnew(i,j);
%visualization of the
output disp(phi);
figure(1);
image(phi*50);colormap('jet(64)');
pcolor(phi);shading flat; axis ('xy');
end
end
end

end