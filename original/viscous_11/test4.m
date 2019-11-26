function test4

clear
clc

nex=3;
ney=4;
nez=5;

nx=nex+1;
ny=ney+1;
nz=nez+1;

nxny=nx*ny;

interactions=interactions_rec_3D(28,ny,nz,nxny);

nodal_matrix=zeros(ny,nx,nz);


k=0;

for z_=1:1:nz
    for x_=1:1:nx
        for y_=1:1:ny
            k=k+1;
            nodal_matrix(y_,x_,z_)=k;
        end
    end
end

nodal_matrix;

%=========================================================================

n=nx*ny*nz;

[irows,icols,isj,col_placement,~]=compact_rec_3D(n,nx,ny,nz);


if issorted(isj)
    display('Sorted')
else
    display('unsorted')
end

isj;

find(isj==960.9556234423);

k = 5;
n = 2^k-1;
[x,y,z] = sphere(n);

c = hadamard(2^k);
c
% size(x)
% size(y)
% size(z)
% size(c)
% 
% x
% y
% z
% c

% figure
surf(x,y,z,c);
% % 
colormap([1  1  0; 0  1  1])
axis equal

x=[1 1 -1 -1 1
    1 1 -1 -1 1];
%     0 0 0 0 0];
y=[-1 1 1 -1 -1
    -1 1 1 -1 -1];
%     0 0 0 0 0];
z=[0 0 0 0 0
    1 1 1 1 1];
%     1 1 1 1 1];

surf(x,y,z)
xlabel('x');ylabel('y');zlabel('z');

x_=[1 1
    -1 -1];
    
y_=[-1 1
    -1 1];

z_=ones(2,2);

hold on

surf(x_,y_,z_)
% xlabel('x');ylabel('y');zlabel('z');

% hold on
%  A = [0 0 0];
%  B = [1 0 0];
%  C = [0 1 0];
%  D = [0 0 1];
%  E = [0 1 1];
%  F = [1 0 1];
%  G = [1 1 0];
%  H = [1 1 1];
%  P = [A;B;F;H;G;C;A;D;E;H;F;D;E;C;G;B];
%  plot3(P(:,1),P(:,2),P(:,3))
%  axis equal

% [x,y,z]=peaks;
% surf(x,y,z)
% pause % Rotate the surface
% hold on
% surf(x,y,-z)


end