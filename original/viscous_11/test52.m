function test52

clear
clc

%--------------------
nx=10;
ny=10;
nz=10;

n=nx*ny*nz;
nex=nx-1;
ney=ny-1;
nexney=nex*ney;

x_coord=linspace(-1,1,nx);
y_coord=linspace(-1,1,ny);
z_coord=linspace(-1,1,nz);
X=zeros(1,n);
Y=zeros(1,n);
Z=zeros(1,n);
node=0;
for zth=1:1:nz
    for yth=1:1:ny
        for xth=1:1:nx
            node=node+1;
            X(node)=x_coord(xth);
            Y(node)=y_coord(yth);
            Z(node)=z_coord(zth);
        end
    end
end
%------------------------------------

b=[...
    0 0 0
    1 0 0
    1 1 0
    0 1 0
    0 0 1
    1 0 1
    1 1 1
    0 1 1
    0.5 0 0 %9
    1 0.5 0
    0.5 1 0
    0 0.5 0 %12
    0 0 0.5 %13
    1 0 0.5
    1 1 0.5 %--15
    0 1 0.5
    0.5 0 1
    1 0.5 1 %--18
    0.5 1 1 %--19
    0 0.5 1];

b=map_to_sphere_3d(b);

plot3(b(:,1),b(:,2),b(:,3),'rx')
grid on
xlabel('x');ylabel('y');zlabel('z')

XYZ=[X' Y' Z'];
xyz=zeros(n,3);
for inode=1:1:n
    xyz(inode,:)=FEH20(XYZ(inode,:),b);
end
plot3(xyz(:,1),xyz(:,2),xyz(:,3),'ro')
grid on
end

function sol=map_to_sphere_3d(coordinates)
[n,~]=size(coordinates);
sol=zeros(n,3);
for i=1:1:n
    x=coordinates(i,1)*sqrt(1-coordinates(i,2)^2/2-coordinates(i,3)^2/2+(coordinates(i,2)*coordinates(i,3))^2/3);
    y=coordinates(i,2)*sqrt(1-coordinates(i,1)^2/2-coordinates(i,3)^2/2+(coordinates(i,1)*coordinates(i,3))^2/3);
    z=coordinates(i,3)*sqrt(1-coordinates(i,1)^2/2-coordinates(i,2)^2/2+(coordinates(i,1)*coordinates(i,2))^2/3);
    sol(i,:)=[x y z];
end
end

function xyz=FEH20(abg,XYZ)

%Each input and output is a vector

xyz=...
    XYZ(1,:)*(1-abg(1))*(1-abg(2))*(1-abg(3))*(-abg(1)-abg(2)-abg(3)-2)/8+...
    XYZ(2,:)*(1+abg(1))*(1-abg(2))*(1-abg(3))*(abg(1)-abg(2)-abg(3)-2)/8+...
    XYZ(3,:)*(1+abg(1))*(1+abg(2))*(1-abg(3))*(abg(1)+abg(2)-abg(3)-2)/8+...
    XYZ(4,:)*(1-abg(1))*(1+abg(2))*(1-abg(3))*(-abg(1)+abg(2)-abg(3)-2)/8+...
    XYZ(5,:)*(1-abg(1))*(1-abg(2))*(1+abg(3))*(-abg(1)-abg(2)+abg(3)-2)/8+...
    XYZ(6,:)*(1+abg(1))*(1-abg(2))*(1+abg(3))*(abg(1)-abg(2)+abg(3)-2)/8+...
    XYZ(7,:)*(1+abg(1))*(1+abg(2))*(1+abg(3))*(abg(1)+abg(2)+abg(3)-2)/8+...
    XYZ(8,:)*(1-abg(1))*(1+abg(2))*(1+abg(3))*(-abg(1)+abg(2)+abg(3)-2)/8+...
    XYZ(9,:)*(1-abg(1)^2)*(1-abg(2))*(1-abg(3))/4+...
    XYZ(10,:)*(1+abg(1))*(1-abg(2)^2)*(1-abg(3))/4+...
    XYZ(11,:)*(1-abg(1)^2)*(1+abg(2))*(1-abg(3))/4+...
    XYZ(12,:)*(1-abg(1))*(1-abg(2)^2)*(1-abg(3))/4+...
    XYZ(13,:)*(1-abg(1))*(1-abg(2))*(1-abg(3)^2)/4+...
    XYZ(14,:)*(1+abg(1))*(1-abg(2))*(1-abg(3)^2)/4+...
    XYZ(15,:)*(1+abg(1))*(1+abg(2))*(1-abg(3)^2)/4+...
    XYZ(16,:)*(1-abg(1))*(1+abg(2))*(1-abg(3)^2)/4+...
    XYZ(17,:)*(1-abg(1)^2)*(1-abg(2))*(1+abg(3))/4+...
    XYZ(18,:)*(1+abg(1))*(1-abg(2)^2)*(1+abg(3))/4+...
    XYZ(19,:)*(1-abg(1)^2)*(1+abg(2))*(1+abg(3))/4+...
    XYZ(20,:)*(1-abg(1))*(1-abg(2)^2)*(1+abg(3))/4;
end