function test51

clear
clc

a=[...
    -1 -1 -1
    1 -1 -1
    1 1 -1
    -1 1 -1
    -1 -1 1
    1 -1 1
    1 1 1
    -1 1 1];

% plot3(a(:,1),a(:,2),a(:,3),'bo')
% grid on

b=[...
    0 -1 -1
    1 0 -1
    0 1 -1
    -1 0 -1
    -1 -1 0
    1 -1 0
    1 1 0
    -1 1 0
    0 -1 1
    1 0 1
    0 1 1
    -1 0 1];

% hold on
% plot3(b(:,1),b(:,2),b(:,3),'rx')
% hold off

aqua=[a;b];
d=zeros(20,3);


pt=20

% plot3(aqua(:,1),aqua(:,2),aqua(:,3),'rx',aqua(pt,1),aqua(pt,2),aqua(pt,3),'bo')
% grid on
% xlabel('x');ylabel('y');zlabel('z')
% return

for i=1:1:20
    x=aqua(i,1)*sqrt(1-aqua(i,2)^2/2-aqua(i,3)^2/2+(aqua(i,2)*aqua(i,3))^2/3);
    y=aqua(i,2)*sqrt(1-aqua(i,1)^2/2-aqua(i,3)^2/2+(aqua(i,1)*aqua(i,3))^2/3);
    z=aqua(i,3)*sqrt(1-aqua(i,1)^2/2-aqua(i,2)^2/2+(aqua(i,1)*aqua(i,2))^2/3);
    d(i,:)=[x y z];
end
    
% plot3(d(:,1),d(:,2),d(:,3),'rx')
% grid on

% pt=20
% plot3(d(:,1),d(:,2),d(:,3),'rx',d(pt,1),d(pt,2),d(pt,3),'bo')
% grid on
% xlabel('x');ylabel('y');zlabel('z')
% return

%--------------------------------------------

nx=20;
ny=20;
nz=20;

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

% plot3(X,Y,Z,'b.')

xyz=zeros(n,3);

XYZ=[X' Y' Z'];

for inode=1:1:n
    xyz(inode,:)=FEH20(XYZ(inode,:),d);
end

plot3(xyz(:,1),xyz(:,2),xyz(:,3),'ro')


hold on
plot3(d(:,1),d(:,2),d(:,3),'bo')
hold off

%---------------------

test=zeros(20,3);

for i=1:1:20
    test(i,:)=FEH20(aqua(i,:),d);
end

% 1 2 10 12 13 14

pt=20

hold on
plot3(test(:,1),test(:,2),test(:,3),'rx')
plot3(test(pt,1),test(pt,2),test(pt,3),'bx')
hold off

grid on
xlabel('x'),ylabel('y'),zlabel('z')

fprintf('--SEE SO-------\n')

for i=13
    test(i,:)=FEH20(aqua(i,:),d);
end

%-----------------------------------

eTest=1

[eX,eY,eZ]=get_eXYZ_3d(eTest,nex,nx,ny,nexney,xyz(:,1),xyz(:,2),xyz(:,3));

hold on
for i=1:1:8
    plot3(eX(i),eY(i),eZ(i),'b*')
end
hold off


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