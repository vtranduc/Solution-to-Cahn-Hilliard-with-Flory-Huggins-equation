function test50(nx,ny,nz,n)

clear
clc



nx=13;
ny=12;
nz=11;

nex=nx-1;
ney=ny-1;
nez=nz-1;
ne=nex*ney*nez;
nexney=nex*ney;

n=nx*ny*nz;

ne
eTest=610

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

[eX,eY,eZ]=get_eXYZ_3d(eTest,nex,nx,ny,nexney,X,Y,Z);
[eX',eY',eZ'];
subplot(1,2,1)
plot3(X,Y,Z,'r.',eX,eY,eZ,'go')

for i=1:1:n
    x=X(i)*sqrt(1-Y(i)^2/2-Z(i)^2/2+(Y(i)*Z(i))^2/3);
    y=Y(i)*sqrt(1-X(i)^2/2-Z(i)^2/2+(X(i)*Z(i))^2/3);
    z=Z(i)*sqrt(1-X(i)^2/2-Y(i)^2/2+(X(i)*Y(i))^2/3);
    X(i)=x;
    Y(i)=y;
    Z(i)=z;
end

[eX,eY,eZ]=get_eXYZ_3d(eTest,nex,nx,ny,nexney,X,Y,Z);
[eX',eY',eZ'];
subplot(1,2,2)
plot3(X,Y,Z,'r.',eX,eY,eZ,'go')

end