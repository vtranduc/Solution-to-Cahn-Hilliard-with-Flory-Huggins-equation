function [X,Y,Z]=generate_sphere_mesh_3d(nx,ny,nz,n)
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
for i=1:1:n
    x=X(i)*sqrt(1-Y(i)^2/2-Z(i)^2/2+(Y(i)*Z(i))^2/3);
    y=Y(i)*sqrt(1-X(i)^2/2-Z(i)^2/2+(X(i)*Z(i))^2/3);
    z=Z(i)*sqrt(1-X(i)^2/2-Y(i)^2/2+(X(i)*Y(i))^2/3);
    X(i)=x;
    Y(i)=y;
    Z(i)=z;
end
end