function [X,Y,Z]=test39(nx,ny,nz,n)
X=zeros(1,n);
Y=zeros(1,n);
Z=zeros(1,n);
x_coord=linspace(0,1,nx);
y_coord=linspace(0,1,ny);
z_coord=linspace(0,1,nz);
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
end