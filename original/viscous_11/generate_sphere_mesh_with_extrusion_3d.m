function [X,Y,Z]=generate_sphere_mesh_with_extrusion_3d(nx,ny,nz,n,n_eSurf,n_nSurf,ne)
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
box_ratio=(n_eSurf(6)/ne+1)^(-1/3);
X=X*box_ratio;
Y=Y*box_ratio;
Z=Z*box_ratio;
X_=zeros(1,n_nSurf(6));
Y_=zeros(1,n_nSurf(6));
Z_=zeros(1,n_nSurf(6));
yth=0;
zth=1;
for index=1:1:n_nSurf(1)
    yth=yth+1;
    if yth==ny+1
        yth=1;
        zth=zth+1;
    end
    node=get_node_3d(1,yth,zth,nx,ny);
    xyz=extrusion_of_sphere_3d([X(node) Y(node) Z(node)]);
    X_(index)=xyz(1);Y_(index)=xyz(2);Z_(index)=xyz(3);
end
yth=0;
zth=1;
for index=n_nSurf(1)+1:1:n_nSurf(2)
    yth=yth+1;
    if yth==ny+1
        yth=1;
        zth=zth+1;
    end
    node=get_node_3d(nx,yth,zth,nx,ny);
    xyz=extrusion_of_sphere_3d([X(node) Y(node) Z(node)]);
    X_(index)=xyz(1);Y_(index)=xyz(2);Z_(index)=xyz(3);
end
xth=1;
zth=1;
for index=n_nSurf(2)+1:1:n_nSurf(3)
    xth=xth+1;
    if xth==nx
        xth=2;
        zth=zth+1;
    end
    node=get_node_3d(xth,1,zth,nx,ny);
    xyz=extrusion_of_sphere_3d([X(node) Y(node) Z(node)]);
    X_(index)=xyz(1);Y_(index)=xyz(2);Z_(index)=xyz(3);
end
xth=1;
zth=1;
for index=n_nSurf(3)+1:1:n_nSurf(4)
    xth=xth+1;
    if xth==nx
        xth=2;
        zth=zth+1;
    end
    node=get_node_3d(xth,ny,zth,nx,ny);
    xyz=extrusion_of_sphere_3d([X(node) Y(node) Z(node)]);
    X_(index)=xyz(1);Y_(index)=xyz(2);Z_(index)=xyz(3);
end
xth=1;
yth=2;
for index=n_nSurf(4)+1:1:n_nSurf(5)
    xth=xth+1;
    if xth==nx
        xth=2;
        yth=yth+1;
    end
    node=get_node_3d(xth,yth,1,nx,ny);
    xyz=extrusion_of_sphere_3d([X(node) Y(node) Z(node)]);
    X_(index)=xyz(1);Y_(index)=xyz(2);Z_(index)=xyz(3);
end
xth=1;
yth=2;
for index=n_nSurf(5)+1:1:n_nSurf(6)
    xth=xth+1;
    if xth==nx
        xth=2;
        yth=yth+1;
    end
    node=get_node_3d(xth,yth,nz,nx,ny);
    xyz=extrusion_of_sphere_3d([X(node) Y(node) Z(node)]);
    X_(index)=xyz(1);Y_(index)=xyz(2);Z_(index)=xyz(3);
end
X=[X X_];Y=[Y Y_];Z=[Z Z_];
end

function [xyz_]=extrusion_of_sphere_3d(xyz)
len=sqrt(sum(xyz.^2));
xyz_=xyz/len;
end


