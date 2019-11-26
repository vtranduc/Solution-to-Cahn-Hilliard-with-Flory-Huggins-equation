function [x_minus_x,x_minus_y,x_minus_z,x_plus_x,x_plus_y,x_plus_z,...
    y_minus_x,y_minus_y,y_minus_z,y_plus_x,y_plus_y,y_plus_z,...
    z_minus_x,z_minus_y,z_minus_z,z_plus_x,z_plus_y,z_plus_z,...
    x_plot,y_plot,z_plot]=...
    setUp_illustration_curved_structure_3d(nx,ny,nz,X,Y,Z)

z_minus_x=zeros(ny,nx);
z_minus_y=zeros(ny,nx);
z_minus_z=zeros(ny,nx);

for yth=1:1:ny
    for xth=1:1:nx
        node=get_node_3d(xth,yth,1,nx,ny);
        z_minus_x(yth,xth)=X(node);
        z_minus_y(yth,xth)=Y(node);
        z_minus_z(yth,xth)=Z(node);
    end
end

z_plus_x=zeros(ny,nx);
z_plus_y=zeros(ny,nx);
z_plus_z=zeros(ny,nx);

for yth=1:1:ny
    for xth=1:1:nx
        node=get_node_3d(xth,yth,nz,nx,ny);
        z_plus_x(yth,xth)=X(node);
        z_plus_y(yth,xth)=Y(node);
        z_plus_z(yth,xth)=Z(node);
    end
end

y_minus_x=zeros(nz,nx);
y_minus_y=zeros(nz,nx);
y_minus_z=zeros(nz,nx);

for zth=1:1:nz
    for xth=1:1:nx
        node=get_node_3d(xth,1,zth,nx,ny);
        y_minus_x(zth,xth)=X(node);
        y_minus_y(zth,xth)=Y(node);
        y_minus_z(zth,xth)=Z(node);
    end
end

y_plus_x=zeros(nz,nx);
y_plus_y=zeros(nz,nx);
y_plus_z=zeros(nz,nx);

for zth=1:1:nz
    for xth=1:1:nx
        node=get_node_3d(xth,ny,zth,nx,ny);
        y_plus_x(zth,xth)=X(node);
        y_plus_y(zth,xth)=Y(node);
        y_plus_z(zth,xth)=Z(node);
    end
end

x_minus_x=zeros(nz,ny);
x_minus_y=zeros(nz,ny);
x_minus_z=zeros(nz,ny);

for zth=1:1:nz
    for yth=1:1:ny
        node=get_node_3d(1,yth,zth,nx,ny);
        x_minus_x(zth,yth)=X(node);
        x_minus_y(zth,yth)=Y(node);
        x_minus_z(zth,yth)=Z(node);
    end
end

x_plus_x=zeros(nz,ny);
x_plus_y=zeros(nz,ny);
x_plus_z=zeros(nz,ny);

for zth=1:1:nz
    for yth=1:1:ny
        node=get_node_3d(nx,yth,zth,nx,ny);
        x_plus_x(zth,yth)=X(node);
        x_plus_y(zth,yth)=Y(node);
        x_plus_z(zth,yth)=Z(node);
    end
end

x_plot=zeros(ny,nx,nz);
y_plot=zeros(ny,nx,nz);
z_plot=zeros(ny,nx,nz);
inode=0;
for zth=1:1:nz
    for yth=1:1:ny
        for xth=1:1:nx
            inode=inode+1;
            x_plot(yth,xth,zth)=X(inode);
            y_plot(yth,xth,zth)=Y(inode);
            z_plot(yth,xth,zth)=Z(inode);
        end
    end
end

end