function frame_cap=illustrate_and_analyze_curved_structure_3d_(c,nx,ny,nz,time,...
	x_minus_x,x_minus_y,x_minus_z,x_plus_x,x_plus_y,x_plus_z,...
    y_minus_x,y_minus_y,y_minus_z,y_plus_x,y_plus_y,y_plus_z,...
    z_minus_x,z_minus_y,z_minus_z,z_plus_x,z_plus_y,z_plus_z)

z_minus_c=zeros(ny,nx);

for yth=1:1:ny
    for xth=1:1:nx
        node=get_node_3d(xth,yth,1,nx,ny);
        z_minus_c(yth,xth)=c(8*(node-1)+1);
    end
end

z_plus_c=zeros(ny,nx);

for yth=1:1:ny
    for xth=1:1:nx
        node=get_node_3d(xth,yth,nz,nx,ny);
        z_plus_c(yth,xth)=c(8*(node-1)+1);
    end
end

y_minus_c=zeros(nz,nx);

for zth=1:1:nz
    for xth=1:1:nx
        node=get_node_3d(xth,1,zth,nx,ny);
        y_minus_c(zth,xth)=c(8*(node-1)+1);
    end
end

y_plus_c=zeros(nz,nx);

for zth=1:1:nz
    for xth=1:1:nx
        node=get_node_3d(xth,ny,zth,nx,ny);
        y_plus_c(zth,xth)=c(8*(node-1)+1);
    end
end

x_minus_c=zeros(nz,ny);

for zth=1:1:nz
    for yth=1:1:ny
        node=get_node_3d(1,yth,zth,nx,ny);
        x_minus_c(zth,yth)=c(8*(node-1)+1);
    end
end

x_plus_c=zeros(nz,ny);

for zth=1:1:nz
    for yth=1:1:ny
        node=get_node_3d(nx,yth,zth,nx,ny);
        x_plus_c(zth,yth)=c(8*(node-1)+1);
    end
end

surf(z_minus_x,z_minus_y,z_minus_z,z_minus_c)
hold on
surf(z_plus_x,z_plus_y,z_plus_z,z_plus_c)
surf(y_minus_x,y_minus_y,y_minus_z,y_minus_c)
surf(y_plus_x,y_plus_y,y_plus_z,y_plus_c)
surf(x_minus_x,x_minus_y,x_minus_z,x_minus_c)
surf(x_plus_x,x_plus_y,x_plus_z,x_plus_c)
hold off

colorbar

str=sprintf('Nonlinear Cahn Hilliard in 3D geometry\n t=%d',time);
suptitle(str)
frame_cap=getframe(gcf);
end