function [x,y,z,x_coord,y_coord,z_coord]=rectangular_meshing_3D(nex,ney,nez,nx,ny,nz,n,xlen,ylen,zlen)

x=zeros(1,n);y=zeros(1,n);z=zeros(1,n);

x_space=xlen/nex;
for i=1:nx
    x_coord=x_space*(i-1);
    for j=1+ny*(i-1):nx*ny:1+ny*(i-1)+nx*ny*(nz-1)
        x(j:j+ny-1)=x_coord;
    end
end

y_space=ylen/ney;
for i=1:ny
    y_coord=y_space*(i-1);
    for j=1:nz
        for k=i+nx*ny*(j-1):ny:i+nx*ny*(j-1)+ny*(nx-1)
            y(k)=y_coord;
        end
    end
end

z_space=zlen/nez;
for i=1:nz
    z_coord=z_space*(i-1);
    z(nx*ny*(i-1)+1:nx*ny*(i-1)+nx*ny)=z_coord;
end

x_coord=linspace(0,xlen,nx);
y_coord=linspace(0,ylen,ny);
z_coord=linspace(0,zlen,nz);

end