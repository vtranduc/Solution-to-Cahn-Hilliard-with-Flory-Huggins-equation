function [X,Y,Z]=generate_sphere_mesh_with_ratio(nx,ny,nz,n)
x_coord=linspace(-1,1,nx);
y_coord=linspace(-1,1,ny);
z_coord=linspace(-1,1,nz);
xyz_sphere=zeros(n,3);
node=0;
for zth=1:1:nz
    for yth=1:1:ny
        for xth=1:1:nx
            node=node+1;
            xyz_sphere(node,1)=x_coord(xth);
            xyz_sphere(node,2)=y_coord(yth);
            xyz_sphere(node,3)=z_coord(zth);
        end
    end
end
for i=1:1:n
    xyz_sphere(i,:)=map_to_sphere_with_ratio(xyz_sphere(i,:));
end
X=xyz_sphere(:,1);
Y=xyz_sphere(:,2);
Z=xyz_sphere(:,3);
X=X';Y=Y';Z=Z';
end

function sol=map_to_sphere_with_ratio(xyz)
ratio=max(abs(xyz));
dist2wall=dist_from_origin(xyz/ratio);
stretch=1/dist2wall;
sol=xyz*stretch;
end

function sol=dist_from_origin(xyz)
sol=sqrt(sum(xyz.^2));
end