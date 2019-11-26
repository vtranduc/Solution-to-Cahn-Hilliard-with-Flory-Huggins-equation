function test53

clear
clc

%------------------------------

nx=20;
ny=20;
nz=20;

nex=nx-1;
ney=ny-1;
nez=nz-1;
ne=nex*ney*nez;
nexney=nex*ney;

n=nx*ny*nz;

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

%------------------------------

xyz_cube=[X' Y' Z'];
xyz_sphere=zeros(n,3);
for i=1:1:n
    xyz_sphere(i,:)=map_to_sphere_with_ratio(xyz_cube(i,:));
end
plot3(xyz_sphere(:,1),xyz_sphere(:,2),xyz_sphere(:,3),'bo')



eTest=1;

[eX,eY,eZ]=get_eXYZ_3d(eTest,nex,nx,ny,nexney,xyz_sphere(:,1),xyz_sphere(:,2),xyz_sphere(:,3));

hold on
for i=1:1:8
    plot3(eX(i),eY(i),eZ(i),'r*')
end
hold off


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