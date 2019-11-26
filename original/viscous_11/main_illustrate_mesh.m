function main_illustrate_mesh

clear
clc


coeffs=get_Hermitian_pol_coeffs_3d();


%------3D-------------

nx=101;
ny=101;
nz=101;

x=linspace(0,1,nx);
y=linspace(0,1,ny);
z=linspace(0,1,nz);

v=zeros(nx,ny,nz);

for ix=1:1:nx
    for iy=1:1:ny
        for iz=1:1:nz
            v(iy,ix,iz)=compute_weight_specific(x(ix),y(iy),z(iz),1,1,1,coeffs);
        end
    end
end


p=patch(isosurface(x,y,z,v,0));








end