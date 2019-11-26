function frame_cap=visual_surf_3D(c,nx,ny,nz,x_coord,y_coord,z_coord)

result=zeros(ny,nx,nz);

l=-7;
for k=1:1:nz
    for n=1:1:nx
        for m=1:1:ny
            l=l+8;
            result(m,n,k)=c(l);
        end
    end
end

subplot(3,2,1)
surface=result(:,:,1);
surface=surface';
contourf(y_coord,x_coord,surface);
xlabel('y');ylabel('x');
title('xy facing -z')
colorbar
grid on

subplot(3,2,2)
surface=result(:,:,nz);
contourf(x_coord,y_coord,surface);
xlabel('x');ylabel('y');
title('xy facing +z')
colorbar

subplot(3,2,3)
surface=zeros(nz,nx);
for i=1:1:nx
    for j=1:1:nz
        surface(j,i)=result(1,i,j);
    end
end
contourf(x_coord,z_coord,surface);
xlabel('x');ylabel('z');
title('xz facing -y')
colorbar

subplot(3,2,4)
surface=zeros(nx,nz);
for i=1:1:nz
    for j=1:1:nx
        surface(j,i)=result(ny,j,i);
    end
end
contourf(z_coord,x_coord,surface);
xlabel('z');ylabel('x');
title('xz facing +y')
colorbar

subplot(3,2,5)
surface=zeros(ny,nz);
for i=1:1:nz
    for j=1:1:ny
        surface(j,i)=result(j,1,i);
    end
end
contourf(z_coord,y_coord,surface);
xlabel('z');ylabel('y');
title('yz facing -x')
colorbar

subplot(3,2,6)
surface=zeros(nz,ny);
for i=1:1:ny
    for j=1:1:nz
        surface(j,i)=result(i,nx,j);
    end
end
contourf(y_coord,z_coord,surface);
xlabel('y');ylabel('z');
title('yz facing +x')
colorbar

frame_cap=getframe(gcf);

end