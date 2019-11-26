function frame_cap=visual_rec_3D(c,ci_ave,nx,ny,nz,x_coord,y_coord,z_coord)

result=extract_abs_result_rec_3D(c,nx,ny,nz);

% result=simplify_2_regions(result,ci_ave);

%============================== z-
surface=result(:,:,1);
z_=ones(ny,nx)*min(z_coord);
surf(x_coord,y_coord,z_,surface);
hold on
%============================== z+
z_=ones(ny,nx)*max(z_coord);
surface=result(:,:,nz);
surf(x_coord,y_coord,z_,surface);
%============================== y-
x_=zeros(nz,nx);
z_=zeros(nz,nx);
for i=1:1:nz
    x_(i,:)=x_coord;
    z_(i,:)=ones(1,nx)*z_coord(i);
end
y_=ones(nz,nx)*min(y_coord);
surface=zeros(nz,nx);
for i=1:1:nx
    for j=1:1:nz
        surface(j,i)=result(1,i,j);
    end
end
surf(x_,y_,z_,surface)
%============================== y+
y_=ones(nz,nx)*max(y_coord);
surface=zeros(nz,nx);
for i=1:1:nx
    for j=1:1:nz
        surface(j,i)=result(ny,i,j);
    end
end
surf(x_,y_,z_,surface)
%============================== x-
x_=ones(nz,ny)*min(x_coord);
y_=zeros(nz,ny);
z_=zeros(nz,ny);
for i=1:1:nz
    y_(i,:)=y_coord;
    z_(i,:)=ones(1,ny)*z_coord(i);
end
surface=zeros(nz,ny);
for i=1:1:ny
    for j=1:1:nz
        surface(j,i)=result(i,1,j);
    end
end
surf(x_,y_,z_,surface)
%============================== x+
x_=ones(nz,ny)*max(x_coord);
y_=zeros(nz,ny);
z_=zeros(nz,ny);
for i=1:1:nz
    y_(i,:)=y_coord;
    z_(i,:)=ones(1,ny)*z_coord(i);
end
surface=zeros(nz,ny);
for i=1:1:ny
    for j=1:1:nz
        surface(j,i)=result(i,nx,j);
    end
end
surf(x_,y_,z_,surface)
%==============================

hold off
% map = [0, 0, 0.3
%     0, 0, 0.4
%     0, 0, 0.5
%     0, 0, 0.6
%     0, 0, 0.8
%     0, 0, 1.0];
% n=10000;
% map=zeros(n,3);
% mapper=linspace(0,1,n);
% for i=1:1:n
%     map(i,3)=mapper(i);
% end
% % colormap(parula(5))
% colormap(map)
%===================================
% caxis([0 1])
%===================================
% caxis([ci_ave-0.01 ci_ave+0.01])
colormap default
colorbar
title('Cubic illustration')
xlabel('x');ylabel('y');zlabel('z');
% colorbar
frame_cap=getframe(gcf);

end