function frame_cap=visual_all_3D(...
    c,isovals,x_coord,y_coord,z_coord,face_colors,edge_color,...
    nx,ny,nz,view_angle,transparency,t,minz,maxz)

result=extract_abs_result_rec_3D(c,nx,ny,nz);

%==========================================================================

subplot(3,3,1)

l=-7;
for k=1:1:nz
    for n=1:1:nx
        for m=1:1:ny
            l=l+8;
            result(m,n,k)=c(l);
        end
    end
end

for i=1:1:length(isovals)
    p=patch(isosurface(x_coord,y_coord,z_coord,result,isovals(i)));
    p.FaceColor=char(face_colors(i));
    p.EdgeColor=edge_color;
end

axis([min(x_coord) max(x_coord) min(y_coord) max(y_coord) min(z_coord) max(z_coord)])
view(view_angle)
alpha(transparency)

camlight 
lighting gouraud

xlabel('x')
ylabel('y')
zlabel('z')

grid on

%==========================================================================

subplot(3,3,2)

%============================== z-
surface1=result(:,:,1);
z_1=ones(ny,nx)*min(z_coord);
surf(x_coord,y_coord,z_1,surface1);
hold on
%============================== z+
z_2=ones(ny,nx)*max(z_coord);
surface2=result(:,:,nz);
surf(x_coord,y_coord,z_2,surface2);
%============================== y-
x_3=zeros(nz,nx);
z_3=zeros(nz,nx);
for i=1:1:nz
    x_3(i,:)=x_coord;
    z_3(i,:)=ones(1,nx)*z_coord(i);
end
y_3=ones(nz,nx)*min(y_coord);
surface3=zeros(nz,nx);
for i=1:1:nx
    for j=1:1:nz
        surface3(j,i)=result(1,i,j);
    end
end
surf(x_3,y_3,z_3,surface3)
%============================== y+
y_4=ones(nz,nx)*max(y_coord);
surface4=zeros(nz,nx);
for i=1:1:nx
    for j=1:1:nz
        surface4(j,i)=result(ny,i,j);
    end
end
surf(x_3,y_4,z_3,surface4)
%============================== x-
x_5=ones(nz,ny)*min(x_coord);
y_5=zeros(nz,ny);
z_5=zeros(nz,ny);
for i=1:1:nz
    y_5(i,:)=y_coord;
    z_5(i,:)=ones(1,ny)*z_coord(i);
end
surface5=zeros(nz,ny);
for i=1:1:ny
    for j=1:1:nz
        surface5(j,i)=result(i,1,j);
    end
end
surf(x_5,y_5,z_5,surface5)
%============================== x+
x_6=ones(nz,ny)*max(x_coord);
y_6=zeros(nz,ny);
z_6=zeros(nz,ny);
for i=1:1:nz
    y_6(i,:)=y_coord;
    z_6(i,:)=ones(1,ny)*z_coord(i);
end
surface6=zeros(nz,ny);
for i=1:1:ny
    for j=1:1:nz
        surface6(j,i)=result(i,nx,j);
    end
end
surf(x_6,y_6,z_6,surface6)
%==============================
hold off
%===================================
% caxis([0 1])
%===================================
colormap default
colorbar
title('Cubic illustration')
xlabel('x');ylabel('y');zlabel('z');
view(view_angle)
% colorbar
% camlight 
% lighting gouraud

%==========================================================================

subplot(3,3,3)
surf(x_coord,y_coord,z_1,surface1);
hold on
surf(x_coord,y_coord,z_2,surface2);
surf(x_3,y_3,z_3,surface3);
surf(x_3,y_4,z_3,surface4)
surf(x_5,y_5,z_5,surface5)
surf(x_6,y_6,z_6,surface6)
hold off
caxis([minz maxz])
title('Cubic illustration')
xlabel('x');ylabel('y');zlabel('z');
view(view_angle)
colorbar
%==========================================================================

subplot(3,3,4)
surf(x_coord,y_coord,surface1)
xlabel('x');ylabel('y');zlabel('c');
title('Facing -z')
% colorbar('southoutside')
axis([0 1 0 1 minz maxz])

subplot(3,3,5)
surf(x_coord,y_coord,surface2)
xlabel('x');ylabel('y');zlabel('c');
title('Facing +z')
% colorbar('southoutside')
axis([0 1 0 1 minz maxz])

subplot(3,3,6)
surf(x_3,z_3,surface3)
xlabel('x');ylabel('z');zlabel('c');
title('Facing -y')
% colorbar('southoutside')
axis([0 1 0 1 minz maxz])

subplot(3,3,7)
surf(x_3,z_3,surface4)
xlabel('x');ylabel('z');zlabel('c');
title('Facing +y')
% colorbar('southoutside')
axis([0 1 0 1 minz maxz])

subplot(3,3,8)
surf(y_5,z_5,surface5)
xlabel('y');ylabel('z');zlabel('c');
title('Facing -x')
% colorbar('southoutside')
axis([0 1 0 1 minz maxz])

subplot(3,3,9)
surf(y_6,z_6,surface6)
xlabel('y');ylabel('z');zlabel('c');
title('Facing +x')
% colorbar('southoutside')
axis([0 1 0 1 minz maxz])

%==========================================================================

str=sprintf('Nonlinear Cahn Hilliard in 3D\n t=%2.8d',t);
suptitle(str)

frame_cap=getframe(gcf);

%==========================================================

end