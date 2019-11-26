function frame_cap=visual3D(c,isovals,x_coord,y_coord,z_coord,face_colors,edge_color,nx,ny,nz,view_angle,transparency,t)

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

for i=1:1:length(isovals)
    p=patch(isosurface(x_coord,y_coord,z_coord,result,isovals(i)));
    p.FaceColor=char(face_colors(i));
    p.EdgeColor=edge_color;
end

axis([min(x_coord) max(x_coord) min(y_coord) max(y_coord) min(z_coord) max(z_coord)])
view(view_angle)
alpha(transparency)

xlabel('x')
ylabel('y')
zlabel('z')

str=sprintf('Nonlinear Cahn Hilliard in 3D\n t=%2.8d',t);
suptitle(str)

grid on

frame_cap=getframe(gcf);

end