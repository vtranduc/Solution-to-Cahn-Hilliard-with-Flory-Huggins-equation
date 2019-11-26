function isosurf(x,y,z,V,view_angle,transparency,isovals,face_colors,edge_color)

%You can only specify one edge color
%face_colors, if multiple isovals is specified, must be in {} backet,
%in which each element is character vector in ''
%For example, {'red','blue'}

for i=1:1:length(isovals)
    p=patch(isosurface(x,y,z,V,isovals(i)));
    p.FaceColor=char(face_colors(i));
    p.EdgeColor=edge_color;
end

axis([min(x) max(x) min(y) max(y) min(z) max(z)])
view(view_angle)
alpha(transparency)

grid on

end