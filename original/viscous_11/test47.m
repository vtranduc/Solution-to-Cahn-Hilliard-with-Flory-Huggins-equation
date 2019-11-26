function frame_cap=test47(c,isovals,fig_domain,...
    face_colors,edge_color,nx,ny,nz,view_angle,transparency,time,...
    x_minus_x,x_minus_y,x_minus_z,x_plus_x,x_plus_y,x_plus_z,...
    y_minus_x,y_minus_y,y_minus_z,y_plus_x,y_plus_y,y_plus_z,...
    z_minus_x,z_minus_y,z_minus_z,z_plus_x,z_plus_y,z_plus_z,...
    x_plot,y_plot,z_plot,domain_EdgeColor,domain_FaceColor,...
    domain_face_transparency,domain_edge_transparency)

result=extract_abs_result_rec_3d(c,nx,ny,nz);

for i=1:1:length(isovals)
    p=patch(isosurface(x_plot,y_plot,z_plot,result,isovals(i)));
    p.FaceColor=char(face_colors(i));
    p.EdgeColor=edge_color;
end
axis(fig_domain)
view(view_angle)
alpha(transparency)
camlight 
lighting gouraud

xlabel('x')
ylabel('y')
zlabel('z')
grid on

hold on
surf(z_minus_x,z_minus_y,z_minus_z,'EdgeColor',domain_EdgeColor,'FaceColor',domain_FaceColor,'FaceAlpha',domain_face_transparency,'EdgeAlpha',domain_edge_transparency)
surf(z_plus_x,z_plus_y,z_plus_z,'EdgeColor',domain_EdgeColor,'FaceColor',domain_FaceColor,'FaceAlpha',domain_face_transparency,'EdgeAlpha',domain_edge_transparency)
surf(y_minus_x,y_minus_y,y_minus_z,'EdgeColor',domain_EdgeColor,'FaceColor',domain_FaceColor,'FaceAlpha',domain_face_transparency,'EdgeAlpha',domain_edge_transparency)
surf(y_plus_x,y_plus_y,y_plus_z,'EdgeColor',domain_EdgeColor,'FaceColor',domain_FaceColor,'FaceAlpha',domain_face_transparency,'EdgeAlpha',domain_edge_transparency)
surf(x_minus_x,x_minus_y,x_minus_z,'EdgeColor',domain_EdgeColor,'FaceColor',domain_FaceColor,'FaceAlpha',domain_face_transparency,'EdgeAlpha',domain_edge_transparency)
surf(x_plus_x,x_plus_y,x_plus_z,'EdgeColor',domain_EdgeColor,'FaceColor',domain_FaceColor,'FaceAlpha',domain_face_transparency,'EdgeAlpha',domain_edge_transparency)
hold off

str=sprintf('Nonlinear Cahn Hilliard in 3D\n t=%d',time);
suptitle(str)

frame_cap=getframe(gcf);

end