function test54

clear
clc

nex=5;
ney=4;
nez=3;

nx=nex+1;
ny=ney+1;
nz=nez+1;
n=nx*ny*nz;
ne=nex*ney*nez;

frame_pos=[0 0];
frame_size=[720 720];
ifig_simulation=2;

face_colors={'yellow','blue','red'};
% The order is {Mean Lower Higher}
view_angle=3;
transparency=0.9;
edge_color='none';
ifig_conclusion=3;

domain_EdgeColor='green';
domain_FaceColor='white';
domain_face_transparency=0.1;
domain_edge_transparency=0.8;



[X,Y,Z]=generate_sphere_mesh_3d(nx,ny,nz,n);

subplot(1,2,1)
plot3(X,Y,Z,'.')
grid on
xlabel('x');ylabel('y');zlabel('z');

[x_minus_x,x_minus_y,x_minus_z,x_plus_x,x_plus_y,x_plus_z,...
    y_minus_x,y_minus_y,y_minus_z,y_plus_x,y_plus_y,y_plus_z,...
    z_minus_x,z_minus_y,z_minus_z,z_plus_x,z_plus_y,z_plus_z,...
    x_plot,y_plot,z_plot]=...
    setUp_illustration_curved_structure_3d(nx,ny,nz,X,Y,Z);


subplot(1,2,2)
surf(z_minus_x,z_minus_y,z_minus_z,'EdgeColor',domain_EdgeColor,'FaceColor',domain_FaceColor,'FaceAlpha',domain_face_transparency,'EdgeAlpha',domain_edge_transparency)
hold on
surf(z_plus_x,z_plus_y,z_plus_z,'EdgeColor',domain_EdgeColor,'FaceColor',domain_FaceColor,'FaceAlpha',domain_face_transparency,'EdgeAlpha',domain_edge_transparency)
surf(y_minus_x,y_minus_y,y_minus_z,'EdgeColor',domain_EdgeColor,'FaceColor',domain_FaceColor,'FaceAlpha',domain_face_transparency,'EdgeAlpha',domain_edge_transparency)
surf(y_plus_x,y_plus_y,y_plus_z,'EdgeColor',domain_EdgeColor,'FaceColor',domain_FaceColor,'FaceAlpha',domain_face_transparency,'EdgeAlpha',domain_edge_transparency)
surf(x_minus_x,x_minus_y,x_minus_z,'EdgeColor',domain_EdgeColor,'FaceColor',domain_FaceColor,'FaceAlpha',domain_face_transparency,'EdgeAlpha',domain_edge_transparency)
surf(x_plus_x,x_plus_y,x_plus_z,'EdgeColor',domain_EdgeColor,'FaceColor',domain_FaceColor,'FaceAlpha',domain_face_transparency,'EdgeAlpha',domain_edge_transparency)
hold off


e=1;

% [eX,eY,eZ]=get_eXYZ_3d(e,nex,nx,ny,nex*ney,X,Y,Z);
% hold on
% plot3(eX,eY,eZ,'ro')
% hold off

n_outer_pts=get_n_outer_pts_sphere_3d(nx,ny,nz)


rs=0.2;


index=0;
for node=1:1:n
    [xth,yth,zth]=get_xyzth_3d(node,nx,ny);
    if isEdge_3d(xth,yth,zth,nx,ny,nz)==1
        index=index+1;
        x=X(node);
        y=Y(node);
        z=Z(node);
        
        [X_extra(index),Y_extra(index),Z_extra(index)]=radial_stretch_3d(x,y,z,rs);
        
        
        
    end

end

hold on
plot3(X_extra,Y_extra,Z_extra,'r.')
hold off


end

function [x_,y_,z_]=radial_stretch_3d(x,y,z,rs)
len=sqrt(x^2+y^2+z^2);
x_=x*rs/len+x;
y_=y*rs/len+y;
z_=z*rs/len+z;
end

function status=isEdge_3d(xth,yth,zth,nx,ny,nz)
if xth==1 || xth==nx || yth==1 || yth==ny || zth==1 || zth==nz
    status=1;
else
    status=0;
end
end

function sol=get_n_outer_pts_sphere_3d(nx,ny,nz)
nx_=nx-2;ny_=ny-2;nz_=nz-2;
sol=2*(nx_*ny_+nx_*nz_+ny_*nz_)+4*(nx_+ny_+nz_)+8;
end