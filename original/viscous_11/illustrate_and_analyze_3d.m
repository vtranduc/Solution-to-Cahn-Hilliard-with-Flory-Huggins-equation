function frame_cap=illustrate_and_analyze_3d(c,isovals,x_coord,y_coord,z_coord,...
    face_colors,edge_color,nx,ny,nz,view_angle,transparency,time)
subplot(1,2,1)
result=extract_abs_result_rec_3d(c,nx,ny,nz);

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

str=sprintf('Nonlinear Cahn Hilliard in 3D\n t=%d',time);
suptitle(str)



%-----------------------------

% alpha_=200;
% beta=8;
% fc=0;
% 
% if fc==0
%     fc=1/(2*min([1/(nx-1) 1/(ny-1)]));
% end
% 
% co=mean(mean(mean(result)));
% 
% c_=zeros(nx,ny);
% 
% take=5;
% 
% for i=1:1:nx
%     for j=1:1:ny
%         c_(i,j)=result(j,i,take);
%     end
% end
% 
% st=structure_factor_2D(alpha_,beta,fc,c_,co);
% 
% subplot(1,2,2)
% 
% 
% freq_domain=linspace(0,fc,alpha_);
% 
% plot(st)

fc=1/(2*min([1/(nx-1) 1/(ny-1)]));

fc=20;

st=structure_factor_3d(100,8,4,fc,result,isovals(2));
subplot(1,2,2)
frequency_domain=linspace(0,fc,100);
plot(frequency_domain,st)

%------------------------------

frame_cap=getframe(gcf);

end