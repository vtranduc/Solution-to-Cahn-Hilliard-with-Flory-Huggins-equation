function [time_domain,mag_sq_log]=test105(...
    c,nx,ny,nz,x_coord,y_coord,z_coord,xlen,ylen,zlen,ci_ave,ci_fluc,...
    alpha,beta,gamma,fc,nworkers,time_domain_,mag_sq_log_,time)


result=extract_abs_result_serendipity(c,nx,ny,nz);

% result=simplify_2_regions(result,ci_ave);

first=subplot(2,2,1);

s1=slice(x_coord,y_coord,z_coord,result,[0 xlen],[0 ylen],[0 zlen]);
set(s1,'FaceColor','interp')
colormap(first,'gray')
caxis([ci_ave-ci_fluc ci_ave+ci_fluc])
xlabel('x');ylabel('y');zlabel('z')

second=subplot(2,2,2);

s2=slice(x_coord,y_coord,z_coord,result,xlen/2,ylen/2,zlen/2);
set(s2,'FaceColor','interp')
colormap(second,'gray')
caxis([ci_ave-ci_fluc ci_ave+ci_fluc])
xlabel('x');ylabel('y');zlabel('z')

subplot(2,2,3)

% disp('Computing')
[frequency_domain,magnitude]=fourier_analysis_3d(result,ci_ave,alpha,beta,gamma,fc,nworkers);
% disp('good')

plot(frequency_domain,magnitude)
xlabel('Frequency domain, k')
ylabel('Structure factor, S')
grid on
grid minor

subplot(2,2,4)

mag_sq_log=[mag_sq_log_ log(max(magnitude))];
time_domain=[time_domain_ time];

plot(time_domain,mag_sq_log,'*')
xlabel('Time, t')
ylabel('Log of intensity')
grid on
grid minor

drawnow

return

%============================== z-
surface=result(:,:,1);
z_=ones(ny,nx)*min(z_coord);
redChannel = surface > ci_ave;
greenChannel = surface < ci_ave;
blueChannel = greenChannel;
surface = double(cat(3, redChannel, greenChannel, blueChannel));
surf(x_coord,y_coord,z_,surface);
hold on
%============================== z+
z_=ones(ny,nx)*max(z_coord);
surface=result(:,:,nz);
redChannel = surface > ci_ave;
greenChannel = surface < ci_ave;
blueChannel = greenChannel;
surface = double(cat(3, redChannel, greenChannel, blueChannel));
surf(x_coord,y_coord,z_,surface);
shading interp
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
redChannel = surface > ci_ave;
greenChannel = surface < ci_ave;
blueChannel = greenChannel;
surface = double(cat(3, redChannel, greenChannel, blueChannel));
surf(x_,y_,z_,surface)
%============================== y+
y_=ones(nz,nx)*max(y_coord);
surface=zeros(nz,nx);
for i=1:1:nx
    for j=1:1:nz
        surface(j,i)=result(ny,i,j);
    end
end
redChannel = surface > ci_ave;
greenChannel = surface < ci_ave;
blueChannel = greenChannel;
surface = double(cat(3, redChannel, greenChannel, blueChannel));
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
redChannel = surface > ci_ave;
greenChannel = surface < ci_ave;
blueChannel = greenChannel;
surface = double(cat(3, redChannel, greenChannel, blueChannel));
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
redChannel = surface > ci_ave;
greenChannel = surface < ci_ave;
blueChannel = greenChannel;
surface = double(cat(3, redChannel, greenChannel, blueChannel));
surf(x_,y_,z_,surface)
%==============================

hold off
drawnow
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
% colormap default
% colorbar
% title('Cubic illustration')
% xlabel('x');ylabel('y');zlabel('z');
% % colorbar
% frame_cap=getframe(gcf);

%=====================================

% face=result(:,:,2);
% 
% redChannel = face > ci_ave;
% greenChannel = face < ci_ave;
% blueChannel = greenChannel;
% colors = double(cat(3, redChannel, greenChannel, blueChannel));
% surf(x,y,zlen*ones(ny,nx),colors)
% axis([0 xlen 0 ylen 0 zlen])
% hold on
% face=zeros(nz,nx);
% for i=1:1:nx
%     for j=1:1:nz
%         face(j,i)=result(1,i,j);
%     end
% end
% redChannel = face > ci_ave;
% greenChannel = face < ci_ave;
% blueChannel = greenChannel;
% colors = double(cat(3, redChannel, greenChannel, blueChannel));
% 
% pos=zeros(nz,nx);
% for i=2:1:nz
%     pos(i,:)=z(i)*ones(1,nx);
% end
% 
% pos
% 
% surf(x,z,pos,colors)
% 
% 
% hold off
end