function test

clear
clc

% x=0:1:10;
% y=x;
% z=x;
% 
% % generate some data
% x = rand(1,100);
% y = rand(1,100);
% z = rand(1,100);
% i = rand(1,100)*200;
% 
% % specify the indexed color for each point
% icolor = ceil((i/max(i))*256);
% 
% figure;
% scatter3(x,y,z,i,icolor);
% figure;
% scatter3(x,y,z,i,icolor,'filled');

[x y z v] = flow;
       p = patch(isosurface(x, y, z, v, -3));
       %isonormals(x,y,z,v, p)
       p.FaceColor = 'red';
       %p.FaceColor=[0,0.5,0]
       %p.EdgeColor = 'none';
       %daspect([1 1 1])
       alpha(0.3)
       view(3)
       %camlight; lighting phong
       grid on
       
       p = patch(isosurface(x, y, z, v, -1));
       %isonormals(x,y,z,v, p)
       p.FaceColor = 'blue';
       alpha(0.3)
       
       size(x)
       size(y)
       size(z)
       size(v)
       


% x=linspace(0,1,100);
% 
% plot(x,basis(x,0,0),x,basis(x,1,0),x,basis(x,0,1),x,basis(x,1,1))
% 
% grid on
% 
% figure(2)
% 
% basis(x,0,0)+basis(x,1,0)
% 
% plot(x,basis(x,0,1)+basis(x,1,1))
% 
% grid on

% x=linspace(0,1,10);
% y=linspace(0,1,10);
% z=linspace(0,1,10);
% 
% Z=zeros(10,10);
% 
% for i=1:10
%     for j=1:10
%         Z(j,i)=basis(x(i),0,0)*basis(y(j),0,0);
%     end
% end
% 
% surf(x,y,Z)
% xlabel('x'),ylabel('y')


end