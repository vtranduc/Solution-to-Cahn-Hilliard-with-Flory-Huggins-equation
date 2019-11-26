function test106

clear
clc

[x,y,z,v] = flow;
figure;
xslice = 5;
yslice = 0;
zslice = 0;

size(x)
size(y)
size(z)
size(v)

min(min(min(x)))

% slice(x,y,z,v,xslice,yslice,zslice);
% view(3);
% axis on;
% grid on;
% light;
% lighting phong;
% camlight('left');
% shading interp;

return

vert = [0 0 0;1 0 0;1 1 0;0 1 0;0 0 1;1 0 1;1 1 1;0 1 1];
vert = 2*pi*vert;
fac = [1 2 6 5;2 3 7 6;3 4 8 7;4 1 5 8;1 2 3 4;5 6 7 8];
C = rand(8,3);
% C = hsv(8);
patch('Vertices',vert,'Faces',fac,'FaceVertexCData',C,'EdgeColor','none','FaceColor','none')
view(3)
axis vis3d
colorbar

% return


% vert = [0 0 0;1 0 0;1 1 0;0 1 0;0 0 1;1 0 1;1 1 1;0 1 1];
% vert = 2*pi*vert;
% fac = [1 2 6 5;2 3 7 6;3 4 8 7;4 1 5 8;1 2 3 4;5 6 7 8];
% C = rand(8,3);
% % C = hsv(8);
% patch('Vertices',vert,'Faces',fac,'FaceVertexCData',C,'FaceColor','interp')
% view(3)
% axis vis3d
% colorbar


clf
X = [0 3; 0 3; 5 6];
Y = [0 3; 5 6; 0 3];
C = [5 4; 2 0; 6 3];
p = patch(X,Y,C,'FaceColor','interp');
colorbar


clf
x=1:10;
f=x.^2;
g=x.^2+1;
patch([x x(end:-1:1)], [f g(end:-1:1)],'y')

clf

data=rand(5,4);
p=patch('Vertices',vert,'Faces',fac,'FaceVertexCData',contourf(data))



end