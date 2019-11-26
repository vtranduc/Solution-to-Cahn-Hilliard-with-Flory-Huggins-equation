function test46

% [x,y,z] = sphere
% 
% x'
% y'
clear
clc

[X,Y,Z] = meshgrid(-10:10,-10:10,-10:10);
V = sqrt(X.^2+Y.^2+Z.^2);


X(1,:,:);

% plot3(X,Y,Z,'go')

size(X)

% fv = isosurface(X,Y,Z,V);
% p = patch(fv);
% set(p,'EdgeColor','g')



end