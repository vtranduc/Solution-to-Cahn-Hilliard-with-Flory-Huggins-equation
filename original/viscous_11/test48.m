function test48

clear
clc

% [x y z v] = flow;
% q = z./x.*y.^3;
% p = patch(isosurface(x, y, z, q, -.08, v));
% isonormals(x,y,z,q, p)
% p.FaceColor = 'interp';
% p.EdgeColor = 'none';
% daspect([1 1 1]); axis tight; 
% colormap(prism(28))
% camup([1 0 0 ]); campos([25 -55 5]) 
% camlight; lighting phong


[x y z v] = flow;
p = patch(isosurface(x, y, z, v, -3));
isonormals(x,y,z,v, p)
p.FaceColor = 'red';
p.EdgeColor = 'green';
daspect([1 1 1])
view(3)
camlight; lighting phong



size(x)
size(y)
size(z)
size(v)

end