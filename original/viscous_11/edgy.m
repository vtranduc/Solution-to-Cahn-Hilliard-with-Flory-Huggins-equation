function [edge_xdir,edge_ydir]=edgy(node,ny,n)
if mod(node-1,ny)==0
    edge_ydir=-1;
elseif mod(node,ny)==0
    edge_ydir=1;
else
    edge_ydir=0;
end
if node<=ny
    edge_xdir=-1;
elseif node>n-ny;
    edge_xdir=1;
else
    edge_xdir=0;
end
end