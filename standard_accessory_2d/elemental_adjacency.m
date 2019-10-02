function [adjacency,inx,iny]=elemental_adjacency(node,ny,n)
iny=mod(node,ny);
if iny==0
    iny=ny;
end
inx=(node-iny)/ny+1;
adjacency=ones(1,4);
[edge_xdir,edge_ydir]=edgy(node,ny,n);
if edge_xdir==-1
    adjacency(1:1:2)=[0 0];
elseif edge_xdir==1
    adjacency(3:1:4)=[0 0];
end
if edge_ydir==-1
    adjacency(1)=0;
    adjacency(3)=0;
elseif edge_ydir==1
    adjacency(2)=0;
    adjacency(4)=0;
end
ney=ny-1;
if adjacency(1)==1
    adjacency(1)=(inx-2)*ney+iny-1;
end
if adjacency(2)==1
    adjacency(2)=(inx-2)*ney+iny;
end
if adjacency(3)==1
    adjacency(3)=(inx-1)*ney+iny-1;
end
if adjacency(4)==1
    adjacency(4)=(inx-1)*ney+iny;
end
end