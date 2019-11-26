function [eX,eY,eZ]=get_eXYZ_3d(e,nex,nx,ny,nexney,X,Y,Z)
eX=zeros(1,8);
eY=zeros(1,8);
eZ=zeros(1,8);
[n1xth,n1yth,n1zth]=get_n1xyzth_3d(e,nex,nexney);
nodes(1)=get_node_3d(n1xth,n1yth,n1zth,nx,ny);
nodes(2)=nodes(1)+1;
nodes(3)=get_node_3d(n1xth,n1yth+1,n1zth,nx,ny);
nodes(4)=nodes(3)+1;
nodes(5)=get_node_3d(n1xth,n1yth,n1zth+1,nx,ny);
nodes(6)=nodes(5)+1;
nodes(7)=get_node_3d(n1xth,n1yth+1,n1zth+1,nx,ny);
nodes(8)=nodes(7)+1;
for i=1:1:8
    eX(i)=X(nodes(i));
    eY(i)=Y(nodes(i));
    eZ(i)=Z(nodes(i));
end
end