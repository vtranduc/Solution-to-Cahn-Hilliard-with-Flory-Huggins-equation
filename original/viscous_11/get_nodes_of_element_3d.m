function nodes=get_nodes_of_element_3d(e,nex,nexney,nxny,nx)
nodes=zeros(1,8);
[n1xth,n1yth,n1zth]=get_n1xyzth_3d(e,nex,nexney);
nodes(1)=(n1zth-1)*nxny+nx*(n1yth-1)+n1xth;
nodes(2)=nodes(1)+1;
nodes(3)=nodes(1)+nx;
nodes(4)=nodes(3)+1;
nodes(5)=nodes(1)+nxny;
nodes(6)=nodes(2)+nxny;
nodes(7)=nodes(3)+nxny;
nodes(8)=nodes(4)+nxny;
end