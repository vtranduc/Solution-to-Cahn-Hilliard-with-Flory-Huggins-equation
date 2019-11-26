function nodes_adjacent=get_node_adjacent_3d(node,nx,ny,nz)
%If loop can be shortened, so that less statements will be checked.
%This can be fixed later if time permitted
nodes_adjacent=zeros(1,27);
[xth,yth,zth]=get_xyzth_3d(node,nx,ny);
if zth~=1
    node_minus_nxny=node-nx*ny;
    if xth~=1 && yth~=1
        nodes_adjacent(1)=node_minus_nxny-nx-1;
    end
    if yth~=1
        nodes_adjacent(2)=node_minus_nxny-nx;
    end
    if xth~=nx && yth~=1
        nodes_adjacent(3)=node_minus_nxny-nx+1;
    end
    if xth~=1
        nodes_adjacent(4)=node_minus_nxny-1;
    end

	nodes_adjacent(5)=node_minus_nxny;

    if xth~=nx
        nodes_adjacent(6)=node_minus_nxny+1;
    end
    if xth~=1 && yth~=ny
        nodes_adjacent(7)=node_minus_nxny+nx-1;
    end
    if yth~=ny
        nodes_adjacent(8)=node_minus_nxny+nx;
    end
    if xth~=nx && yth~=ny
        nodes_adjacent(9)=node_minus_nxny+nx+1;
    end
end
if xth~=1 && yth~=1
    nodes_adjacent(10)=node-nx-1;
end
if yth~=1
    nodes_adjacent(11)=node-nx;
end
if xth~=nx && yth~=1
    nodes_adjacent(12)=node-nx+1;
end
if xth~=1
    nodes_adjacent(13)=node-1;
end

nodes_adjacent(14)=node;

if xth~=nx
    nodes_adjacent(15)=node+1;
end
if xth~=1 && yth~=ny
    nodes_adjacent(16)=node+nx-1;
end
if yth~=ny
    nodes_adjacent(17)=node+nx;
end
if xth~=nx && yth~=ny
    nodes_adjacent(18)=node+nx+1;
end
if zth~=nz
    node_plus_nxny=node+nx*ny;
    if xth~=1 && yth~=1
        nodes_adjacent(19)=node_plus_nxny-nx-1;
    end
    if yth~=1
        nodes_adjacent(20)=node_plus_nxny-nx;
    end
    if xth~=nx && yth~=1
        nodes_adjacent(21)=node_plus_nxny-nx+1;
    end
    if xth~=1
        nodes_adjacent(22)=node_plus_nxny-1;
    end

	nodes_adjacent(23)=node_plus_nxny;

    if xth~=nx
        nodes_adjacent(24)=node_plus_nxny+1;
    end
    if xth~=1 && yth~=ny
        nodes_adjacent(25)=node_plus_nxny+nx-1;
    end
    if yth~=ny
        nodes_adjacent(26)=node_plus_nxny+nx;
    end
    if xth~=nx && yth~=ny
        nodes_adjacent(27)=node_plus_nxny+nx+1;
    end
end
end