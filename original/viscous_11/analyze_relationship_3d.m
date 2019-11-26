function analyze_relationship_3d(xth,yth,zth,nodeth,relationth,...
    nex,ney,nx,ny,nz)

%THIS FUNCTION IS INCOMPLETE AND IS NO LONGER NEEDED

%YOU CAN DELETE THIS FILE IF YOU WANT

adjacency=get_elements_3d(xth,yth,zth,nex,ney,nx,ny,nz);

node_existence=zeros(1,27);

if adjacency(1)~=0
    node_existence(1)=1;
end
if adjacency(1)~=0 || adjacency(2)~=0
    node_existence(2)=1;
end
%====
if adjacency(3)~=0
    node_existence(3)=1;
end
if adjacency(1)~=0 || adjacency(3)~=0
    node_existence(4)=1;
end
if zth~=1
    node_existence(5)=1;
end
if adjacency(2)~=0 || adjacency(4)~=0
    node_existence(6)=1;
end
if adjacency(3)~=0
    node_existence(7)=1;
end
if adjacency(3)~=0 || adjacency(4)~=0
    node_existence(8)=1;
end
if adjacency(4)~=0
    node_existence(9)=1;
end
%===
if adjacency(1)~=0 || adjacency(5)~=0
    node_existence(10)=1;
end
if yth~=1
    node_existence(11)=1;
end
if adjacency(2)~=0 || adjacency(6)~=0
    node_existence(12)=1;
end
%====
if adjacency(1)~=0 || adjacency(2)~=0
    node_existence(13)=1;
end
if adjacency(1)~=0 || adjacency(2)~=0
    node_existence(14)=1;
end
if adjacency(1)~=0 || adjacency(2)~=0
    node_existence(15)=1;
end
if adjacency(1)~=0 || adjacency(2)~=0
    node_existence(16)=1;
end
if adjacency(1)~=0 || adjacency(2)~=0
    node_existence(17)=1;
end
if adjacency(1)~=0 || adjacency(2)~=0
    node_existence(18)=1;
end
if adjacency(1)~=0 || adjacency(2)~=0
    node_existence(19)=1;
end
if adjacency(1)~=0 || adjacency(2)~=0
    node_existence(20)=1;
end
if adjacency(1)~=0 || adjacency(2)~=0
    node_existence(21)=1;
end
if adjacency(1)~=0 || adjacency(2)~=0
    node_existence(22)=1;
end
if adjacency(1)~=0 || adjacency(2)~=0
    node_existence(23)=1;
end
if adjacency(1)~=0 || adjacency(2)~=0
    node_existence(24)=1;
end
if adjacency(1)~=0 || adjacency(2)~=0
    node_existence(25)=1;
end
if adjacency(1)~=0 || adjacency(2)~=0
    node_existence(26)=1;
end
if adjacency(1)~=0 || adjacency(2)~=0
    node_existence(27)=1;
end

end