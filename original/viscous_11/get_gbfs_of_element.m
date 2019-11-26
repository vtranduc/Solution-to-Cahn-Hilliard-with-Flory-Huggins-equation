function solution=get_gbfs_of_element(element,nexney,nex,nx,ny)
%nexney=nex*ney
solution=zeros(1,64);
xythFirst=mod(element,nexney);
if xythFirst==0
    xythFirst=nexney;
end
zthFirst=(element-xythFirst)/nexney+1;
xthFirst=mod(element,nex);
if xthFirst==0
    xthFirst=nex;
end
ythFirst=(xythFirst-xthFirst)/nex+1;

solution(1:8)=get_gbs_3d(get_node_3d(...
    xthFirst,ythFirst,zthFirst,nx,ny));

% solution(9:16)=get_gbs_3d(get_node_3d(...
%     xthFirst+1,ythFirst,zthFirst,nx,ny));

solution(17:24)=get_gbs_3d(get_node_3d(...
    xthFirst,ythFirst+1,zthFirst,nx,ny));

% solution(25:32)=get_gbs_3d(get_node_3d(...
%     xthFirst+1,ythFirst+1,zthFirst,nx,ny));

index1=8;
index2=24;

for i=1:1:8
    solution(index1+i)=solution(i)+8;
    solution(index2+i)=solution(16+i)+8;
end

solution(33:64)=solution(1:32)+ones(1,32)*nx*ny*8;
end