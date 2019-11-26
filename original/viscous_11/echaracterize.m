function [gbasis,dx,dy,dz]=echaracterize(e,nex,ney,~,nx,ny,~,x,y,z)

nodes=zeros(1,8);
gbasis=zeros(1,64);

nxy=nx*ny;

z1=mod(e,ney);

if z1==0
    z1=ney;
end

z2=(e-z1)/ney;
z3=mod(z2,nex);
z4=(z2-z3)/nex;
            
nodes(1)=nx*ny*z4+ny*z3+z1;
nodes(2)=nodes(1)+1;
nodes(3)=nodes(1)+ny;
nodes(4)=nodes(3)+1;
nodes(5)=nodes(1)+nxy;
nodes(6)=nodes(5)+1;
nodes(7)=nodes(5)+ny;
nodes(8)=nodes(7)+1;
            
for i=1:1:8
    z1=(i-1)*8+1;
    gbasis(z1:z1+7)=nodes(i)*8-7:1:nodes(i)*8;
end
            
dx=x(nodes(3))-x(nodes(1));
dy=y(nodes(2))-y(nodes(1));
dz=z(nodes(5))-z(nodes(1));

end