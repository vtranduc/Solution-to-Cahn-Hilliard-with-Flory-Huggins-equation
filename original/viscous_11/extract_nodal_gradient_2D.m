function [d_dx,d_dy]=extract_nodal_gradient_2D(c,nx,ny)

n=nx*ny;

d_dx=zeros(1,n);
k=0;
for i=2:4:4*n
    k=k+1;
    d_dx(k)=c(i);
end
d_dx=reshape(d_dx,[ny,nx]);

d_dy=zeros(1,n);
k=0;
for i=3:4:4*n
    k=k+1;
    d_dy(k)=c(i);
end
d_dy=reshape(d_dy,[ny,nx]);

end