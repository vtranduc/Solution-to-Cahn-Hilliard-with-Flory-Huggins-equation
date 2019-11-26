function [DX,DY]=extract_gradient_2D(c,nx,ny)
n=nx*ny;
DX=zeros(1,n);
k=0;
for i=2:4:4*n-2
    k=k+1;
    DX(k)=c(i);
end
DX=reshape(DX,[ny,nx]);
DY=zeros(1,n);
k=0;
for i=3:4:4*n-1
    k=k+1;
    DY(k)=c(i);
end
DY=reshape(DY,[ny,nx]);
end