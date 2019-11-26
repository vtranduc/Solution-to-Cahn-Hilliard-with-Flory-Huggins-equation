function abs_w=extract_nodal_weights_2D(c,nx,ny)

n=nx*ny;
abs_w=zeros(1,n);
k=0;
for i=1:4:4*n-3
    k=k+1;
    abs_w(k)=c(i);
end

abs_w=reshape(abs_w,[ny,nx]);

end