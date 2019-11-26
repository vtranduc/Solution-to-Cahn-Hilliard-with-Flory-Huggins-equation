function co=create_sine_wave_2d(nx,ny,nfour,A,w,x_coord,y_coord,ci)
co=zeros(1,nfour);
index=-3;
term=2*pi*w;
for ix=1:1:nx
    for iy=1:1:ny
        index=index+4;
        co(index)=A*sin(term*x_coord(ix))*sin(term*y_coord(iy))+ci;
    end
end
end