function laplacian=approximate_laplacian_2d(DX,DY,nx,ny,nex,ney)
DXX=zeros(ny,nx);
for inx=2:1:nx-1
    for iny=1:1:ny
        DXX(iny,inx)=DX(iny,inx+1)-DX(iny,inx-1);
    end
end
nx1=nx-1;nx2=nx-2;
for iny=1:1:ny
    DXX(iny,1)=4*DX(iny,2)-DX(iny,3);
    DXX(iny,nx)=-4*DX(iny,nx1)+DX(iny,nx2);
end
DXX=DXX*nex;
DYY=zeros(ny,nx);
for iny=2:1:ny-1
    for inx=1:1:nx
        DYY(iny,inx)=DY(iny+1,inx)-DY(iny-1,inx); 
    end
end
ny1=ny-1;ny2=ny-2;
for inx=1:1:nx
    DYY(1,inx)=4*DY(2,inx)-DY(3,inx);
    DYY(ny,inx)=-4*DY(ny1,inx)+DY(ny2,inx);
end
DYY=DYY*ney;
laplacian=(DXX+DYY)/2;
end