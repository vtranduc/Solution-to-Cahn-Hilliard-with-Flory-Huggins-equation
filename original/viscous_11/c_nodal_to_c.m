function sol=c_nodal_to_c(c_nodal,nx,ny,nex,ney,dx,dy)

sol=zeros(4*nx*ny,1);
    
[fx,fy]=gradient(c_nodal,dx,dy);
[fxy1,~]=gradient(fy,dx,dy);
[~,fxy2]=gradient(fx,dx,dy);
fxy=(fxy1+fxy2)/2;
index=0;
for xth=1:1:nx
    for yth=1:1:ny
        index=index+1;
        sol(index)=c_nodal(yth,xth);
        index=index+1;
        sol(index)=fx(yth,xth);
        index=index+1;
        sol(index)=fy(yth,xth);
        index=index+1;
        sol(index)=fxy(yth,xth);
    end
end
ll=2;li=4;lu=2+4*ney;
bl=3;bi=4*ny;bu=3+4*nex*ny;
tl=3+4*ney;ti=4*ny;tu=3+4*ney+4*nex*ny;
rl=2+4*nex*ny;ri=4;ru=2+4*ney+4*nex*ny;
sol=bc(sol,ll,li,lu,bl,bi,bu,rl,ri,ru,tl,ti,tu);

end