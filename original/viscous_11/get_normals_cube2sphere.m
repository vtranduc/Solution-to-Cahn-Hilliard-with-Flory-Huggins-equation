function sol=get_normals_cube2sphere(nex,ney,nez,nx,ny,nz)

nexney=nex*ney;
nexnez=nex*nez;
neynez=ney*nez;

sol=zeros(2*(nexney+nexnez+neynez),3,3,3);


xth=1;
yth=1;
zth=1;

e=get_elements_3d(xth,yth,zth,nex,ney,nx,ny,nz,8)

index=0;

for zth=1:1:nez
    for yth=1:1:

end