function weights=generateWeights_serendipity_3d(xlen,ylen,zlen,nex,ney,nez)

gp=[(-sqrt(3/5)+1)/2 0.5 (sqrt(3/5)+1)/2];
dx=xlen/nex;
dy=ylen/ney;
dz=zlen/nez;

eX=[0 dx 0 dx 0 dx 0 dx];
eY=[0 0 dy dy 0 0 dy dy];
eZ=[0 0 0 0 dz dz dz dz];

weights=compute_elemental_weights_serendipity(gp,get_phis(),eX,eY,eZ);



end