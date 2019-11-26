function [xth,yth,zth]=get_root_of_extrusion_sph(node_,n_nSurf,n,nx,ny,nz,nx_)
node=node_-n;
if node<=n_nSurf(1)
    [xth,yth,zth]=get_xyzth_s1(node,ny);
elseif node<=n_nSurf(2)
    [xth,yth,zth]=get_xyzth_s2(node,nx,ny,n_nSurf);
elseif node<=n_nSurf(3)
    [xth,yth,zth]=get_xyzth_s3(node,nx_,n_nSurf);
elseif node<=n_nSurf(4)
    [xth,yth,zth]=get_xyzth_s4(node,ny,nx_,n_nSurf);
elseif node<=n_nSurf(5)
    [xth,yth,zth]=get_xyzth_s5(node,nx_,n_nSurf);
elseif node<=n_nSurf(6)
    [xth,yth,zth]=get_xyzth_s6(node,nz,nx_,n_nSurf);
end
end

function [xth,yth,zth]=get_xyzth_s1(node,ny)
xth=1;
yth=mod_sp(node,ny);
zth=(node-yth)/ny+1;
end

function [xth,yth,zth]=get_xyzth_s2(node_,nx,ny,n_nSurf)
xth=nx;
node=node_-n_nSurf(1);
yth=mod_sp(node,ny);
zth=(node-yth)/ny+1;
end

function [xth,yth,zth]=get_xyzth_s3(node_,nx_,n_nSurf)
yth=1;
node=node_-n_nSurf(2);
xth=mod_sp(node,nx_);
zth=(node-xth)/nx_+1;
xth=xth+1;
end

function [xth,yth,zth]=get_xyzth_s4(node_,ny,nx_,n_nSurf)
yth=ny;
node=node_-n_nSurf(3);
xth=mod_sp(node,nx_);
zth=(node-xth)/nx_+1;
xth=xth+1;
end

function [xth,yth,zth]=get_xyzth_s5(node_,nx_,n_nSurf)
zth=1;
node=node_-n_nSurf(4);
xth=mod_sp(node,nx_);
yth=(node-xth)/nx_+2;
xth=xth+1;
end

function [xth,yth,zth]=get_xyzth_s6(node_,nz,nx_,n_nSurf)
zth=nz;
node=node_-n_nSurf(5);
xth=mod_sp(node,nx_);
yth=(node-xth)/nx_+2;
xth=xth+1;
end