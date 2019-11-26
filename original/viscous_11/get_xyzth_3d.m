function [xth,yth,zth]=get_xyzth_3d(node,nx,ny)
nxny=nx*ny;
xyth=mod(node,nxny);
if xyth==0
    xyth=nxny;
end
zth=(node-xyth)/nxny+1;
xth=mod(xyth,nx);
if xth==0
    xth=nx;
end
yth=(xyth-xth)/nx+1;
end