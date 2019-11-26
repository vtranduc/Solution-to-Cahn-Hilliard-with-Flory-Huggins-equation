function [n1xth,n1yth,n1zth]=get_n1xyzth_3d(e,nex,nexney)
n1xth=mod(e,nex);
if n1xth==0
    n1xth=nex;
end
n1xyth=mod(e,nexney);
if n1xyth==0
    n1xyth=nexney;
end
n1zth=(e-n1xyth)/nexney+1;
n1yth=(n1xyth-n1xth)/nex+1;
end