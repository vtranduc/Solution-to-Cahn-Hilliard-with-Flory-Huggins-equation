function [node,type]=analyze_gbs_3d(gbs)
type=mod(gbs,8);
if type==0
    type=8;
end
node=(gbs-type)/8+1;
end