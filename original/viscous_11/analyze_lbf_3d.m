function [orientation,type]=analyze_lbf_3d(lbf)
type=mod(lbf,8);
if type==0
    type=8;
end
orientation=(lbf-type)/8+1;
end