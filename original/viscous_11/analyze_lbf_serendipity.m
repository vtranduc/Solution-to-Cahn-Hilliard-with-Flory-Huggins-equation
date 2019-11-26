function [orientation,type]=analyze_lbf_serendipity(lbf)
type=mod(lbf,4);
if type==0
    type=4;
end
orientation=(lbf-type)/4+1;
end