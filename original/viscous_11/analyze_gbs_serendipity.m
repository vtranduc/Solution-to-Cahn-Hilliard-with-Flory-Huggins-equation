function [node,type]=analyze_gbs_serendipity(gbf)
type=mod(gbf,4);
if type==0
    type=4;
end
node=(gbf-type)/4+1;
end