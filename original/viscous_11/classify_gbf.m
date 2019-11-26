function [nodeID,class]=classify_gbf(gbf)
class=mod(gbf,4);
if class==0
    class=4;
end
nodeID=(gbf-class)/4+1;
end