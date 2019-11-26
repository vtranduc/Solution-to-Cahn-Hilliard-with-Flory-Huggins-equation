function [chi,diffT]=chi_diffT_extractor(adjusted_chi,adjusted_diffT,inx)
chi=adjusted_chi(inx,:);
diffT=adjusted_diffT(inx,:);
end