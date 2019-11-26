function c_adjusted=bc_for_c_serendipity_wetting_3d(c,bc_wetting)
c_adjusted=c;
for i=1:1:length(bc_wetting)
    c_adjusted(bc_wetting(1,i))=bc_wetting(3,i)+c(bc_wetting(2,i))*bc_wetting(4,i);
end
end