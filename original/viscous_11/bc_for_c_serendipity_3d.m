function c_adjusted=bc_for_c_serendipity_3d(c,iOnes,gbfs1)
c_adjusted=c;
for i=1:1:length(iOnes)
    c_adjusted(gbfs1(iOnes(i)))=0;
end
end