function [sf,sj]=bc_traditional_3d(sf,sj,gbfs1,gbfs2,iZeros,iOnes)
for i=1:1:length(iZeros)
    sj(gbfs1(iZeros(i)),gbfs2(iZeros(i)))=0.0;
end
for i=1:1:length(iOnes)
    sj(gbfs1(iOnes(i)),gbfs2(iOnes(i)))=1.0;
    sf(gbfs1(iOnes(i)))=0.0;
end
end