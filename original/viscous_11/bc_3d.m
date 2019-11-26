function [sf,sj]=bc_3d(sf,sj,iZeros,iOnes,gbfs1)
for i=iZeros
    sj(i)=0;
end
for i=iOnes
    sj(i)=1;
    sf(gbfs1(i))=0;
end
end