function [sf,sj]=bc_serendipity_3d(sf_,sj_,iZeros,iOnes,gbfs1)
sf=sf_;
sj=sj_;
for i=1:1:length(iZeros)
    sj(iZeros(i))=0;
end
for i=1:1:length(iOnes)
    sj(iOnes(i))=1;
    sf(gbfs1(iOnes(i)))=0;
end

end