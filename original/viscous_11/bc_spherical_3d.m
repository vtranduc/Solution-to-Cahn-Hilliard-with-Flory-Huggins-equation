function [sf,sj]=bc_spherical_3d(sf,sj,iOnes,iZeros,surface_derivatives,gbfs1)
for i=iZeros
    sj(i)=0;
end
index=0;
for i=iOnes
    index=index+1;
    sj(i)=1;
    sf(gbfs1(i))=surface_derivatives(index);
end
end