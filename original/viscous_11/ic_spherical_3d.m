function co=ic_spherical_3d(co,iOnes,surface_derivatives,gbfs1)
index=0;
for i=iOnes
    index=index+1;
    co(gbfs1(i))=surface_derivatives(index);
end
end