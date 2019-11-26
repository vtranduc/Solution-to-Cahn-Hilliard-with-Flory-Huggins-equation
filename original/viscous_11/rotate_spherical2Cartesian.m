function cr=rotate_spherical2Cartesian(c,neight,n_nSurf,inverse_rotation)

cr=c;
first_index=neight-7;
for i=1:1:n_nSurf(6)
    first_index=first_index+8;
    for j=1:1:3
        cr(first_index+j)=...
            inverse_rotation(i,j,1)*c(first_index+1)...
            +inverse_rotation(i,j,2)*c(first_index+2)...
            +inverse_rotation(i,j,3)*c(first_index+3);
    end
end
end