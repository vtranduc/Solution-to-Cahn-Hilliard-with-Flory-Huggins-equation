function xyz=spherical_to_cartesian(spherical_coord)
xyz=zeros(1,3);
[x,y,z] = sph2cart(spherical_coord(2),pi/2-spherical_coord(3),spherical_coord(1));
xyz(1)=x;
xyz(2)=y;
xyz(3)=z;
end