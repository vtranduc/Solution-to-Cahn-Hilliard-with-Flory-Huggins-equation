function spherical_coord=cartesian_to_spherical(xyz)
[azimuth,elevation,r] = cart2sph(xyz(1),xyz(2),xyz(3));
spherical_coord=[r,azimuth,pi/2-elevation];
end