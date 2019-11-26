function sol=rotate_sph(n_nSurf,n,X,Y,Z)
%Cartesian to Spherical
sol=zeros(n_nSurf(6),3,3);
parfor i=1:1:n_nSurf(6)
    spherical_coord=cartesian_to_spherical(get_xyz(i+n,X,Y,Z));
    sol(i,:,:)=rotation_matrix(spherical_coord);
end
end

function sol=rotation_matrix(spherical_coord)
sol=zeros(3,3);
sol(1,1)=cos(spherical_coord(2))*sin(spherical_coord(3));
sol(1,2)=sin(spherical_coord(2))*sin(spherical_coord(3));
sol(1,3)=cos(spherical_coord(3));
sol(2,1)=-spherical_coord(1)*sin(spherical_coord(2))*sin(spherical_coord(3));
sol(2,2)=spherical_coord(1)*cos(spherical_coord(2))*sin(spherical_coord(3));
sol(3,1)=spherical_coord(1)*cos(spherical_coord(2))*cos(spherical_coord(3));
sol(3,2)=spherical_coord(1)*sin(spherical_coord(2))*cos(spherical_coord(3));
sol(3,3)=-spherical_coord(1)*sin(spherical_coord(3));
end

% function sol=rotation_matrix(spherical_coord)
% sol=zeros(3,3);
% sol(1,1)=cos(spherical_coord(2))*sin(spherical_coord(3));
% sol(1,2)=sin(spherical_coord(2))*sin(spherical_coord(3));
% sol(1,3)=cos(spherical_coord(3));
% sol(2,1)=-sin(spherical_coord(2))/(spherical_coord(1)*sin(spherical_coord(3)));
% sol(2,2)=cos(spherical_coord(2))/(spherical_coord(1)*sin(spherical_coord(3)));
% sol(3,1)=cos(spherical_coord(2))*cos(spherical_coord(3))/spherical_coord(1);
% sol(3,2)=sin(spherical_coord(2))*cos(spherical_coord(3))/spherical_coord(1);
% sol(3,3)=-sin(spherical_coord(3))/spherical_coord(1);
% end