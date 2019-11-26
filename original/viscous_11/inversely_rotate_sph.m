function sol=inversely_rotate_sph(n_nSurf,n,X,Y,Z)
%Spherical to Cartesian
sol=zeros(n_nSurf(6),3,3);
parfor i=1:1:n_nSurf(6)
    spherical_coord=cartesian_to_spherical(get_xyz(i+n,X,Y,Z));
    sol(i,:,:)=rotation_matrix(spherical_coord);
end
end

function sol=rotation_matrix(spherical_coord)
sol=zeros(3,3);
sol(1,1)=cos(spherical_coord(2))*sin(spherical_coord(3));
sol(2,1)=sin(spherical_coord(2))*sin(spherical_coord(3));
sol(3,1)=cos(spherical_coord(3));
sol(1,2)=-sin(spherical_coord(2))/(spherical_coord(1)*sin(spherical_coord(3)));
sol(2,2)=cos(spherical_coord(2))/(spherical_coord(1)*sin(spherical_coord(3)));
sol(1,3)=cos(spherical_coord(2))*cos(spherical_coord(3))/spherical_coord(1);
sol(2,3)=sin(spherical_coord(2))*cos(spherical_coord(3))/spherical_coord(1);
sol(3,3)=-sin(spherical_coord(3))/spherical_coord(1);
end

% function sol=rotation_matrix(spherical_coord)
% sol=zeros(3,3);
% sol(1,1)=cos(spherical_coord(2))*sin(spherical_coord(3));
% sol(1,2)=-spherical_coord(1)*sin(spherical_coord(2))*sin(spherical_coord(3));
% sol(1,3)=spherical_coord(1)*cos(spherical_coord(2))*cos(spherical_coord(3));
% sol(2,1)=sin(spherical_coord(2))*sin(spherical_coord(3));
% sol(2,2)=spherical_coord(1)*cos(spherical_coord(2))*sin(spherical_coord(3));
% sol(3,3)=spherical_coord(1)*sin(spherical_coord(2))*cos(spherical_coord(3));
% sol(3,1)=-cos(spherical_coord(3));
% sol(3,3)=-spherical_coord(1)*sin(spherical_coord(3));
% end

% function sol=rotation_matrix(spherical_coord)
% sol=zeros(3,3);
% sol(1,1)=cos(spherical_coord(2))*sin(spherical_coord(3));
% sol(1,2)=-spherical_coord(1)*sin(spherical_coord(2))*sin(spherical_coord(3));
% sol(1,3)=spherical_coord(1)*cos(spherical_coord(2))*cos(spherical_coord(3));
% sol(2,1)=sin(spherical_coord(2))*sin(spherical_coord(3));
% sol(2,2)=spherical_coord(1)*cos(spherical_coord(2))*sin(spherical_coord(3));
% sol(3,3)=spherical_coord(1)*sin(spherical_coord(2))*cos(spherical_coord(3));
% sol(3,1)=-cos(spherical_coord(3));
% sol(3,3)=-spherical_coord(1)*sin(spherical_coord(3));
% end


% function sol=inversely_rotate_sph(rotation_mx)
% size_=size(rotation_mx);
% size_=size_(1);
% sol=zeros(size_,3,3);
% parfor node=1:1:size_
%     sol(node,:,:)=inv(dim_reducer1_2d(rotation_mx(node,:,:)));
% end
% end
% 
% function sol=dim_reducer1_2d(mx)
% sol=zeros(3,3);
% for i=1:1:3
%     for j=1:1:3
%         sol(i,j)=mx(1,i,j);
%     end
% end
% end