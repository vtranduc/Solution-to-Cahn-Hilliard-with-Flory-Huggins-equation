function [X,Y,Z,weights]=...
    retrieve_isoparametrically_mapped_weights_3d(...
    shape,nex,ney,nez,parallel_computing)

if shape==1
    x_title=sprintf('/Meshes/mesh_cube_X_%i_%i_%i.dat',nex,ney,nez);
    X=dlmread([pwd x_title]);
    y_title=sprintf('/Meshes/mesh_cube_Y_%i_%i_%i.dat',nex,ney,nez);
    Y=dlmread([pwd y_title]);
    z_title=sprintf('/Meshes/mesh_cube_Z_%i_%i_%i.dat',nex,ney,nez);
    Z=dlmread([pwd z_title]);
    weight_title=sprintf('/Meshes/mesh_cube_weights_%i_%i_%i.dat',nex,ney,nez);
    weights=dlmread([pwd weight_title]);
    weights=reformat_stored_weights_3d(weights,parallel_computing);
% elseif shape==123456789
%     x_title=sprintf('/Meshes/mesh_sphere_X_%i_%i_%i.dat',nex,ney,nez);
%     X=dlmread([pwd x_title]);
%     y_title=sprintf('/Meshes/mesh_sphere_Y_%i_%i_%i.dat',nex,ney,nez);
%     Y=dlmread([pwd y_title]);
%     z_title=sprintf('/Meshes/mesh_sphere_Z_%i_%i_%i.dat',nex,ney,nez);
%     Z=dlmread([pwd z_title]);
%     weight_title=sprintf('/Meshes/mesh_sphere_weights_%i_%i_%i.dat',nex,ney,nez);
%     weights=dlmread([pwd weight_title]);
%     weights=reformat_stored_weights_3d(weights,parallel_computing);
elseif shape==2
    x_title=sprintf('/Meshes/mesh_sphere_with_extrusion_X_%i_%i_%i.dat',nex,ney,nez);
    X=dlmread([pwd x_title]);
    y_title=sprintf('/Meshes/mesh_sphere_with_extrusion_Y_%i_%i_%i.dat',nex,ney,nez);
    Y=dlmread([pwd y_title]);
    z_title=sprintf('/Meshes/mesh_sphere_with_extrusion_Z_%i_%i_%i.dat',nex,ney,nez);
    Z=dlmread([pwd z_title]);
    weight_title=sprintf('/Meshes/mesh_sphere_with_extrusion_weights_%i_%i_%i.dat',nex,ney,nez);
    weights=dlmread([pwd weight_title]);
    weights=reformat_stored_weights_3d(weights,parallel_computing);
elseif shape==10
    x_title=sprintf('/Meshes/mesh_sphere_serendipity_X_%i_%i_%i.dat',nex,ney,nez);
    X=dlmread([pwd x_title]);
    y_title=sprintf('/Meshes/mesh_sphere_serendipity_Y_%i_%i_%i.dat',nex,ney,nez);
    Y=dlmread([pwd y_title]);
    z_title=sprintf('/Meshes/mesh_sphere_serendipity_Z_%i_%i_%i.dat',nex,ney,nez);
    Z=dlmread([pwd z_title]);
    weight_title=sprintf('/Meshes/mesh_sphere_serendipity_weights_%i_%i_%i.dat',nex,ney,nez);
    weights=dlmread([pwd weight_title]);
    weights=reformat_stored_weights_serendipity(weights,parallel_computing);
else
%     error('Only cube (shape=1), sphere (shape=2), or sphere with extrusion (shape=3) are available now')
    error('Only cube (shape=1) or sphere with extrusion (shape=2), or sphere with serendipity mesh (shape=10) are available now')
end

end