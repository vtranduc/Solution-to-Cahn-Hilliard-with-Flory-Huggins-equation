function main_isoparametric_mapping_3d

clear
clc

% === Start the timer =====================================================

counter_init=tic;

% === User specification ==================================================

shape=10;
% So far, we have the following available
% 1 - Cube
% 2- Sphere

% If you wanna add more shapes later, be sure to modify the following
% functions

% main_curved_surface_3d (Just the comments)
% retrieve_isoparametrically_mapped_weights_3d
% create_temperature_gradient_3d

nex=10;
ney=10;
nez=10;

parallel_computing=1;

sig_dig=17;
% 17 is the normal maximum in matlab at the time this script was written.

% === Set Up ==============================================================

nx=nex+1;
ny=ney+1;
nz=nez+1;
n=nx*ny*nz;
ne=nex*ney*nez;

if shape==2
    nx_=nx-2;
    ny_=ny-2;
    nxny=nx*ny;
    nexney=nex*ney;
    [n_nSurf,n_eSurf]=get_n_nSurf_n_eSurf(nex,ney,nez,ny,nz,nx_,ny_);
    neTotal=ne+n_eSurf(6);
elseif shape==10
    nexney=nex*ney;
end

% === Validation ==========================================================

% if shape==2
%     if parallel_computing~=1
%         error('Only parallel computing is available for shape 3')
%     end
% end

% === Compute the coordinates =============================================

fprintf('Generating the mesh\n')
if shape==1
    [X,Y,Z]=generate_cube_mesh_3d(nx,ny,nz,n);
elseif shape==123456789
    [X,Y,Z]=generate_sphere_mesh_3d(nx,ny,nz,n);
%     [X,Y,Z]=generate_sphere_mesh_with_ratio(nx,ny,nz,n);
elseif shape==2
    [X,Y,Z]=generate_sphere_mesh_with_extrusion_3d(nx,ny,nz,n,n_eSurf,n_nSurf,ne);
elseif shape==10
    [X,Y,Z]=generate_sphere_mesh_3d(nx,ny,nz,n);
end

% === Compute the weights =================================================

fprintf('Computing the weights\n')
if shape==1 || shape==123456789
    weights=isoparametrically_adjust_weights_3d(...
        ne,nex,nx,ny,nex*ney,X,Y,Z,parallel_computing);
elseif shape==2
    
    
%     %----------------------------
    disp('Start on surface')
    weights_surface=isoparametrically_adjust_weights_surface_sph(...
        ne,nex,ney,n,nx,ny,nz,nxny,nexney,X,Y,Z,parallel_computing,...
        n_eSurf,n_nSurf,nx_);
    disp('FIN')
%     %-------------------------------------
    
%     return
    
    weights=isoparametrically_adjust_weights_sph(...
        ne,nex,ney,n,nx,ny,nz,nxny,nexney,X,Y,Z,parallel_computing,...
        n_eSurf,n_nSurf,nx_,neTotal);
    
elseif shape==10
    
    weights=isoparametrically_adjust_weights_serendipity(...
        ne,nex,nx,ny,nexney,X,Y,Z,parallel_computing);
    
    %----------------------------
    
    
    
    %----------------------------
    
end

%=----------------------------------------------------------
% take=size(weights);
% ne*8*4*7*3*3*3;
% nnz(weights);
% 
% e=3;
% 
% [eX,eY,eZ]=get_eXYZ_3d(e,nex,nx,ny,nexney,X,Y,Z);
% gps=[(-sqrt(3/5)+1)/2 0.5 (sqrt(3/5)+1)/2];
% phis=get_phis();
% take=compute_elemental_weights_serendipity(gps,phis,eX,eY,eZ);
% 
% take
% 
% plot3(X,Y,Z,'o')
% 
% hold on
% for i=1:1:8
%     plot3(eX(i),eY(i),eZ(i),'rx')
% end
% hold off
% grid on
% xlabel('x');ylabel('y');zlabel('z');
% 
% 
% % return
% 
% 
% for e=1:1:3
%     for orientation=1:1:8
%         for type=1:1:4
%             for order=1:1:7
%                 for ix=1:1:3
%                     for iy=1:1:3
%                         for iz=1:1:3
%                             if weights(e,orientation,type,order,ix,iy,iz)==0 && order~=5
%                                 [e,orientation,type,order,ix,iy,iz]
%                                 warning('just stop')
%                             end
%                         end
%                     end
%                 end
%             end
%         end
%     end
% end
% 
% size(take)
% 
% summation=0;
% coeffs=get_Hermitian_pol_coeffs_3d();
% c=rand(1,32)
% conc_=compute_conc_serendipity(c,take);
% for ix=1:1:3
%     for iy=1:1:3
%         for iz=1:1:3
%             take=compute_specific_c_3d(c,gps(ix),gps(iy),gps(iz),eX,eY,eZ,coeffs);
%             for order=1:1:7
%                 if take(order)==0 %|| conc_(order,ix,iy,iz)==0
%                     [ix,iy,iz,order]
%                     error('dfadgagdag')
%                 end
%                 diff=abs(take(order)-conc_(order,ix,iy,iz));
%                 summation=summation+diff;
%             end
%         end
%     end
% end
% 
% summation
% 
% return
%=----------------------------------------------------------

toc(counter_init)

% === Reshape weights for storage =========================================

% nnz(weights)
% size(weights)

fprintf('Reshaping weights\n')
if shape==10
    weights=compact_weight_for_storage_serendipity(weights,parallel_computing);
else
    weights=compact_weight_for_storage_3d(weights,parallel_computing);
end

toc(counter_init)

% nnz(weights)
% size(weights)
% return

% === Obtain data =========================================================

fprintf('Exporting the data\n')
if ~exist([pwd '/Meshes/'],'dir')
    mkdir([pwd '/Meshes/'])
end

if shape==1
    x_title=sprintf('/Meshes/mesh_cube_X_%i_%i_%i.dat',nex,ney,nez);
    dlmwrite([pwd x_title],X,'precision',sig_dig);
    y_title=sprintf('/Meshes/mesh_cube_Y_%i_%i_%i.dat',nex,ney,nez);
    dlmwrite([pwd y_title],Y,'precision',sig_dig);
    z_title=sprintf('/Meshes/mesh_cube_Z_%i_%i_%i.dat',nex,ney,nez);
    dlmwrite([pwd z_title],Z,'precision',sig_dig);
    weight_title=sprintf('/Meshes/mesh_cube_weights_%i_%i_%i.dat',nex,ney,nez);
    dlmwrite([pwd weight_title],weights,'precision',sig_dig);
elseif shape==123456789
    x_title=sprintf('/Meshes/mesh_sphere_X_%i_%i_%i.dat',nex,ney,nez);
    dlmwrite([pwd x_title],X,'precision',sig_dig);
    y_title=sprintf('/Meshes/mesh_sphere_Y_%i_%i_%i.dat',nex,ney,nez);
    dlmwrite([pwd y_title],Y,'precision',sig_dig);
    z_title=sprintf('/Meshes/mesh_sphere_Z_%i_%i_%i.dat',nex,ney,nez);
    dlmwrite([pwd z_title],Z,'precision',sig_dig);
    weight_title=sprintf('/Meshes/mesh_sphere_weights_%i_%i_%i.dat',nex,ney,nez);
    dlmwrite([pwd weight_title],weights,'precision',sig_dig);
elseif shape==2
    x_title=sprintf('/Meshes/mesh_sphere_with_extrusion_X_%i_%i_%i.dat',nex,ney,nez);
    dlmwrite([pwd x_title],X,'precision',sig_dig);
    y_title=sprintf('/Meshes/mesh_sphere_with_extrusion_Y_%i_%i_%i.dat',nex,ney,nez);
    dlmwrite([pwd y_title],Y,'precision',sig_dig);
    z_title=sprintf('/Meshes/mesh_sphere_with_extrusion_Z_%i_%i_%i.dat',nex,ney,nez);
    dlmwrite([pwd z_title],Z,'precision',sig_dig);
    weight_title=sprintf('/Meshes/mesh_sphere_with_extrusion_weights_%i_%i_%i.dat',nex,ney,nez);
    dlmwrite([pwd weight_title],weights,'precision',sig_dig);
    
    
    %---- Testing ----------------------------
    
    weights_surface=compact_weight_surface_for_storage_3d(weights_surface,parallel_computing);
    dlmwrite('test_weights_surface.dat',weights_surface,'precision',sig_dig)
    
    %------------------------------
elseif shape==10
    x_title=sprintf('/Meshes/mesh_sphere_serendipity_X_%i_%i_%i.dat',nex,ney,nez);
    dlmwrite([pwd x_title],X,'precision',sig_dig);
    y_title=sprintf('/Meshes/mesh_sphere_serendipity_Y_%i_%i_%i.dat',nex,ney,nez);
    dlmwrite([pwd y_title],Y,'precision',sig_dig);
    z_title=sprintf('/Meshes/mesh_sphere_serendipity_Z_%i_%i_%i.dat',nex,ney,nez);
    dlmwrite([pwd z_title],Z,'precision',sig_dig);
    weight_title=sprintf('/Meshes/mesh_sphere_serendipity_weights_%i_%i_%i.dat',nex,ney,nez);
    dlmwrite([pwd weight_title],weights,'precision',sig_dig);
    
    %---- Testing ----------------------------
    
    
    %------------------------------
    
end

% === Conclusion ==========================================================

fprintf('Computation has finished\n')
toc(counter_init)

end