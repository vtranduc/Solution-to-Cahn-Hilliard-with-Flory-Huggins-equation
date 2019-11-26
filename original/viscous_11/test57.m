function test57

clear
clc

gps=[(-sqrt(3/5)+1)/2 0.5 (sqrt(3/5)+1)/2];
gps_plus=zeros(3,3);
for i=1:1:3
    gps_plus(i,:)=gps*gps(i);
end
w=[5/18 4/9 5/18];
orientation_list=...
    [0 0 0
    1 0 0
    0 1 0
    1 1 0
    0 0 1
    1 0 1
    0 1 1
    1 1 1];

sol=zeros(8,8,7,3,3);

eX=[0 1 0 1 0 1 0 1];
eY=[0 0 1 1 0 0 1 1];
eZ=[0 0 0 0 1 1 1 1];

% sol=type2_transformer_surface_sph(gps,gps_plus,w,eX,eY,eZ,orientation_list,sol);

nex=5;
ney=4;
nez=3;

nx=nex+1;
ny=ney+1;
nz=nez+1;
n=nx*ny*nz;
ne=nex*ney*nez;

if 1==1
    nx_=nx-2;
    ny_=ny-2;
    nxny=nx*ny;
    nexney=nex*ney;
    [n_nSurf,n_eSurf]=get_n_nSurf_n_eSurf(nex,ney,nez,ny,nz,nx_,ny_);
    neTotal=ne+n_eSurf(6);
end

% nodes=get_nodes_of_elements_sph(e,ne,nex,ney,nexney,nx,ny,nz,n,nxny,n_eSurf,n_nSurf,nx_);

% weights=dlmread([pwd '/Meshes/mesh_sphere_with_extrusion_weights_5_4_3.dat']);

X=dlmread([pwd '/Meshes/mesh_sphere_with_extrusion_X_5_4_3.dat']);
Y=dlmread([pwd '/Meshes/mesh_sphere_with_extrusion_Y_5_4_3.dat']);
Z=dlmread([pwd '/Meshes/mesh_sphere_with_extrusion_Z_5_4_3.dat']);

e=15;

% weights_test=zeros(8,8,7,3,3,3);
% 
% size(weights)
% 
% weights=reformat_stored_weights_3d(weights,1);
% 
% for iorientation=1:1:8
%     for itype=1:1:8
%         for iorder=1:1:7
%             for ix=1:1:3
%                 for iy=1:1:3
%                     for iz=1:1:3
%                         weights_test(iorientation,itype,iorder,ix,iy,iz)=...
%                             weights(e,iorientation,itype,iorder,ix,iy,iz);
%                     end
%                 end
%             end
%         end
%     end
% end

% weights=weights(e,:);
% dlmwrite('to_be_deleted_test57',weights,'precision',17);

weights=dlmread('to_be_deleted_test57');

weights=reformat_stored_weights_3d(weights,1);
% 
% size(weights)
% 
nodes=get_nodes_of_element_sph(15,ne,nex,ney,nexney,nx,ny,nz,n,nxny,n_eSurf,n_nSurf,nx_);
eX=zeros(1,8);
eY=zeros(1,8);
eZ=zeros(1,8);

for i=1:1:8
    eX(i)=X(nodes(i));
    eY(i)=Y(nodes(i));
    eZ(i)=Z(nodes(i));
end
% 
% 
% orientation=[0 0 0];
% 
% weights(1,1,2,2,1,2,3)
% 
% test=type2_transformer_specific_3d(gps(1),gps(2),gps(3),gps,w,eX,eY,eZ,orientation);
% 
% test=type3_transformer_specific_3d(gps(1),gps(2),gps(3),gps,w,eX,eY,eZ,orientation);
% weights(1,1,3,2,1,2,3)
% test(2)
% 
% test=type4_transformer_specific_3d(gps(1),gps(2),gps(3),gps,w,eX,eY,eZ,orientation);
% weights(1,1,4,2,1,2,3)
% test(2)
% 
% test=type5_transformer_specific_3d(gps(1),gps(2),gps(3),gps,w,eX,eY,eZ,orientation);
% weights(1,1,5,2,1,2,3)
% test(2)
% 
% test=type6_transformer_specific_3d(gps(1),gps(2),gps(3),gps,w,eX,eY,eZ,orientation);
% weights(1,1,6,2,1,2,3)
% test(2)
% 
% test=type7_transformer_specific_3d(gps(1),gps(2),gps(3),gps,w,eX,eY,eZ,orientation);
% weights(1,1,7,2,1,2,3)
% test(2)
% 
% test=type8_transformer_specific_3d(gps(1),gps(2),gps(3),gps,w,eX,eY,eZ,orientation);
% weights(1,1,8,2,1,2,3)
% test(2)
% 
% test=type1_transformer_specific_3d(gps(1),gps(2),gps(3),eX,eY,eZ,orientation);
% weights(1,1,1,2,1,2,3)
% test(2)
% 
% test101=zeros(1,8,8,7,3,3,3);
% 
% % for iorientation=1:1:8
% %     for ix=1:1:3
% %         for iy=1:1:3
% %             for iz=1:1:3
% %                 test=type1_transformer_specific_3d(gps(ix),gps(iy),gps(iz),eX,eY,eZ,orientation_list(iorientation,:));
% %                 for i=1:1:7
% %                     test101(1,iorientation,1,i,ix,iy,iz)=test(i);
% %                 end
% %                 test=type2_transformer_specific_3d(gps(ix),gps(iy),gps(iz),gps,w,eX,eY,eZ,orientation_list(iorientation,:));
% %                 for i=1:1:7
% %                     test101(1,iorientation,2,i,ix,iy,iz)=test(i);
% %                 end
% %                 test=type3_transformer_specific_3d(gps(ix),gps(iy),gps(iz),gps,w,eX,eY,eZ,orientation_list(iorientation,:));
% %                 for i=1:1:7
% %                     test101(1,iorientation,3,i,ix,iy,iz)=test(i);
% %                 end
% %                 test=type4_transformer_specific_3d(gps(ix),gps(iy),gps(iz),gps,w,eX,eY,eZ,orientation_list(iorientation,:));
% %                 for i=1:1:7
% %                     test101(1,iorientation,4,i,ix,iy,iz)=test(i);
% %                 end
% %                 test=type5_transformer_specific_3d(gps(ix),gps(iy),gps(iz),gps,w,eX,eY,eZ,orientation_list(iorientation,:));
% %                 for i=1:1:7
% %                     test101(1,iorientation,5,i,ix,iy,iz)=test(i);
% %                 end
% %                 test=type6_transformer_specific_3d(gps(ix),gps(iy),gps(iz),gps,w,eX,eY,eZ,orientation_list(iorientation,:));
% %                 for i=1:1:7
% %                     test101(1,iorientation,6,i,ix,iy,iz)=test(i);
% %                 end
% %                 test=type7_transformer_specific_3d(gps(ix),gps(iy),gps(iz),gps,w,eX,eY,eZ,orientation_list(iorientation,:));
% %                 for i=1:1:7
% %                     test101(1,iorientation,7,i,ix,iy,iz)=test(i);
% %                 end
% %                 test=type8_transformer_specific_3d(gps(ix),gps(iy),gps(iz),gps,w,eX,eY,eZ,orientation_list(iorientation,:));
% %                 for i=1:1:7
% %                     test101(1,iorientation,8,i,ix,iy,iz)=test(i);
% %                 end
% %             end
% %         end
% %     end
% % end
% 
% disp('==========================================')
% 
% 
% difference=abs(test101-weights);
% 
% max(max(max(max(max(max(max(max(max(max(difference))))))))))
% 
% % for iorientation=1:1:8
% %     for itype=1:1:7
%         

test102=zeros(1,8,8,7,3,3,3);

for ix=1:1:3
    for iy=1:1:3
        for iz=1:1:3
            take=compute_weights_specific_3d(gps(ix),gps(iy),gps(iz),gps,w,eX,eY,eZ,orientation_list);
            for iorientation=1:1:8
                for itype=1:1:8
                    for iorder=1:1:7
                        test102(1,iorientation,itype,iorder,ix,iy,iz)=take(iorientation,itype,iorder);
                    end
                end
            end
        end
    end
end

difference=abs(test102-weights);

max(max(max(max(max(max(max(max(max(max(difference))))))))))

sum(sum(sum(sum(sum(sum(sum(sum(sum(sum(difference))))))))))

end



%=============================================================================
%=============================================================================
%=============================================================================
%=============================================================================
%=============================================================================
%=============================================================================
%=============================================================================
%=============================================================================
%=============================================================================

function sol=type1_transformer_specific_3d(alpha,beta,gamma,eX,eY,eZ,orientation)
sol=zeros(1,7);
sol(1)=...
    local_basis_3d(alpha,beta,gamma,...
    orientation,[0 0 0],[0 0 0]);
df=get_fg_dervs_3d(alpha,beta,gamma,...
    orientation,[0 0 0],eX,eY,eZ,1:1:6);
for idf=1:1:6
    sol(idf+1)=df(idf);
end
end

function sol=type2_transformer_specific_3d(alpha,beta,gamma,gps,w,eX,eY,eZ,orientation)
sol=zeros(1,7);
gps_plus=alpha*gps;
summation1=0;
summation2=0;
summation3=0;
summation4=0;
summation5=0;
for ialpha=1:1:3
    summation1=summation1+w(ialpha)...
        *local_basis_3d(gps_plus(ialpha),beta,gamma,...
        orientation,[1 0 0],[1 0 0]);
    df=get_fg_dervs_3d(gps_plus(ialpha),beta,gamma,...
        orientation,[1 0 0],eX,eY,eZ,[2 3 5 6]);
    summation2=summation2+w(ialpha)*df(1);
    summation3=summation3+w(ialpha)*df(2);
    summation4=summation4+w(ialpha)*df(3);
    summation5=summation5+w(ialpha)*df(4);
end
mult=alpha*map_to_global_1_compt_3d(...
    alpha,beta,gamma,eX,[1 0 0]);
sol(1)=summation1*mult;
sol(3)=summation2*mult;
sol(4)=summation3*mult;
sol(6)=summation4*mult;
sol(7)=summation5*mult;
sol(2)=local_basis_3d(alpha,beta,gamma,...
    orientation,[1 0 0],[1 0 0]);
sol(5)=get_fg_dervs_3d(alpha,beta,gamma,...
    orientation,[1 0 0],eX,eY,eZ,1);
end

function sol=type3_transformer_specific_3d(alpha,beta,gamma,gps,w,eX,eY,eZ,orientation)
sol=zeros(1,7);
gps_plus=beta*gps;
summation1=0;
summation2=0;
summation3=0;
summation4=0;
summation5=0;
for ibeta=1:1:3
    summation1=summation1+w(ibeta)...
        *local_basis_3d(alpha,gps_plus(ibeta),gamma,...
        orientation,[0 1 0],[0 1 0]);
    df=get_fg_dervs_3d(alpha,gps_plus(ibeta),gamma,...
        orientation,[0 1 0],eX,eY,eZ,[1 3 4 6]);
    summation2=summation2+w(ibeta)*df(1);
    summation3=summation3+w(ibeta)*df(2);
    summation4=summation4+w(ibeta)*df(3);
    summation5=summation5+w(ibeta)*df(4);
end
mult=beta*map_to_global_1_compt_3d(...
    alpha,beta,gamma,eY,[0 1 0]);
sol(1)=summation1*mult;
sol(2)=summation2*mult;
sol(4)=summation3*mult;
sol(5)=summation4*mult;
sol(7)=summation5*mult;
sol(3)=local_basis_3d(alpha,beta,gamma,...
    orientation,[0 1 0],[0 1 0]);
sol(6)=get_fg_dervs_3d(alpha,beta,gamma,...
    orientation,[0 1 0],eX,eY,eZ,2);
end

function sol=type4_transformer_specific_3d(alpha,beta,gamma,gps,w,eX,eY,eZ,orientation)
sol=zeros(1,7);
gps_plus=gamma*gps;
summation1=0;
summation2=0;
summation3=0;
summation4=0;
summation5=0;
for igamma=1:1:3
    summation1=summation1+w(igamma)...
        *local_basis_3d(alpha,beta,gps_plus(igamma),...
        orientation,[0 0 1],[0 0 1]);
    df=get_fg_dervs_3d(alpha,beta,gps_plus(igamma),...
        orientation,[0 0 1],eX,eY,eZ,[1 2 4 5]);
    summation2=summation2+w(igamma)*df(1);
    summation3=summation3+w(igamma)*df(2);
    summation4=summation4+w(igamma)*df(3);
    summation5=summation5+w(igamma)*df(4);
end
mult=gamma*map_to_global_1_compt_3d(...
    alpha,beta,gamma,eZ,[0 0 1]);
sol(1)=summation1*mult;
sol(2)=summation2*mult;
sol(3)=summation3*mult;
sol(5)=summation4*mult;
sol(6)=summation5*mult;
sol(4)=local_basis_3d(alpha,beta,gamma,...
    orientation,[0 0 1],[0 0 1]);
sol(7)=get_fg_dervs_3d(alpha,beta,gamma,...
    orientation,[0 0 1],eX,eY,eZ,3);
end

function sol=type5_transformer_specific_3d(alpha,beta,gamma,gps,w,eX,eY,eZ,orientation)
sol=zeros(1,7);
gps_plus_alpha=alpha*gps;
gps_plus_beta=beta*gps;
summation=0;
for ibeta=1:1:3
    summation=summation+w(ibeta)...
        *local_basis_3d(alpha,gps_plus_beta(ibeta),gamma,...
        orientation,[1 1 0],[1 1 0]);
end
sol(2)=summation*beta...
    *map_to_global_1_compt_3d(alpha,beta,gamma,eY,[0 1 0]);
summation=0;
for ialpha=1:1:3
    summation=summation+w(ialpha)...
        *local_basis_3d(gps_plus_alpha(ialpha),beta,gamma,...
        orientation,[1 1 0],[1 1 0]);
end
sol(3)=summation*alpha...
    *map_to_global_1_compt_3d(alpha,beta,gamma,eX,[1 0 0]);
summation1=0;
summation2=0;
summation3=0;
for ialpha=1:1:3
    for ibeta=1:1:3
        term1=map_to_global_1_compt_3d(gps_plus_alpha(ialpha),gps_plus_beta(ibeta),gamma,eX,[1 0 0]);
        term2=map_to_global_1_compt_3d(gps_plus_alpha(ialpha),gps_plus_beta(ibeta),gamma,eY,[0 1 0]);
        term3=map_to_global_1_compt_3d(gps_plus_alpha(ialpha),gps_plus_beta(ibeta),gamma,eX,[0 1 0]);
        term4=map_to_global_1_compt_3d(gps_plus_alpha(ialpha),gps_plus_beta(ibeta),gamma,eY,[1 0 0]);
        determinant=term1*term2-term3*term4;
        fz=get_fg_dervs_3d(gps_plus_alpha(ialpha),gps_plus_beta(ibeta),gamma,...
            orientation,[1 1 0],eX,eY,eZ,3);
        summation1=summation1+w(ialpha)*w(ibeta)*determinant*fz;
        fzz=get_fg_dervs_3d(gps_plus_alpha(ialpha),gps_plus_beta(ibeta),gamma,...
            orientation,[1 1 0],eX,eY,eZ,6);
        summation2=summation2+w(ialpha)*w(ibeta)*determinant*fzz;
        summation3=summation3+w(ialpha)*w(ibeta)*determinant...
            *local_basis_3d(gps_plus_alpha(ialpha),gps_plus_beta(ibeta),gamma,...
            orientation,[1 1 0],[1 1 0]);
    end
end
sol(4)=summation1*alpha*beta;
sol(7)=summation2*alpha*beta;
sol(1)=summation3*alpha*beta;
summation=0;
for ibeta=1:1:3
    fx=get_fg_dervs_3d(alpha,gps_plus_beta(ibeta),gamma,...
        orientation,[1 1 0],eX,eY,eZ,1);
    summation=summation+w(ibeta)*fx;
end
sol(5)=summation*beta...
    *map_to_global_1_compt_3d(alpha,beta,gamma,eY,[0 1 0]);
summation=0;
for ialpha=1:1:3
    fy=get_fg_dervs_3d(gps_plus_alpha(ialpha),beta,gamma,...
        orientation,[1 1 0],eX,eY,eZ,2);
    summation=summation+w(ialpha)*fy;
end
sol(6)=summation*alpha...
    *map_to_global_1_compt_3d(alpha,beta,gamma,eX,[1 0 0]);
end

function sol=type6_transformer_specific_3d(alpha,beta,gamma,gps,w,eX,eY,eZ,orientation)
sol=zeros(1,7);
gps_plus_alpha=alpha*gps;
gps_plus_gamma=gamma*gps;
summation1=0;
summation2=0;
summation3=0;
for ialpha=1:1:3
    for igamma=1:1:3
        term1=map_to_global_1_compt_3d(gps_plus_alpha(ialpha),beta,gps_plus_gamma(igamma),eX,[1 0 0]);
        term2=map_to_global_1_compt_3d(gps_plus_alpha(ialpha),beta,gps_plus_gamma(igamma),eZ,[0 0 1]);
        term3=map_to_global_1_compt_3d(gps_plus_alpha(ialpha),beta,gps_plus_gamma(igamma),eX,[0 0 1]);
        term4=map_to_global_1_compt_3d(gps_plus_alpha(ialpha),beta,gps_plus_gamma(igamma),eZ,[1 0 0]);
        determinant=term1*term2-term3*term4;
        fy=get_fg_dervs_3d(gps_plus_alpha(ialpha),beta,gps_plus_gamma(igamma),...
            orientation,[1 0 1],eX,eY,eZ,2);
        summation1=summation1+w(ialpha)*w(igamma)*determinant*fy;
        fyy=get_fg_dervs_3d(gps_plus_alpha(ialpha),beta,gps_plus_gamma(igamma),...
            orientation,[1 0 1],eX,eY,eZ,5);
        summation2=summation2+w(ialpha)*w(igamma)*determinant*fyy;
        summation3=summation3+w(ialpha)*w(igamma)*determinant...
            *local_basis_3d(gps_plus_alpha(ialpha),beta,gps_plus_gamma(igamma),...
            orientation,[1 0 1],[1 0 1]);
    end
end
mult=alpha*gamma;
sol(3)=summation1*mult;
sol(6)=summation2*mult;
sol(1)=summation3*mult;
summation1=0;
summation2=0;
for igamma=1:1:3
    summation1=summation1+w(igamma)...
        *local_basis_3d(alpha,beta,gps_plus_gamma(igamma),...
        orientation,[1 0 1],[1 0 1]);
    fx=get_fg_dervs_3d(alpha,beta,gps_plus_gamma(igamma),...
        orientation,[1 0 1],eX,eY,eZ,1);
    summation2=summation2+w(igamma)*fx;
end
mult=gamma*map_to_global_1_compt_3d(...
    alpha,beta,gamma,eZ,[0 0 1]);
sol(2)=summation1*mult;
sol(5)=summation2*mult;
summation1=0;
summation2=0;
for ialpha=1:1:3
    summation1=summation1+w(ialpha)...
        *local_basis_3d(gps_plus_alpha(ialpha),beta,gamma,...
        orientation,[1 0 1],[1 0 1]);
    fz=get_fg_dervs_3d(gps_plus_alpha(ialpha),beta,gamma,...
        orientation,[1 0 1],eX,eY,eZ,3);
    summation2=summation2+w(ialpha)*fz;
end
mult=alpha*map_to_global_1_compt_3d(...
    alpha,beta,gamma,eX,[1 0 0]);
sol(4)=summation1*mult;
sol(7)=summation2*mult;
end

function sol=type7_transformer_specific_3d(alpha,beta,gamma,gps,w,eX,eY,eZ,orientation)
sol=zeros(1,7);
gps_plus_beta=beta*gps;
gps_plus_gamma=gamma*gps;
summation1=0;
summation2=0;
summation3=0;
for ibeta=1:1:3
    for igamma=1:1:3
        term1=map_to_global_1_compt_3d(alpha,gps_plus_beta(ibeta),gps_plus_gamma(igamma),eY,[0 1 0]);
        term2=map_to_global_1_compt_3d(alpha,gps_plus_beta(ibeta),gps_plus_gamma(igamma),eZ,[0 0 1]);
        term3=map_to_global_1_compt_3d(alpha,gps_plus_beta(ibeta),gps_plus_gamma(igamma),eY,[0 0 1]);
        term4=map_to_global_1_compt_3d(alpha,gps_plus_beta(ibeta),gps_plus_gamma(igamma),eZ,[0 1 0]);
        determinant=term1*term2-term3*term4;
        fx=get_fg_dervs_3d(alpha,gps_plus_beta(ibeta),gps_plus_gamma(igamma),...
            orientation,[0 1 1],eX,eY,eZ,1);
        summation1=summation1+w(ibeta)*w(igamma)*determinant*fx;
        fxx=get_fg_dervs_3d(alpha,gps_plus_beta(ibeta),gps_plus_gamma(igamma),...
            orientation,[0 1 1],eX,eY,eZ,4);
        summation2=summation2+w(ibeta)*w(igamma)*determinant*fxx;
        summation3=summation3+w(ibeta)*w(igamma)*determinant...
            *local_basis_3d(alpha,gps_plus_beta(ibeta),gps_plus_gamma(igamma),...
            orientation,[0 1 1],[0 1 1]);
    end
end
mult=beta*gamma;
sol(2)=summation1*mult;
sol(5)=summation2*mult;
sol(1)=summation3*mult;
summation1=0;
summation2=0;
for igamma=1:1:3
    summation1=summation1+w(igamma)...
        *local_basis_3d(alpha,beta,gps_plus_gamma(igamma),...
        orientation,[0 1 1],[0 1 1]);
    fy=get_fg_dervs_3d(alpha,beta,gps_plus_gamma(igamma),...
        orientation,[0 1 1],eX,eY,eZ,2);
    summation2=summation2+w(igamma)*fy;
end
mult=gamma*map_to_global_1_compt_3d(...
    alpha,beta,gamma,eZ,[0 0 1]);
sol(3)=summation1*mult;
sol(6)=summation2*mult;
summation1=0;
summation2=0;
for ibeta=1:1:3
    summation1=summation1+w(ibeta)...
        *local_basis_3d(alpha,gps_plus_beta(ibeta),gamma,...
        orientation,[0 1 1],[0 1 1]);
    fz=get_fg_dervs_3d(alpha,gps_plus_beta(ibeta),gamma,...
        orientation,[0 1 1],eX,eY,eZ,3);
    summation2=summation2+w(ibeta)*fz;
end
mult=beta*map_to_global_1_compt_3d(...
    alpha,beta,gamma,eY,[0 1 0]);
sol(4)=summation1*mult;
sol(7)=summation2*mult;
end

function sol=type8_transformer_specific_3d(alpha,beta,gamma,gps,w,eX,eY,eZ,orientation)
sol=zeros(1,7);
gps_plus_alpha=alpha*gps;
gps_plus_beta=beta*gps;
gps_plus_gamma=gamma*gps;
summation=0;
J=zeros(3,3);
for ialpha=1:1:3
    for ibeta=1:1:3
        for igamma=1:1:3
            J(1,1)=map_to_global_1_compt_3d(gps_plus_alpha(ialpha),gps_plus_beta(ibeta),gps_plus_gamma(igamma),eX,[1 0 0]);
            J(1,2)=map_to_global_1_compt_3d(gps_plus_alpha(ialpha),gps_plus_beta(ibeta),gps_plus_gamma(igamma),eX,[0 1 0]);
            J(1,3)=map_to_global_1_compt_3d(gps_plus_alpha(ialpha),gps_plus_beta(ibeta),gps_plus_gamma(igamma),eX,[0 0 1]);
            J(2,1)=map_to_global_1_compt_3d(gps_plus_alpha(ialpha),gps_plus_beta(ibeta),gps_plus_gamma(igamma),eY,[1 0 0]);
            J(2,2)=map_to_global_1_compt_3d(gps_plus_alpha(ialpha),gps_plus_beta(ibeta),gps_plus_gamma(igamma),eY,[0 1 0]);
            J(2,3)=map_to_global_1_compt_3d(gps_plus_alpha(ialpha),gps_plus_beta(ibeta),gps_plus_gamma(igamma),eY,[0 0 1]);
            J(3,1)=map_to_global_1_compt_3d(gps_plus_alpha(ialpha),gps_plus_beta(ibeta),gps_plus_gamma(igamma),eZ,[1 0 0]);
            J(3,2)=map_to_global_1_compt_3d(gps_plus_alpha(ialpha),gps_plus_beta(ibeta),gps_plus_gamma(igamma),eZ,[0 1 0]);
            J(3,3)=map_to_global_1_compt_3d(gps_plus_alpha(ialpha),gps_plus_beta(ibeta),gps_plus_gamma(igamma),eZ,[0 0 1]);
            summation=summation+w(ialpha)*w(ibeta)*w(igamma)*det(J)...
                *local_basis_3d(gps_plus_alpha(ialpha),gps_plus_beta(ibeta),gps_plus_gamma(igamma),...
                orientation,[1 1 1],[1 1 1]);
        end
    end
end
sol(1)=summation*alpha*beta*gamma;
summation=0;
for ibeta=1:1:3
    for igamma=1:1:3
        term1=map_to_global_1_compt_3d(alpha,gps_plus_beta(ibeta),gps_plus_gamma(igamma),eY,[0 1 0]);
        term2=map_to_global_1_compt_3d(alpha,gps_plus_beta(ibeta),gps_plus_gamma(igamma),eZ,[0 0 1]);
        term3=map_to_global_1_compt_3d(alpha,gps_plus_beta(ibeta),gps_plus_gamma(igamma),eY,[0 0 1]);
        term4=map_to_global_1_compt_3d(alpha,gps_plus_beta(ibeta),gps_plus_gamma(igamma),eZ,[0 1 0]);
        determinant=term1*term2-term3*term4;
        summation=summation+w(ibeta)*w(igamma)*determinant...
            *local_basis_3d(alpha,gps_plus_beta(ibeta),gps_plus_gamma(igamma),orientation,[1 1 1],[1 1 1]);
    end
end
sol(2)=summation*beta*gamma;
summation=0;
for ialpha=1:1:3
    for igamma=1:1:3
        term1=map_to_global_1_compt_3d(gps_plus_alpha(ialpha),beta,gps_plus_gamma(igamma),eX,[1 0 0]);
        term2=map_to_global_1_compt_3d(gps_plus_alpha(ialpha),beta,gps_plus_gamma(igamma),eZ,[0 0 1]);
        term3=map_to_global_1_compt_3d(gps_plus_alpha(ialpha),beta,gps_plus_gamma(igamma),eX,[0 0 1]);
        term4=map_to_global_1_compt_3d(gps_plus_alpha(ialpha),beta,gps_plus_gamma(igamma),eZ,[1 0 0]);
        determinant=term1*term2-term3*term4;
        summation=summation+w(ialpha)*w(igamma)*determinant...
            *local_basis_3d(gps_plus_alpha(ialpha),beta,gps_plus_gamma(igamma),orientation,[1 1 1],[1 1 1]);
    end
end
sol(3)=summation*alpha*gamma;
summation=0;
for ialpha=1:1:3
    for ibeta=1:1:3
        term1=map_to_global_1_compt_3d(gps_plus_alpha(ialpha),gps_plus_beta(ibeta),gamma,eX,[1 0 0]);
        term2=map_to_global_1_compt_3d(gps_plus_alpha(ialpha),gps_plus_beta(ibeta),gamma,eY,[0 1 0]);
        term3=map_to_global_1_compt_3d(gps_plus_alpha(ialpha),gps_plus_beta(ibeta),gamma,eX,[0 1 0]);
        term4=map_to_global_1_compt_3d(gps_plus_alpha(ialpha),gps_plus_beta(ibeta),gamma,eY,[1 0 0]);
        determinant=term1*term2-term3*term4;
        summation=summation+w(ialpha)*w(ibeta)*determinant...
            *local_basis_3d(gps_plus_alpha(ialpha),gps_plus_beta(ibeta),gamma,orientation,[1 1 1],[1 1 1]);
    end
end
sol(4)=summation*alpha*beta;
summation=0;
for ibeta=1:1:3
    for igamma=1:1:3
        term1=map_to_global_1_compt_3d(alpha,gps_plus_beta(ibeta),gps_plus_gamma(igamma),eY,[0 1 0]);
        term2=map_to_global_1_compt_3d(alpha,gps_plus_beta(ibeta),gps_plus_gamma(igamma),eZ,[0 0 1]);
        term3=map_to_global_1_compt_3d(alpha,gps_plus_beta(ibeta),gps_plus_gamma(igamma),eY,[0 0 1]);
        term4=map_to_global_1_compt_3d(alpha,gps_plus_beta(ibeta),gps_plus_gamma(igamma),eZ,[0 1 0]);
        determinant=term1*term2-term3*term4;
        fx=get_fg_dervs_3d(alpha,gps_plus_beta(ibeta),gps_plus_gamma(igamma),orientation,[1 1 1],eX,eY,eZ,1);
        summation=summation+w(ibeta)*w(igamma)*determinant*fx;
    end
end
sol(5)=summation*beta*gamma;
summation=0;
for ialpha=1:1:3
    for igamma=1:1:3
        term1=map_to_global_1_compt_3d(gps_plus_alpha(ialpha),beta,gps_plus_gamma(igamma),eX,[1 0 0]);
        term2=map_to_global_1_compt_3d(gps_plus_alpha(ialpha),beta,gps_plus_gamma(igamma),eZ,[0 0 1]);
        term3=map_to_global_1_compt_3d(gps_plus_alpha(ialpha),beta,gps_plus_gamma(igamma),eX,[0 0 1]);
        term4=map_to_global_1_compt_3d(gps_plus_alpha(ialpha),beta,gps_plus_gamma(igamma),eZ,[1 0 0]);
        determinant=term1*term2-term3*term4;
        fy=get_fg_dervs_3d(gps_plus_alpha(ialpha),beta,gps_plus_gamma(igamma),orientation,[1 1 1],eX,eY,eZ,2);
        summation=summation+w(ialpha)*w(igamma)*determinant*fy;
    end
end
sol(6)=summation*alpha*gamma;
summation=0;
for ialpha=1:1:3
    for ibeta=1:1:3
        term1=map_to_global_1_compt_3d(gps_plus_alpha(ialpha),gps_plus_beta(ibeta),gamma,eX,[1 0 0]);
        term2=map_to_global_1_compt_3d(gps_plus_alpha(ialpha),gps_plus_beta(ibeta),gamma,eY,[0 1 0]);
        term3=map_to_global_1_compt_3d(gps_plus_alpha(ialpha),gps_plus_beta(ibeta),gamma,eX,[0 1 0]);
        term4=map_to_global_1_compt_3d(gps_plus_alpha(ialpha),gps_plus_beta(ibeta),gamma,eY,[1 0 0]);
        determinant=term1*term2-term3*term4;
        fz=get_fg_dervs_3d(gps_plus_alpha(ialpha),gps_plus_beta(ibeta),gamma,orientation,[1 1 1],eX,eY,eZ,3);
        summation=summation+w(ialpha)*w(ibeta)*determinant*fz;
    end
end
sol(7)=summation*alpha*beta;
end