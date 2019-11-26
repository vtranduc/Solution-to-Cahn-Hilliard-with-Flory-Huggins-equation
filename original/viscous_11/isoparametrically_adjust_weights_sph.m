function sol=isoparametrically_adjust_weights_sph(ne,nex,ney,n,nx,ny,nz,nxny,nexney,X,Y,Z,parallel_computing,n_eSurf,n_nSurf,nx_,neTotal)

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
sol=zeros(neTotal,8,8,7,3,3,3);
if parallel_computing==1
    parfor e=1:1:neTotal
        [eX,eY,eZ]=get_eXYZ_sph(e,ne,nex,ney,nexney,nx,ny,nz,n,nxny,n_eSurf,n_nSurf,nx_,X,Y,Z);
        sol(e,:,:,:,:,:,:)=transformer_3d(gps,gps_plus,w,eX,eY,eZ,orientation_list);
    end
elseif parallel_computing==0
    for e=1:1:neTotal
        [eX,eY,eZ]=get_eXYZ_sph(e,ne,nex,ney,nexney,nx,ny,nz,n,nxny,n_eSurf,n_nSurf,nx_,X,Y,Z);
        sol(e,:,:,:,:,:,:)=transformer_3d(gps,gps_plus,w,eX,eY,eZ,orientation_list);
    end
end
end

function sol=transformer_3d(gps,gps_plus,w,eX,eY,eZ,orientation_list)
sol=zeros(8,8,7,3,3,3);
sol=type1_transformer_3d(gps,eX,eY,eZ,orientation_list,sol);
sol=type2_transformer_3d(gps,gps_plus,w,eX,eY,eZ,orientation_list,sol);
sol=type3_transformer_3d(gps,gps_plus,w,eX,eY,eZ,orientation_list,sol);
sol=type4_transformer_3d(gps,gps_plus,w,eX,eY,eZ,orientation_list,sol);
sol=type5_transformer_3d(gps,gps_plus,w,eX,eY,eZ,orientation_list,sol);
sol=type6_transformer_3d(gps,gps_plus,w,eX,eY,eZ,orientation_list,sol);
sol=type7_transformer_3d(gps,gps_plus,w,eX,eY,eZ,orientation_list,sol);
sol=type8_transformer_3d(gps,gps_plus,w,eX,eY,eZ,orientation_list,sol);
end