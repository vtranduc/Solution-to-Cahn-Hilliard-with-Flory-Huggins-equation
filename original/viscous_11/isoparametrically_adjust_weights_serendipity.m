function weights=isoparametrically_adjust_weights_serendipity(...
    ne,nex,nx,ny,nexney,X,Y,Z,parallel_computing)
weights=zeros(ne,8,4,7,3,3,3);
gps=[(-sqrt(3/5)+1)/2 0.5 (sqrt(3/5)+1)/2];
phis=get_phis();
if parallel_computing==1
    parfor e=1:1:ne
        [eX,eY,eZ]=get_eXYZ_3d(e,nex,nx,ny,nexney,X,Y,Z);
        weights(e,:,:,:,:,:,:)=...
            compute_elemental_weights_serendipity(gps,phis,eX,eY,eZ);
    end
elseif parallel_computing==0
    for e=1:1:ne
        [eX,eY,eZ]=get_eXYZ_3d(e,nex,nx,ny,nexney,X,Y,Z);
        weights(e,:,:,:,:,:,:)=...
            compute_elemental_weights_serendipity(gps,phis,eX,eY,eZ);
    end
end
end

