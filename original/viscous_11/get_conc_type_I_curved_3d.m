function conc_=get_conc_type_I_curved_3d(ne,c,weights,nexney,nex,nx,ny)
conc_=zeros(ne,3,3,3);

parfor e=1:1:ne
    conc_(e,:,:,:)=compute_elemental_weights_type_I_curved_3d(e,c,weights,nexney,nex,nx,ny);
end
end