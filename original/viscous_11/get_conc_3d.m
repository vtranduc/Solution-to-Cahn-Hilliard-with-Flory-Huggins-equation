function solution=get_conc_3d(ne,c,weights,nexney,nex,nx,ny)
solution=zeros(ne,7,3,3,3);
parfor e=1:1:ne %parfor must be feasible!
    solution(e,:,:,:,:)=compute_elemental_weights_3d(e,c,weights,nexney,nex,nx,ny);
end
end