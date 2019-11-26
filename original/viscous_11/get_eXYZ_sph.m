function [eX,eY,eZ]=get_eXYZ_sph(e,ne,nex,ney,nexney,nx,ny,nz,n,nxny,n_eSurf,n_nSurf,nx_,X,Y,Z)
nodes=get_nodes_of_element_sph(e,ne,nex,ney,nexney,nx,ny,nz,n,nxny,n_eSurf,n_nSurf,nx_);
eX=zeros(1,8);
eY=zeros(1,8);
eZ=zeros(1,8);
for i=1:1:8
    eX(i)=X(nodes(i));
    eY(i)=Y(nodes(i));
    eZ(i)=Z(nodes(i));
end
end