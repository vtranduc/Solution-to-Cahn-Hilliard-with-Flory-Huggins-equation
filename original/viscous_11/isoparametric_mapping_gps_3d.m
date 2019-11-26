function isoparametrically_mapped_gps=isoparametric_mapping_gps_3d(ne,nex,nexney,nx,ny,X,Y,Z,parallel_computing)
bases=get_isoparametric_bases_3d();
isoparametrically_mapped_gps=zeros(ne,3,3,3,3);
if parallel_computing==1
    parfor e=1:1:ne
        xyz_mx=get_e_coords_3d(e,nex,nexney,nx,ny,X,Y,Z);
        isoparametrically_mapped_gps(e,:,:,:,:)=get_transformed_weight_points_3d(xyz_mx,bases);
    end
elseif parallel_computing==0
    for e=1:1:ne
        xyz_mx=get_e_coords_3d(e,nex,nexney,nx,ny,X,Y,Z);
        isoparametrically_mapped_gps(e,:,:,:,:)=get_transformed_weight_points_3d(xyz_mx,bases);
    end
end
end

function xyz_mx=get_e_coords_3d(e,nex,nexney,nx,ny,X,Y,Z)
xyz_mx=zeros(8,3);
[n1xth,n1yth,n1zth]=get_n1xyzth_3d(e,nex,nexney);
nodes=zeros(1,8);
nodes(1)=get_node_3d(n1xth,n1yth,n1zth,nx,ny);
nodes(2)=nodes(1)+1;
nodes(3)=get_node_3d(n1xth,n1yth+1,n1zth,nx,ny);
nodes(4)=nodes(3)+1;
nodes(5)=get_node_3d(n1xth,n1yth,n1zth+1,nx,ny);
nodes(6)=nodes(5)+1;
nodes(7)=get_node_3d(n1xth,n1yth+1,n1zth+1,nx,ny);
nodes(8)=nodes(7)+1;
for i=1:1:8
    xyz_mx(i,1)=X(nodes(i));
    xyz_mx(i,2)=Y(nodes(i));
    xyz_mx(i,3)=Z(nodes(i));
end
end

function bases=get_isoparametric_bases_3d()
gps=[(-sqrt(3/5)+1)/2 0.5 (sqrt(3/5)+1)/2];
w_1comp=zeros(2,3);
w_1comp(2,:)=gps;
for i=1:1:3
    w_1comp(1,i)=1-gps(i);
end
bases=zeros(3,3,3,8);
for gammath=1:1:3
    for betath=1:1:3
        for alphath=1:1:3
            bases(alphath,betath,gammath,1)=w_1comp(1,alphath)*w_1comp(1,betath)*w_1comp(1,gammath);
            bases(alphath,betath,gammath,2)=w_1comp(2,alphath)*w_1comp(1,betath)*w_1comp(1,gammath);
            bases(alphath,betath,gammath,3)=w_1comp(1,alphath)*w_1comp(2,betath)*w_1comp(1,gammath);
            bases(alphath,betath,gammath,4)=w_1comp(2,alphath)*w_1comp(2,betath)*w_1comp(1,gammath);
            bases(alphath,betath,gammath,5)=w_1comp(1,alphath)*w_1comp(1,betath)*w_1comp(2,gammath);
            bases(alphath,betath,gammath,6)=w_1comp(2,alphath)*w_1comp(1,betath)*w_1comp(2,gammath);
            bases(alphath,betath,gammath,7)=w_1comp(1,alphath)*w_1comp(2,betath)*w_1comp(2,gammath);
            bases(alphath,betath,gammath,8)=w_1comp(2,alphath)*w_1comp(2,betath)*w_1comp(2,gammath);
        end
    end
end
end