function gbs_nodal=get_gbs_3d(node)
pre_gbs=(node-1)*8;
gbs_nodal=pre_gbs+1:1:pre_gbs+8;
end