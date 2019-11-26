function gbs_nodal=get_gbs_serendipity(node)
pre_gbs=(node-1)*4;
gbs_nodal=pre_gbs+1:1:pre_gbs+4;
end