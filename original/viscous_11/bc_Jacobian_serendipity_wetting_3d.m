function sj=bc_Jacobian_serendipity_wetting_3d(sj_,g_bc)
sj=sj_;
for i=1:1:length(g_bc)
    sj(g_bc(1,i))=g_bc(2,i);
end
end