function [sf,sj]=bcsfsj(sf,sj,nx,ny,special_zero,ll,li,lu,bl,bi,bu,rl,ri,ru,tl,ti,tu)
%Apply initial conditions to concentration

for i=ll:li:lu
    sf(i)=0.0;
    sj=set_specific_row_sj(sj,special_zero,i,nx,ny);
    sj(i,i)=1.0;
    sf(i+2)=0.0;
    sj=set_specific_row_sj(sj,special_zero,i+2,nx,ny);
    sj(i+2,i+2)=1.0;
end
for i=bl:bi:bu
    sf(i)=0.0;
    sj=set_specific_row_sj(sj,special_zero,i,nx,ny);
    sj(i,i)=1.0;
    sf(i+1)=0.0;
    sj=set_specific_row_sj(sj,special_zero,i+1,nx,ny);
    sj(i+1,i+1)=1.0;
end
for i=rl:ri:ru
    sf(i)=0.0;
    sj=set_specific_row_sj(sj,special_zero,i,nx,ny);
    sj(i,i)=1.0;
    sf(i+2)=0.0;
    sj=set_specific_row_sj(sj,special_zero,i+2,nx,ny);
    sj(i+2,i+2)=1.0;
end
for i=tl:ti:tu
    sf(i)=0.0;
    sj=set_specific_row_sj(sj,special_zero,i,nx,ny);
    sj(i,i)=1.0;
    sf(i+1)=0.0;
    sj=set_specific_row_sj(sj,special_zero,i+1,nx,ny);
    sj(i+1,i+1)=1.0;
end

end