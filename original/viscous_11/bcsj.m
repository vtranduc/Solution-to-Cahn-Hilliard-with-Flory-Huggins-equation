function sj=bcsj(sj,nx,ny,ll,li,lu,bl,bi,bu,rl,ri,ru,tl,ti,tu)

for i=ll:li:lu
    sj=set_specific_row_sj(sj,0,i,nx,ny);
end
for i=bl:bi:bu
    sj=set_specific_row_sj(sj,0,i,nx,ny);
end
for i=rl:ri:ru
    sj=set_specific_row_sj(sj,0,i,nx,ny);
end
for i=tl:ti:tu
    sj=set_specific_row_sj(sj,0,i,nx,ny);
end

end