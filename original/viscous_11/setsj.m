function sj=setsj(sj,value,nfour,nx,ny)

for i=1:1:nfour
    sj=set_specific_row_sj(sj,value,i,nx,ny);
end

end