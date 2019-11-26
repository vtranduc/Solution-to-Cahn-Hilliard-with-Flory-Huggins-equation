function co=generate_co(ci,ci_fluc,ll,li,lu,bl,bi,bu,rl,ri,ru,tl,ti,tu,zi,zu,bny,neight)

co = zeros(1,neight);

for i=1:neight
    co(i)=ci+(ci_fluc*(2*rand(1,1)-1));
end

co=icc3D(co,ll,li,lu,bl,bi,bu,rl,ri,ru,tl,ti,tu,zi,zu,bny);

end