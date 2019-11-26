function co=icc3D(co,ll,li,lu,bl,bi,bu,rl,ri,ru,tl,ti,tu,zi,zu,nxy)

for i=ll:li:lu
    for j=i:zi:i+zu
        co(j)=0.0;
    end
end
for i=bl:bi:bu
    for j=i:zi:i+zu
        co(j)=0.0;
    end
end
for i=rl:ri:ru
    for j=i:zi:i+zu
        co(j)=0.0;
    end
end
for i=tl:ti:tu
    for j=i:zi:i+zu
        co(j)=0.0;
    end
end

for i=4:8:nxy-4
    co(i)=0.0;
    co(i+zu)=0.0;
end

end