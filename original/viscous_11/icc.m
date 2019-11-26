function co=icc(co,ll,li,lu,bl,bi,bu,rl,ri,ru,tl,ti,tu)
%Apply initial conditions to concentration

for i=ll:li:lu
    co(i)=0.0;
    co(i+2)=0.0;
end
for i=bl:bi:bu
    co(i)=0.0;
    co(i+1)=0.0;
end
for i=rl:ri:ru
    co(i)=0.0;
    co(i+2)=0.0;
end
for i=tl:ti:tu
    co(i)=0.0;
    co(i+1)=0.0;
end

end