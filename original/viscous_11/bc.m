function c=bc(c,ll,li,lu,bl,bi,bu,rl,ri,ru,tl,ti,tu)
%Apply initial conditions to concentration

for i=ll:li:lu
    c(i)=0.0;
    c(i+2)=0.0;
end
for i=bl:bi:bu
    c(i)=0.0;
    c(i+1)=0.0;
end
for i=rl:ri:ru
    c(i)=0.0;
    c(i+2)=0.0;
end
for i=tl:ti:tu
    c(i)=0.0;
    c(i+1)=0.0;
end

end