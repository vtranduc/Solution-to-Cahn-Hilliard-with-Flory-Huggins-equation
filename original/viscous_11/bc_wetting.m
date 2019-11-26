function c=bc_wetting(c,ll,li,lu,bl,bi,bu,rl,ri,ru,tl,ti,tu,...
    gl,gb,gr,gt,hl,hb,hr,ht)
%Apply initial conditions to concentration

index=0;

% [ll,li,lu
%     bl,bi,bu
%     rl,ri,ru
%     tl,ti,tu]
% error('mikata')

for i=ll:li:lu
%     c(i)=0.0;
    index=index+1;
    c(i)=hl(index)+gl(index)*c(i-1);
    
%     c(i+2)=0.0;
end
index=0;
for i=bl:bi:bu
%     c(i)=0.0;
    index=index+1;
    c(i)=hb(index)+gb(index)*c(i-2);
    
%     c(i+1)=0.0;
end
index=0;
for i=rl:ri:ru
%     c(i)=0.0;
    index=index+1;
    c(i)=hr(index)+gr(index)*c(i-1);
    
%     c(i+2)=0.0;
end
index=0;
for i=tl:ti:tu
%     c(i)=0.0;
    index=index+1;
    c(i)=ht(index)+gt(index)*c(i-2);
    
%     c(i+1)=0.0;
end

end