function [sf,sj]=parallel_bc(sj_,iZeros,iOnes,...
    sf_,ll,li,lu,bl,bi,bu,rl,ri,ru,tl,ti,tu)
sj=sj_;
for i=iZeros(1,1):iZeros(1,2):iZeros(1,3)
    sj(i)=0.0;
end
for i=iZeros(2,1):iZeros(2,2):iZeros(2,3)
    for j=0:1:23
        sj(i+j)=0.0;
    end
end
for i=iZeros(3,1):iZeros(3,2):iZeros(3,3)
    for j=0:1:23
        sj(i+j)=0.0;
    end
end
for i=iZeros(4,1):iZeros(4,2):iZeros(4,3)
    sj(i)=0.0;
end
for i=iZeros(5,1):iZeros(5,2):iZeros(5,3)
    for j=0:1:47
        sj(i+j)=0.0;
    end
end
for i=iZeros(6,1):iZeros(6,2):iZeros(6,3)
    for j=0:1:47
        sj(i+j)=0.0;
    end
end
for i=iZeros(7,1):iZeros(7,2):iZeros(7,3)
    sj(i)=0.0;
end
for i=iZeros(8,1):iZeros(8,2):iZeros(8,3)
    for j=0:1:23
        sj(i+j)=0.0;
    end
end
for i=iZeros(9,1):iZeros(9,2):iZeros(9,3)
    for j=0:1:23
        sj(i+j)=0.0;
    end
end
for i=iZeros(10,1):iZeros(10,2):iZeros(10,3)
    sj(i)=0.0;
end
for i=iOnes
    sj(i)=1.0;
end
%---Switch on sf-------
sf=sf_;

for i=ll:li:lu
    sf(i)=0.0;
    sf(i+2)=0.0;
end
for i=bl:bi:bu
    sf(i)=0.0;
    sf(i+1)=0.0;
end
for i=rl:ri:ru
    sf(i)=0.0;
    sf(i+2)=0.0;
end
for i=tl:ti:tu
    sf(i)=0.0;
    sf(i+1)=0.0;
end

end