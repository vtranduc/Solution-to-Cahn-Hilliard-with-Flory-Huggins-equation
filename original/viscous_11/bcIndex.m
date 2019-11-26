function [iZeros,iOnes]=bcIndex(irow,icol,ll,li,lu,bl,bi,bu,rl,ri,ru,tl,ti,tu,nx,ny)

z=4*(nx+ny);
iZeros=zeros(z,2);
iOnes=zeros(1,z);

index=0;

for i=ll:li:lu
    identifyZeros=find(irow==i);
    for j=identifyZeros
        if icol(j)==i
            identifyOne=j;
            break
        end
    end
    index=index+1;
    iZeros(index,:)=[identifyZeros(1) identifyZeros(length(identifyZeros))];
    iOnes(index)=identifyOne;
    
    identifyZeros=find(irow==i+2);
    for j=identifyZeros
        if icol(j)==i+2
            identifyOne=j;
            break
        end
    end
    index=index+1;
    iZeros(index,:)=[identifyZeros(1) identifyZeros(length(identifyZeros))];
    iOnes(index)=identifyOne;
end

for i=bl:bi:bu
    identifyZeros=find(irow==i);
    for j=identifyZeros
        if icol(j)==i
            identifyOne=j;
            break
        end
    end
    index=index+1;
    iZeros(index,:)=[identifyZeros(1) identifyZeros(length(identifyZeros))];
    iOnes(index)=identifyOne;
    
    identifyZeros=find(irow==i+1);
    for j=identifyZeros
        if icol(j)==i+1
            identifyOne=j;
            break
        end
    end
    index=index+1;
    iZeros(index,:)=[identifyZeros(1) identifyZeros(length(identifyZeros))];
    iOnes(index)=identifyOne;
end

for i=rl:ri:ru
    identifyZeros=find(irow==i);
    for j=identifyZeros
        if icol(j)==i
            identifyOne=j;
            break
        end
    end
    index=index+1;
    iZeros(index,:)=[identifyZeros(1) identifyZeros(length(identifyZeros))];
    iOnes(index)=identifyOne;
    
    identifyZeros=find(irow==i+2);
    for j=identifyZeros
        if icol(j)==i+2
            identifyOne=j;
            break
        end
    end
    index=index+1;
    iZeros(index,:)=[identifyZeros(1) identifyZeros(length(identifyZeros))];
    iOnes(index)=identifyOne;
end

for i=tl:ti:tu
    identifyZeros=find(irow==i);
    for j=identifyZeros
        if icol(j)==i
            identifyOne=j;
            break
        end
    end
    index=index+1;
    iZeros(index,:)=[identifyZeros(1) identifyZeros(length(identifyZeros))];
    iOnes(index)=identifyOne;
    
    identifyZeros=find(irow==i+1);
    for j=identifyZeros
        if icol(j)==i+1
            identifyOne=j;
            break
        end
    end
    index=index+1;
    iZeros(index,:)=[identifyZeros(1) identifyZeros(length(identifyZeros))];
    iOnes(index)=identifyOne;
end

end