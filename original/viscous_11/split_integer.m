function sol=split_integer(val,nSplits)
z=mod(val,nSplits);
if z==0
    sol=(val/nSplits)*ones(1,nSplits);
else
    sol=(val-z)/nSplits*ones(1,nSplits);
    for i=1:1:z
        sol(i)=sol(i)+1;
    end
end
end