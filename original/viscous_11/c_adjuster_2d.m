function [sol,range_check]=c_adjuster_2d(c,c_correct,n,nfour)
sol=c;
range_check=1;
summation=0;
for i=1:4:nfour-3
    summation=summation+c(i);
end
deviation=summation/n-c_correct;
if deviation~=0
    for i=1:4:nfour-3
        sol(i)=c(i)-deviation;
        if sol(i)>=1 || sol(i)<=0
            range_check=0;
            return
        end
    end
end
end