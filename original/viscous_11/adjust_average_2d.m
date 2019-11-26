function c_adjusted=adjust_average_2d(c_unadjusted,adjusted_average,n,nfour)
c_adjusted=c_unadjusted;
summation=0;
for i=1:4:nfour-3
    summation=summation+c_unadjusted(i);
end
deviation=summation/n-adjusted_average;
if deviation==0
    return
else
    for i=1:4:nfour-3
        c_adjusted(i)=c_adjusted(i)-deviation;
    end
end
end