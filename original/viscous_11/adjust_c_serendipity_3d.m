function [c_adjusted,range_check]=adjust_c_serendipity_3d(c,ci_ave,nx,ny,nz,nfour)
result=extract_abs_result_serendipity(c,nx,ny,nz);
deviation=mean(mean(mean(result)))-ci_ave;
if max(max(max(result)))-deviation>=1 || min(min(min(result)))-deviation<=0
    range_check=0;
    c_adjusted=NaN;
    return
else
    range_check=1;
end
c_adjusted=c;
for i=1:4:nfour-3
    c_adjusted(i)=c_adjusted(i)-deviation;
end

end