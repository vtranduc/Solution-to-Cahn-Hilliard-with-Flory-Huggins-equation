function [coef_n2,n2_new]=vary_n2(type,spec,n2,time,...
    n2_distribution_type,ne,n2_distribution_spec,nex,ney,xlen,ylen)

%type==2 both ends go from n2(1) to n2(2) at different linear rates

if length(n2)~=2
    error('the length of n2 must be 2')
end

if type==2
    
    if length(spec)~=2
        error('Spec=[time1, time2], where time is the time it takes to reach n2(2) at each end')
    end
    
    n2_variation=(n2(2)-n2(1))*time;
    n2_new=zeros(1,2);
    for i=1:1:2
        n2_new(i)=n2_variation/spec(i)+n2(1);
        if n2_new(i)>n2(2)
            n2_new(i)=n2(2);
        end
    end
    
    coef_n2=n2_characterization(n2_distribution_type,ne,n2_new,...
        n2_distribution_spec,nex,ney,xlen,ylen);
    
elseif type==3
    
    if length(spec)~=1
        error('spec=time to reach n(2) at one end')
    end
    
    n2_new=[n2(1) (n2(2)-n2(1))*time/spec+n2(1)];
%     n2_new=zeros(1,2);
%     for i=1:1:2
%         n2_new(i)=n2_variation/spec(i)+n2(1);
%         if n2_new(i)>n2(2)
%             n2_new(i)=n2(2);
%         end
%     end
    
    coef_n2=n2_characterization(n2_distribution_type,ne,n2_new,...
        n2_distribution_spec,nex,ney,xlen,ylen);
    
else
    error('Specified type is not available')
end


end