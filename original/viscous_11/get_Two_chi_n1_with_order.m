function solution=get_Two_chi_n1_with_order(T,entropy,n1,order)
n=length(T);
solution=zeros(1,n);
if order==0
    for i=1:1:n
        solution(i)=0.5-entropy*(1-1/T(i));
    end
elseif order==1
    for i=1:1:n
        solution(i)=-entropy/(T(i)^2);
    end
elseif order==2
    for i=1:1:n
        solution(i)=2*entropy/(T(i)^3);
    end
else
    error('Only order of up to 2 is available!')
end
solution=solution*2/n1;
end