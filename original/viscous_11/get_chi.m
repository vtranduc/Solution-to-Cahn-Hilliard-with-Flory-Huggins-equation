function solution=get_chi(T_,entropy,T_theta)
n=length(T_);
solution=zeros(1,n);
for i=1:1:n
    T=T_(i);
    solution(i)=0.5-entropy*(1-T_theta/T);
end
end