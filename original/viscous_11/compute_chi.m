function solution=compute_chi(T_,entropy)
n=length(T_);
solution=zeros(1,n);
for i=1:1:n
    T=T_(i);
    solution(i)=0.5-entropy*(1-1/T);
end
end