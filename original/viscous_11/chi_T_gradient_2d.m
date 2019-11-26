function [Two_chi_n1,coef_T]=chi_T_gradient_2d(T,dx,entropy,n1)
gp=[(-sqrt(3/5)+1)/2 0.5 (sqrt(3/5)+1)/2];
adjusted_gp=gp*dx;
n=length(T)-1;
Two_chi_n1=zeros(n,3);
coef_T=zeros(n,3);
for i=1:1:n
    T_=get_local_T(adjusted_gp,dx,T(i),T(i+1));
    Two_chi_n1(i,:)=get_Two_chi_n1_with_order(T_,entropy,n1,0);
    coef_T(i,:)=T_;
end
end

function solution=get_local_T(adjusted_gp,dx,T_lower,T_higher)
if T_lower==T_higher
    solution=[T_lower T_lower T_lower];
else
    solution=zeros(1,3);
    slope=(T_higher-T_lower)/dx;
    for i=1:1:3
        solution(i)=slope*adjusted_gp(i)+T_lower;
    end
end
end