function [adjusted_chi,adjusted_diffT]=chi_diffT_adjuster(T,dx,entropy,T_theta,diff)
gp=[(-sqrt(3/5)+1)/2 0.5 (sqrt(3/5)+1)/2];
adjusted_gp=gp*dx;
n=length(T)-1;
adjusted_chi=zeros(n,3);
adjusted_diffT=zeros(n,3);
for i=1:1:n
    T_lower=T(i);
    T_higher=T(i+1);
    T_=get_local_T(adjusted_gp,dx,T_lower,T_higher);
    adjusted_chi(i,:)=get_chi(T_,entropy,T_theta);
    adjusted_diffT(i,:)=T_*diff;
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