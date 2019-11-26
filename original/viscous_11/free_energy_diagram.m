function [cs,f,potential]=free_energy_diagram(n1,n2,T,entropy,T_theta)

kbt_v=1;

n=100;

cs=linspace(0,1,n);

% n1=1;
% n2=10;
% T=0.35;
% 
% entropy=1;
% T_theta=1;

chi=get_chi(T,entropy,T_theta);

f=free_energy(cs,n1,n2,kbt_v,chi);
potential=chemical_potential(cs,n1,n2,kbt_v,chi);

% subplot(1,2,1)
% plot(cs,f)
% grid on
% subplot(1,2,2)
% plot(cs,potential)
% grid on



end

function solution=free_energy(cs,n1,n2,kbt_v,chi)
n=length(cs);
solution=zeros(1,n);
for i=1:1:n
    solution(i)=(cs(i)/n1)*log(cs(i))+((1-cs(i))/n2)*log(1-cs(i))+chi*cs(i)*(1-cs(i));
end
solution=solution*kbt_v;
end

function solution=chemical_potential(cs,n1,n2,kbt_v,chi)
n=length(cs);
solution=zeros(1,n);
for i=1:1:n
    solution(i)=log(cs(i))/n1+1/n1-log(1-cs(i))/n2-1/n2+chi*(1-cs(i))-chi*cs(i);
end
solution=solution*kbt_v;
end