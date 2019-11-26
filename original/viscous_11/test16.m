function test16

clear
clc

n=100;
domain=linspace(0,1,n);

plot(domain,myfun(domain),domain,log(ones(1,n)-(domain)))
grid on
legend('Taylor','Analytical')

syms x
for n=0:1:10
    myfun2(n)
end


end

function solution=myfun(domain)
n=length(domain);
solution=zeros(1,n);
for i=1:1:n
    x=domain(i);
    solution(i)=-x-0.5*x^2-(1/3)*x^3-0.25*x^4;
end
end

function solution=myfun2(n)
syms x a
solution=log(1-a);
for i=1:1:n
    solution=diff(solution);
end
solution=solution*(x^n)/factorial(n);
solution=subs(solution,a,0);
end