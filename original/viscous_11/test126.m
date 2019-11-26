function test126

clear
clc

x1=0;
x2=1;
x_mid=0.9;

degree=4;

T1=0;
T2=100;
T_mid=(T1+T2)/2;


mx=mx_polynomial_order_constraint_v2(degree,x1,x2,x_mid);

b=zeros(degree+1,1);
b(1)=T1;b(2)=T2;b(3)=T_mid;

size(mx)
size(b)

coef=mx\b;

% mx
% 
% [coef b]

%------Plot--------
n=100;
degree_plus=degree+1;

domain=linspace(x1,x2,n);
range=zeros(1,n);

for i=1:1:n
    summation=0;
    for j=1:1:degree+1
        summation=summation+coef(j)*domain(i)^(degree_plus-j);
    end
    range(i)=summation;
end

% plot(domain,range)
% grid on
% grid minor

%-------Old file test-----------------

% degree=15;
% degree_plus=degree+1;
% mx=mx_polynomial_order_constraint(degree,x1,x2);
% b=zeros(degree_plus,1);
% b(1)=T1;b(2)=T2;
% coef=mx\b;
% 
% for i=1:1:n
%     summation=0;
%     for j=1:1:degree+1
%         summation=summation+coef(j)*domain(i)^(degree_plus-j);
%     end
%     range(i)=summation;
% end
% 
% plot(domain,range)
% grid on
% grid minor

%------------------------- test char

n1=1;
n2=1;
diff=100000;
co=0.7;
T=0.34;

chi=chi_test(T)

n=100;
n2_domain=linspace(1.1,2,n);
f_char=zeros(1,n);

chunk1=-diff*T;
chunk2=2*pi*sqrt(2);
for in2=1:1:n
    n2=n2_domain(in2);
    d2f_dc2=get_d2f_dc2_test(n1,n2,co,T);
    
    f_char(in2)=sqrt(chunk1*d2f_dc2)/chunk2;
end

plot(n2_domain,f_char)

end

function d2f_dc2=get_d2f_dc2_test(n1,n2,co,T)
chi=chi_test(T);
d2f_dc2=1/(co*n1)+1/((1-co)*n2)-2*chi/n1;
end

function chi=chi_test(T)
chi=1/T-0.5;
end