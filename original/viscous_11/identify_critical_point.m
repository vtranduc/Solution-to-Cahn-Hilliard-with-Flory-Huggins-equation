function [c_critical,T_critical]=identify_critical_point(n1,n2,entropy,T_theta)

gamma=2*n1*n2;
A=gamma^2;
B=-2*gamma*(n1+n2);
C=(n2-n1)^2;
Xcutoff1=(-B+sqrt(B^2-4*A*C))/(2*A);
T_critical=(((Xcutoff1-0.5)/entropy+1)^-1)/T_theta;

chi=get_chi(T_critical,entropy,T_theta);
alpha=2*chi*n1*n2;
beta=n1-n2-alpha;
c_critical=-beta/(2*alpha);

end