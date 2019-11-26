function main_derivation

clear
clc

coeffs=get_Hermitian_pol_coeffs_3d();

alpha=1;
beta=0;
gamma=1;
w=6;
t=4;
order=6;

sol=compute_weight_specific(alpha,beta,gamma,w,t,order,coeffs)

end