function thesis_1

clear
clc


%Compute coefficients of polynomials of Hermitians

sols=zeros(32,32);

% Hermitian_polynomial_terms_3d(a,b,c,order);

alphas=[0 1];
betas=[0 1];
gammas=[0 1];

row=0;
for igamma=1:1:2
    for ibeta=1:1:2
        for ialpha=1:1:2
            for order=1:1:4
                row=row+1;
                serendipity_Hermitian=Hermitian_polynomial_terms_3d(alphas(ialpha),betas(ibeta),gammas(igamma),order);
                mx(row,:)=serendipity_Hermitian';
            end
        end
    end
end

for i=1:1:32
    b=zeros(32,1);
    b(i)=1;
    x=mx\b;
    sols(i,:)=x';
end

sols

%Test------------
alpha=1;
beta=1;
gamma=1;

orientation=8;
type=2;

order=2;

val=dot(sols((orientation-1)*4+type,:),Hermitian_polynomial_terms_3d(alpha,beta,gamma,order))

xlswrite('coef',sols)

%------------------


end