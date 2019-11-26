function mx=mx_polynomial_order_constraint_v2(degree,x1,x2,x_mid)

if mod(degree,2)~=0 || degree<2
    error('Degree must be an even number larger or equal to 2')
end
if x_mid<=x1 || x_mid>=x2
    error('The order must be x1<x_mid<x2')
end

degree_plus=degree+1;
mx=zeros(degree_plus,degree_plus);

xs=[x1,x2,x_mid];

for row=1:1:3
    for column=1:1:degree_plus
        deg=degree_plus-column;
        mx(row,column)=xs(row)^deg;
    end
end

for order=2:1:degree/2
    pre_index=3+(order-2)*2;
    for column=1:1:degree_plus
        deg=degree_plus-column-order+1;
        
        if deg>=0
            coef=factorial_sp(degree_plus-column,order-1);
            mx(pre_index+1,column)=coef*x1^deg;
            mx(pre_index+2,column)=coef*x2^deg;
        end
        
        
    end
end

end

function sol=factorial_sp(val,n)
sol=val;
for i=1:1:n-1
    sol=sol*(val-i);
end
end