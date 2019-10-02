function mx=mx_polynomial_order_constraint(degree,x1,x2)
if mod(degree,2)~=1
    error('Degree must be an odd number')
end
degree_plus=degree+1;
mx=zeros(degree_plus,degree_plus);
for order=1:1:degree_plus/2
    for column=1:1:degree_plus
        deg=degree_plus-column-order+1;
        if order==1
            coef=1;
        else
            coef=factorial_sp(degree_plus-column,order-1);
        end
        pre_index=2*order-2;
        if deg>0
            mx(pre_index+1,column)=coef*x1^deg;
            mx(pre_index+2,column)=coef*x2^deg;
        elseif deg==0
            mx(pre_index+1,column)=coef;
            mx(pre_index+2,column)=coef;
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