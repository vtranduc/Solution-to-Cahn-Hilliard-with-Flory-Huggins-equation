function test68(n_nSurf,rotation_mx,inverse_rotation)

for i=1:1:n_nSurf(6)
    A=dim_reducer1_2d(rotation_mx(i,:,:));
    B=dim_reducer1_2d(inverse_rotation(i,:,:));
    C=inv(A);
    D=B-C;
    if sum(sum(D))>10^-15
        sum(sum(D))
        error('bfdasdfgag')
    end
end

return

for i=1:1:n_nSurf(6)
    A=dim_reducer1_2d(rotation_mx(i,:,:));
    B=dim_reducer1_2d(inverse_rotation(i,:,:));
    C=A*B;
    D=C-[1,0,0;0,1,0;0,0,1];
    if sum(sum(D))>10^-15
        abs(det(C)-1)
        error('bfdasdfgag')
    end
end
    
end

function sol=dim_reducer1_2d(mx)
sol=zeros(3,3);
for i=1:1:3
    for j=1:1:3
        sol(i,j)=mx(1,i,j);
    end
end
end