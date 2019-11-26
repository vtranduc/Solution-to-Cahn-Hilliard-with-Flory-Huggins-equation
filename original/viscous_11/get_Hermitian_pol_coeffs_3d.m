function coeffs=get_Hermitian_pol_coeffs_3d()
mx=zeros(32,32);
index=0;
ws=...
    [0 0 0
    1 0 0
    0 1 0
    1 1 0
    0 0 1
    1 0 1
    0 1 1
    1 1 1];
for orientation=1:1:8
    for order=1:1:4
        index=index+1;
        mx(index,:)=Hermitian_polynomial_terms_3d(...
            ws(orientation,1),ws(orientation,2),ws(orientation,3),order);
    end
end
coeffs=zeros(32,32);
index=0;
b=zeros(32,1);
for w=1:1:8
    for type=1:1:4
        index=index+1;
        b(index)=1;
        x=mx\b;
        coeffs(:,index)=x;
        b(index)=0;
    end
end
end