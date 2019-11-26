function sol=get_fg_dervs_3d(alpha,beta,gamma,orientations,types,eX,eY,eZ,terms_to_be_obtained)
invJ=get_invJ(alpha,beta,gamma,eX,eY,eZ);
f_der=get_f_dervs_3d(alpha,beta,gamma,orientations,types);
sol=zeros(1,length(terms_to_be_obtained));
index=0;
for i=terms_to_be_obtained
    val=0;
    for j=1:1:9
        val=val+invJ(i,j)*f_der(j);
    end
    index=index+1;
    sol(index)=val;
end
end