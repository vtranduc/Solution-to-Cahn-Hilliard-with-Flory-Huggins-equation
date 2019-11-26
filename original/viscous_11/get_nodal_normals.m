function sol=get_nodal_normals(n,n_nSurf,X,Y,Z)
sol=zeros(n_nSurf(6),3);
for i=1:1:n_nSurf(6)
    sol(i,:)=convert_to_unit_vector([X(n+i) Y(n+i) Z(n+i)]);
end
end