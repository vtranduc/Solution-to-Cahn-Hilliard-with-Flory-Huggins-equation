function sol=compute_weight_specific(alpha,beta,gamma,orientation,type,order,coeffs)
terms=Hermitian_polynomial_terms_3d(alpha,beta,gamma,order);
coeff=coeffs(:,(orientation-1)*4+type);
sol=dot(terms,coeff');
end