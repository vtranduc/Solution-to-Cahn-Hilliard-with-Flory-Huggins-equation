function f_der=get_f_dervs_3d(alpha,beta,gamma,orientations,types)
f_der=zeros(1,9);
f_der(1)=local_basis_3d(alpha,beta,gamma,orientations,types,[1 0 0]+types);
f_der(2)=local_basis_3d(alpha,beta,gamma,orientations,types,[0 1 0]+types);
f_der(3)=local_basis_3d(alpha,beta,gamma,orientations,types,[0 0 1]+types);
f_der(4)=local_basis_3d(alpha,beta,gamma,orientations,types,[2 0 0]+types);
f_der(5)=local_basis_3d(alpha,beta,gamma,orientations,types,[0 2 0]+types);
f_der(6)=local_basis_3d(alpha,beta,gamma,orientations,types,[0 0 2]+types);
f_der(7)=local_basis_3d(alpha,beta,gamma,orientations,types,[1 1 0]+types);
f_der(8)=local_basis_3d(alpha,beta,gamma,orientations,types,[1 0 1]+types);
f_der(9)=local_basis_3d(alpha,beta,gamma,orientations,types,[0 1 1]+types);
end