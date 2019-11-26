function sol=map_to_global_3_compts_3d(alpha,beta,gamma,eX,eY,eZ,orders)
sol=zeros(1,3);
sol(1)=map_to_global_1_compt_3d(alpha,beta,gamma,eX,orders);
sol(2)=map_to_global_1_compt_3d(alpha,beta,gamma,eY,orders);
sol(3)=map_to_global_1_compt_3d(alpha,beta,gamma,eZ,orders);
end