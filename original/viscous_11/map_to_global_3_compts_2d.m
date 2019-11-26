function sol=map_to_global_3_compts_2d(alpha,beta,eX,eY,orders)
sol=zeros(1,2);
sol(1)=map_to_global_1_compt_2d(alpha,beta,eX,orders);
sol(2)=map_to_global_1_compt_2d(alpha,beta,eY,orders);
end