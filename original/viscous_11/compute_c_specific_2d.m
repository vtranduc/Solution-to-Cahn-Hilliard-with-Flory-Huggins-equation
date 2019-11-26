function sol=compute_c_specific_2d(...
    alpha,beta,eX,eY)

sol=zeros(4,4,3);

%basis(vals,orientation,type,order)
orientation_list=...
    [0 0
    1 0
    0 1
    1 1];

x_alpha=map_to_global_1_compt_2d(alpha,beta,eX,[1 0]);
x_beta=map_to_global_1_compt_2d(alpha,beta,eX,[0 1]);
x_alpha_beta=map_to_global_1_compt_2d(alpha,beta,eX,[1 1]);

y_alpha=map_to_global_1_compt_2d(alpha,beta,eY,[1 0]);
y_beta=map_to_global_1_compt_2d(alpha,beta,eY,[0 1]);
y_alpha_beta=map_to_global_1_compt_2d(alpha,beta,eY,[1 1]);

%-----------------------------------------------
J=[x_alpha,y_alpha;x_beta,y_beta];
% invJ_tiny=inv(invJ_tiny);

%------------------------------------------

f_local=zeros(2,1);

for orientation=1:1:4
    
    phi00=basis_2d(alpha,beta,orientation_list(orientation,:),[0 0],[0 0]);
    phi10=basis_2d(alpha,beta,orientation_list(orientation,:),[1 0],[0 0]);
    phi01=basis_2d(alpha,beta,orientation_list(orientation,:),[0 1],[0 0]);
    phi11=basis_2d(alpha,beta,orientation_list(orientation,:),[1 1],[0 0]);
    sol(orientation,1,1)=phi00;
    sol(orientation,2,1)=phi10*x_alpha+phi01*x_beta+phi11*x_alpha_beta;
    sol(orientation,3,1)=phi10*y_alpha+phi01*y_beta+phi11*y_alpha_beta;
    sol(orientation,4,1)=phi11*(x_beta*y_alpha+y_beta*x_alpha);
    
    phi10_alpha=basis_2d(alpha,beta,orientation_list(orientation,:),[1 0],[1 0]);
    phi01_alpha=basis_2d(alpha,beta,orientation_list(orientation,:),[0 1],[1 0]);
    phi11_alpha=basis_2d(alpha,beta,orientation_list(orientation,:),[1 1],[1 0]);
    
    phi10_beta=basis_2d(alpha,beta,orientation_list(orientation,:),[1 0],[0 1]);
    phi01_beta=basis_2d(alpha,beta,orientation_list(orientation,:),[0 1],[0 1]);
    phi11_beta=basis_2d(alpha,beta,orientation_list(orientation,:),[1 1],[0 1]);
    
    f_local(1)=phi10_alpha*x_alpha+phi01_alpha*x_beta+phi01*x_alpha_beta+phi11_alpha*x_alpha_beta;
    f_local(2)=phi10_beta*x_alpha+phi01_beta*x_beta+phi10*x_alpha_beta+phi11_beta*x_alpha_beta;
    take=J\f_local;
    sol(orientation,2,2)=take(1);
    sol(orientation,2,3)=take(2);
    
    f_local(1)=phi10_alpha*y_alpha+phi01_alpha*y_beta+phi01*y_alpha_beta+phi11_alpha*y_alpha_beta;
    f_local(2)=phi10_beta*y_alpha+phi01_beta*y_beta+phi10*y_alpha_beta+phi11_beta*y_alpha_beta;
    take=J\f_local;
    sol(orientation,3,2)=take(1);
    sol(orientation,3,3)=take(2);
    
    f_local(1)=phi11_alpha*(x_beta*y_alpha+x_alpha*y_beta)+phi11*(x_alpha_beta*y_alpha+x_alpha*y_alpha_beta);
    f_local(2)=phi11_beta*(x_beta*y_alpha+x_alpha*y_beta)+phi11*(x_beta*y_alpha_beta+x_alpha_beta*y_beta);
    take=J\f_local;
    sol(orientation,4,2)=take(1);
    sol(orientation,4,3)=take(2);
    
    
end

end

function sol=basis_2d(alpha,beta,orientations,types,orders)
sol=basis(alpha,orientations(1),types(1),orders(1))...
    *basis(beta,orientations(2),types(2),orders(2));
end