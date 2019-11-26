function invJ=get_invJ(alpha,beta,gamma,eX,eY,eZ)
J1=get_Jacobian_part_1_3d(alpha,beta,gamma,eX,eY,eZ);
J5=get_Jacobian_part_5_3d(J1);
J6=get_Jacobian_part_6_3d(J1);
J7=get_Jacobian_part_7_3d(alpha,beta,gamma,eX,eY,eZ);
J8=get_Jacobian_part_8_3d(J1);
J9=get_Jacobian_part_9_3d(J1);
J=[J1 zeros(3,3) zeros(3,3)
    zeros(3,3) J5 J6
    J7 J8 J9];
invJ=inv(J);
end

function J9=get_Jacobian_part_9_3d(J1)
J9=zeros(3,3);
shift=...
    [1 2
    1 3
    2 3];
for irow=1:1:3
    for icol=1:1:3
        compt=shift(icol,:);
        orders=shift(irow,:);
        J9(irow,icol)=J1(orders(1),compt(1))*J1(orders(2),compt(2))...
            +J1(orders(2),compt(1))*J1(orders(1),compt(2));
    end
end
end

function J8=get_Jacobian_part_8_3d(J1)
J8=zeros(3,3);
for i=1:1:3
    J8(1,i)=J1(1,i)*J1(2,i);
    J8(2,i)=J1(1,i)*J1(3,i);
    J8(3,i)=J1(2,i)*J1(3,i);
end
end

function J7=get_Jacobian_part_7_3d(alpha,beta,gamma,eX,eY,eZ)
J7=zeros(3,3);
order_list=...
    [1 1 0
    1 0 1
    0 1 1];
for i=1:1:3
    J7(i,1)=map_to_global_1_compt_3d(alpha,beta,gamma,eX,order_list(i,:));
end
for i=1:1:3
    J7(i,2)=map_to_global_1_compt_3d(alpha,beta,gamma,eY,order_list(i,:));
end
for i=1:1:3
    J7(i,3)=map_to_global_1_compt_3d(alpha,beta,gamma,eZ,order_list(i,:));
end
end

function J6=get_Jacobian_part_6_3d(J1)
J6=zeros(3,3);
for i=1:1:3
    J6(i,1)=2*J1(i,1)*J1(i,2);
    J6(i,2)=2*J1(i,1)*J1(i,3);
    J6(i,3)=2*J1(i,2)*J1(i,3);
end
end

function J5=get_Jacobian_part_5_3d(J1)
J5=J1.^2.0;
end

function J1=get_Jacobian_part_1_3d(alpha,beta,gamma,eX,eY,eZ)
J1=zeros(3,3);
order_list=...
    [1 0 0
    0 1 0
    0 0 1];
for i=1:1:3
    J1(i,1)=map_to_global_1_compt_3d(alpha,beta,gamma,eX,order_list(i,:));
end
for i=1:1:3
    J1(i,2)=map_to_global_1_compt_3d(alpha,beta,gamma,eY,order_list(i,:));
end
for i=1:1:3
    J1(i,3)=map_to_global_1_compt_3d(alpha,beta,gamma,eZ,order_list(i,:));
end
end