function test38

clear
clc

gps=[(-sqrt(3/5)+1)/2 0.5 (sqrt(3/5)+1)/2];


J1=get_Jacobian_part_1_3d(alpha,beta,gamma,eX,eY,eZ)


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

function sol=map_to_global_1_compt_3d(alpha,beta,gamma,eCompt,orders)
sol=eCompt(1)*isoBasis(alpha,0,orders(1))*isoBasis(beta,0,orders(2))*isoBasis(gamma,0,orders(3))...
    +eCompt(2)*isoBasis(alpha,1,orders(1))*isoBasis(beta,0,orders(2))*isoBasis(gamma,0,orders(3))...
    +eCompt(3)*isoBasis(alpha,0,orders(1))*isoBasis(beta,1,orders(2))*isoBasis(gamma,0,orders(3))...
    +eCompt(4)*isoBasis(alpha,1,orders(1))*isoBasis(beta,1,orders(2))*isoBasis(gamma,0,orders(3))...
    +eCompt(5)*isoBasis(alpha,0,orders(1))*isoBasis(beta,0,orders(2))*isoBasis(gamma,1,orders(3))...
    +eCompt(6)*isoBasis(alpha,1,orders(1))*isoBasis(beta,0,orders(2))*isoBasis(gamma,1,orders(3))...
    +eCompt(7)*isoBasis(alpha,0,orders(1))*isoBasis(beta,1,orders(2))*isoBasis(gamma,1,orders(3))...
    +eCompt(8)*isoBasis(alpha,1,orders(1))*isoBasis(beta,1,orders(2))*isoBasis(gamma,1,orders(3));
end

function sol=isoBasis(val,orientation,order)
%val must be between 0 and 1
if order==0
    if orientation==0
        sol=1-val;
    elseif orientation==1
        sol=val;
    end
elseif order==1
    if orientation==0
        sol=-1;
    elseif orientation==1
        sol=1;
    end
end
end