function [diffT,chi]=T_gradient_3d(T,diff,nex,entropy,T_theta)
if length(T)==2
    gp=[(-sqrt(3/5)+1)/2 0.5 (sqrt(3/5)+1)/2];
    diffT=zeros(nex,3);
    chi=zeros(nex,3);
    nex_increment=(T(2)-T(1))/nex;
    gp_increment=nex_increment*gp;
    for i=1:1:nex
        adder=T(1)+(i-1)*nex_increment;
        for j=1:1:3
            diffT(i,j)=adder+gp_increment(j);
        end
    end
    for i=1:1:nex
        for j=1:1:3
            chi(i,j)=0.5-entropy*(1-T_theta/diffT(i,j));
        end
    end
    diffT=diffT*diff;
elseif length(T)==1
    diffT=diff*T;
    chi=0.5-entropy*(1-T_theta/T);
end
end