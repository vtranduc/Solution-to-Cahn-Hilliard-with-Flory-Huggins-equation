function solution=conc(c_,weights,order_type)
%c_ is a vector containing 16 elements, which are nodal weights of a local
%element
solution=zeros(3,3);
for ixgp=1:1:3
    for iygp=1:1:3
        con=0;
        for ilocal=1:1:16
            con=con+...
                c_(ilocal)*weights(ixgp,iygp,ilocal,order_type);
        end
        solution(ixgp,iygp)=con;
    end
end
end