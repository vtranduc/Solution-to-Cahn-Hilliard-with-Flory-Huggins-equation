function [sol,weights]=compute_c_specific(nodal_vals,orders_to_be_computed,alpha,beta,gamma,gps,w,eX,eY,eZ,orientation_list)

% orders_to_be_computed must be value(s) between 1 and 7.
% Note that weights of all orders will be computed, so it is wasteful to
% compute not all order.
% nodal_vals are nodal weights. It must be a vector containing 64 elements.
len=length(orders_to_be_computed);
sol=zeros(1,len);
weights=compute_weights_specific_3d(...
    alpha,beta,gamma,gps,w,eX,eY,eZ,orientation_list);
iorder=0;
for order=orders_to_be_computed
    iorder=iorder+1;
    inode=0;
    for orientation=1:1:8
        for type=1:1:8
            inode=inode+1;
            sol(iorder)=sol(iorder)+nodal_vals(inode)...
                *weights(orientation,type,order);
        end
    end
end
end