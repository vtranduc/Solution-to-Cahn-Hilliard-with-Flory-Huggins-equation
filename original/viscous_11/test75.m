function [bc_row_refiner,bc_col_refiner,bc_refiner]=...
    test75(normal_nodal,bc_gbfs_r,n_bc_nodes,nfour)

bc_row_refiner=zeros(1,n_bc_nodes*3);
bc_col_refiner=bc_row_refiner;
bc_refiner=bc_row_refiner;

index=0;

for row=1:1:n_bc_nodes
    for compt=1:1:3
        index=index+1;
        bc_row_refiner(index)=row+nfour;
        bc_col_refiner(index)=bc_gbfs_r(row)+compt-1;
        bc_refiner(index)=normal_nodal(row,compt);
    end
end


end