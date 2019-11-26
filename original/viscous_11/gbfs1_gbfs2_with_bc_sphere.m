function [bc_row_refiner,bc_col_refiner,bc_refiner]=...
    gbfs1_gbfs2_with_bc_sphere(neightTotal,neight,n_nSurf,normals)
len=n_nSurf(6)*3;
bc_refiner=zeros(1,len);
bc_row_refiner=zeros(1,len);
bc_col_refiner=zeros(1,len);
row=neightTotal;
index=0;
for i=1:1:n_nSurf(6)
    row=row+1;
    col=neight+(i-1)*8+1;
    for j=1:1:3
        index=index+1;
        bc_row_refiner(index)=row;
        bc_col_refiner(index)=col+j;
        bc_refiner(index)=normals(i,j);
    end
end
end