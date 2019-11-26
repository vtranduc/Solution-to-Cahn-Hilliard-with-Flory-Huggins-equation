function [irows,icols,isj,col_placement,nnz_]=compact_rec_3D(n,nx,ny,nz)

nxny=nx*ny;

nnz_xy=48+36*((nx-2)+(ny-2))+27*(nx-2)*(ny-2);
nnz_xy_extreme=32+24*((nx-2)+(ny-2))+(nx-2)*(ny-2)*18;

nnz_=(nnz_xy)*(nz-2)+nnz_xy_extreme*2;

nnz_=nnz_*64;

irows=zeros(1,nnz_);
icols=zeros(1,nnz_);

isj=zeros(1,nnz_);


col_placement=1;
while col_placement<8*n
    col_placement=col_placement*10;
end

l=0;

for i=1:1:n
    interactions=interactions_rec_3D(i,ny,nz,nxny);
    lll=8*i-8;
    for ii=1:1:8
        llll=lll+ii;
        %============
        for j=1:1:length(interactions)
            ll=8*interactions(j)-8;
            for k=1:1:8
                l=l+1;
                irows(l)=llll;
                icols(l)=ll+k;
                isj(l)=irows(l)+icols(l)/col_placement;
            end
        end
        %====================
    end
end

end