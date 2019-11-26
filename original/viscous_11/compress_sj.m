function [rowi,coli,val]=compress_sj(sj,nrows,ncols)

n=nnz(sj);
display('Setting up zeros')
rowi=zeros(1,n);
coli=zeros(1,n);
val=zeros(1,n);
display('Assigning')
k=0;
for i=1:1:nrows
    for j=1:1:ncols
        if sj(i,j)~=0
            k=k+1;
            rowi(k)=i;
            coli(k)=j;
            val(k)=sj(i,j);
        end
    end
end

end