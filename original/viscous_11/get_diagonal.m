function sol=get_diagonal(diag,M)

[m,n]=size(M);

sol=zeros(1,m);

for i=1:1:m
    col=mod(diag+i-1,n);
    if col==0
        col=n;
    end
    sol(1,i)=M(i,col);
end

end