function simplified=simplify_2_regions(M,border)

[m,n,o]=size(M);
simplified=zeros(m,n,o);

for i=1:1:m
    for j=1:1:n
        for k=1:1:o
            val=M(i,j,k);
            if val>border
                simplified(i,j,k)=1;
            elseif val<border
                simplified(i,j,k)=-1;
            else
                simplified(i,j,k)=0;
            end
        end
    end
end

end