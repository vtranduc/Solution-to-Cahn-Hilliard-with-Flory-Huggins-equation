function mx=generate_traditional_mx(rows,cols,vals)
mx=zeros(max(rows),max(cols));
for i=1:1:length(vals)
    mx(rows(i),cols(i))=vals(i);
end
end