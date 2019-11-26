function result=extract_abs_result_serendipity(c,nx,ny,nz)

result=zeros(ny,nx,nz);

l=-3;
for k=1:1:nz
    for n=1:1:ny
        for m=1:1:nx
            l=l+4;
            result(n,m,k)=c(l);
        end
    end
end

end