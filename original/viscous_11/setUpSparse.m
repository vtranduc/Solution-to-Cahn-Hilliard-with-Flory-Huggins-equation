function [irow,icol]=setUpSparse(nnz_,ny,n)
irow=zeros(1,nnz_);
icol=zeros(1,nnz_);
index=0;
for gbf1=1:1:4
    for gbf2=1:1:8
        index=index+1;
        irow(index)=gbf1;
        icol(index)=gbf2;
    end
    starter=ny*4;
    for i=1:1:8
        index=index+1;
        irow(index)=gbf1;
        icol(index)=starter+i;
    end
end
for e=2:1:ny-1
    starter=(e-2)*4;
    starter_=starter+4*ny;
    for i=1:1:4
        gbf1=(e-1)*4+i;
        for j=1:1:12
            index=index+1;
            irow(index)=gbf1;
            icol(index)=starter+j;
        end
        for j=1:1:12
            index=index+1;
            irow(index)=gbf1;
            icol(index)=starter_+j;
        end
    end
end
for gbf1=4*ny-3:1:ny*4
    for gbf2=(ny-2)*4+1:1:ny*4
        index=index+1;
        irow(index)=gbf1;
        icol(index)=gbf2;
    end
    for gbf2=8*ny-7:1:8*ny
        index=index+1;
        irow(index)=gbf1;
        icol(index)=gbf2;
    end
end
for node=ny+1:1:n-ny
    [~,edge_ydir]=edgy(node,ny,n);
    if edge_ydir==0
        for gbf1=(node-1)*4+1:1:4*node
            for gbf2=4*(node-ny)-7:1:4*(node-ny+1)
                index=index+1;
                irow(index)=gbf1;
                icol(index)=gbf2;
            end
            for gbf2=4*node-7:1:4*(node+1)
                index=index+1;
                irow(index)=gbf1;
                icol(index)=gbf2;
            end
            for gbf2=4*(node+ny)-7:1:4*(node+ny+1)
                index=index+1;
                irow(index)=gbf1;
                icol(index)=gbf2;
            end
        end
    elseif edge_ydir==-1
        for gbf1=(node-1)*4+1:1:4*node
            for gbf2=4*(node-ny)-3:1:4*(node-ny+1)
                index=index+1;
                irow(index)=gbf1;
                icol(index)=gbf2;
            end
            for gbf2=4*node-3:1:4*(node+1)
                index=index+1;
                irow(index)=gbf1;
                icol(index)=gbf2;
            end
            for gbf2=4*(node+ny)-3:1:4*(node+ny+1)
                index=index+1;
                irow(index)=gbf1;
                icol(index)=gbf2;
            end
        end
    elseif edge_ydir==1
        for gbf1=(node-1)*4+1:1:4*node
            for gbf2=4*(node-ny)-7:1:4*(node-ny)
                index=index+1;
                irow(index)=gbf1;
                icol(index)=gbf2;
            end
            for gbf2=4*node-7:1:4*node
                index=index+1;
                irow(index)=gbf1;
                icol(index)=gbf2;
            end
            for gbf2=4*(node+ny)-7:1:4*(node+ny)
                index=index+1;
                irow(index)=gbf1;
                icol(index)=gbf2;
            end
        end
    end 
end
node=node+1;
for gbf1=(node-1)*4+1:1:4*node
    for gbf2=4*(node-ny)-3:1:4*(node-ny+1)
        index=index+1;
        irow(index)=gbf1;
        icol(index)=gbf2;
    end
    for gbf2=4*node-3:1:4*(node+1)
        index=index+1;
        irow(index)=gbf1;
        icol(index)=gbf2;
    end
end
for node=n-ny+2:1:n-1
    for gbf1=(node-1)*4+1:1:4*node
        for gbf2=4*(node-ny)-7:1:4*(node-ny+1)
            index=index+1;
            irow(index)=gbf1;
            icol(index)=gbf2;
        end
        for gbf2=4*node-7:1:4*(node+1)
            index=index+1;
            irow(index)=gbf1;
            icol(index)=gbf2;
        end
    end
end
for gbf1=4*n-3:1:4*n
    for gbf2=4*(n-ny)-7:1:4*(n-ny)
        index=index+1;
        irow(index)=gbf1;
        icol(index)=gbf2;
    end
    for gbf2=4*n-7:1:4*n
        index=index+1;
        irow(index)=gbf1;
        icol(index)=gbf2;
    end
end
end