function [freezing_zeros,freezing_ones]=freezing_setUp(irow,icol,nx,ny)

%FINDS ONES BY BRUTE FORCE HERE!

freezing_zeros=zeros(nx,(ny-2)*144+192);
freezing_ones=zeros(nx,ny*4);


nMiddle=(ny-2)*144+192;
nEdge=(ny-2)*96+128;

index=0;
index_zero=0;
for i=1:1:nEdge
    index=index+1;
    freezing_zeros(1,i)=index;
    if irow(index)==icol(index)
        index_zero=index_zero+1;
        freezing_ones(1,index_zero)=index;
    end
end

for ix=2:1:nx-1
    index_zero=0;
    for i=1:1:nMiddle
        index=index+1;
        freezing_zeros(ix,i)=index;
        if irow(index)==icol(index)
            index_zero=index_zero+1;
            freezing_ones(ix,index_zero)=index;
        end
    end
end

index_zero=0;
for i=1:1:nEdge
    index=index+1;
    freezing_zeros(nx,i)=index;
    if irow(index)==icol(index)
        index_zero=index_zero+1;
        freezing_ones(nx,index_zero)=index;
    end
end

% freezing_zeros';
% 
% freezing_ones';
% 
% 
% %------------------------------
% test_col=4;
% test1=zeros(1,ny*4);
% test2=zeros(1,ny*4);
% test0=zeros(1,ny*4);
% for i=1:1:ny*4
%     index=freezing_ones(test_col,i);
%     test1(i)=irow(index);
%     test2(i)=icol(index);
%     test0(i)=index;
% end
% 
% [test0' test1' test2']
%-------------------------------



end