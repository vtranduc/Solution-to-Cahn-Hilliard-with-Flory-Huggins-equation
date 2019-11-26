function [iZeros,iOnes]=get_bc_keyVals(nx,ny,nnz_,irow,icol)

b=144*(ny-2)+192;
iZeros=...
    [17 1 64
    89 96 89+96*(ny-3)
    137 96 137+96*(ny-3)
    81+96*(ny-2) 1 96*(ny-2)+128
    96*(ny-2)+177 b 96*(ny-2)+177+b*(nx-3)
    240*ny-207 b 240*ny-207+b*(nx-3)
    96*(ny-2)+145+b*(nx-2) 1 96*(ny-2)+192+b*(nx-2)
    96*(ny-2)+b*(nx-2)+217 96 96*(ny-2)+b*(nx-2)+217+96*(ny-3)
    96*(ny-2)+b*(nx-2)+265 96 96*(ny-2)+b*(nx-2)+265+96*(ny-3)
    nnz_-47 1 nnz_];

iOnes=bruteForceiOnesFinder(iZeros,irow,icol,nx,ny);

end

% function iOnes=mar25Test(iZeros,irow,icol,nx,ny)
% %This function is ugly
% iOnes=zeros(4*(nx+ny)-4,2);
% index=0;
% for i=iZeros(1,1):iZeros(1,2):iZeros(1,3)
%     if irow(i)==icol(i)
%         index=index+1;
%         iOnes(index,:)=[i irow(i)];
%     end
% end
% for i=iZeros(2,1):iZeros(2,2):iZeros(2,3)
%     for j=0:1:23
%         if irow(i+j)==icol(i+j)
%             index=index+1;
%             iOnes(index,:)=[i+j irow(i+j)];
%         end
%     end
% end
% for i=iZeros(3,1):iZeros(3,2):iZeros(3,3)
%     for j=0:1:23
%         if irow(i+j)==icol(i+j)
%             index=index+1;
%             iOnes(index,:)=[i+j irow(i+j)];
%         end
%     end
% end
% for i=iZeros(4,1):iZeros(4,2):iZeros(4,3)
%     if irow(i)==icol(i)
%         index=index+1;
%         iOnes(index,:)=[i irow(i)];
%     end
% end
% for i=iZeros(5,1):iZeros(5,2):iZeros(5,3)
%     for j=0:1:47
%         if irow(i+j)==icol(i+j)
%             index=index+1;
%             iOnes(index,:)=[i+j irow(i+j)];
%         end
%     end
% end
% for i=iZeros(6,1):iZeros(6,2):iZeros(6,3)
%     for j=0:1:47
%         if irow(i+j)==icol(i+j)
%             index=index+1;
%             iOnes(index,:)=[i+j irow(i+j)];
%         end
%     end
% end
% for i=iZeros(7,1):iZeros(7,2):iZeros(7,3)
%     if irow(i)==icol(i)
%         index=index+1;
%         iOnes(index,:)=[i irow(i)];
%     end
% end
% for i=iZeros(8,1):iZeros(8,2):iZeros(8,3)
%     for j=0:1:23
%         if irow(i+j)==icol(i+j)
%             index=index+1;
%             iOnes(index,:)=[i+j irow(i+j)];
%         end
%     end
% end
% for i=iZeros(9,1):iZeros(9,2):iZeros(9,3)
%     for j=0:1:23
%         if irow(i+j)==icol(i+j)
%             index=index+1;
%             iOnes(index,:)=[i+j irow(i+j)];
%         end
%     end
% end
% for i=iZeros(10,1):iZeros(10,2):iZeros(10,3)
%     if irow(i)==icol(i)
%         index=index+1;
%         iOnes(index,:)=[i irow(i)];
%     end
% end
% end

function iOnes=bruteForceiOnesFinder(iZeros,irow,icol,nx,ny)
%This algorithm is ugly
iOnes=zeros(1,4*(nx+ny)-4);
index=0;
for i=iZeros(1,1):iZeros(1,2):iZeros(1,3)
    if irow(i)==icol(i)
        index=index+1;
        iOnes(index)=i;
    end
end
for i=iZeros(2,1):iZeros(2,2):iZeros(2,3)
    for j=0:1:23
        if irow(i+j)==icol(i+j)
            index=index+1;
            iOnes(index)=i+j;
        end
    end
end
for i=iZeros(3,1):iZeros(3,2):iZeros(3,3)
    for j=0:1:23
        if irow(i+j)==icol(i+j)
            index=index+1;
            iOnes(index)=i+j;
        end
    end
end
for i=iZeros(4,1):iZeros(4,2):iZeros(4,3)
    if irow(i)==icol(i)
        index=index+1;
        iOnes(index)=i;
    end
end
for i=iZeros(5,1):iZeros(5,2):iZeros(5,3)
    for j=0:1:47
        if irow(i+j)==icol(i+j)
            index=index+1;
            iOnes(index)=i+j;
        end
    end
end
for i=iZeros(6,1):iZeros(6,2):iZeros(6,3)
    for j=0:1:47
        if irow(i+j)==icol(i+j)
            index=index+1;
            iOnes(index)=i+j;
        end
    end
end
for i=iZeros(7,1):iZeros(7,2):iZeros(7,3)
    if irow(i)==icol(i)
        index=index+1;
        iOnes(index)=i;
    end
end
for i=iZeros(8,1):iZeros(8,2):iZeros(8,3)
    for j=0:1:23
        if irow(i+j)==icol(i+j)
            index=index+1;
            iOnes(index)=i+j;
        end
    end
end
for i=iZeros(9,1):iZeros(9,2):iZeros(9,3)
    for j=0:1:23
        if irow(i+j)==icol(i+j)
            index=index+1;
            iOnes(index)=i+j;
        end
    end
end
for i=iZeros(10,1):iZeros(10,2):iZeros(10,3)
    if irow(i)==icol(i)
        index=index+1;
        iOnes(index)=i;
    end
end
end