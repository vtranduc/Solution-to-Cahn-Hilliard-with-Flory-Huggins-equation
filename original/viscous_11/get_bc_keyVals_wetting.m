function [iZeros,iOnes,ig]=get_bc_keyVals_wetting(...
    nx,ny,nnz_,irow,icol,hl,hb,hr,ht,gl,gb,gr,gt)

b=144*(ny-2)+192;

% [17 1 32
%     89 96 89+96*(ny-3)
%     81+96*(ny-2) 1 83+96*(ny-2)]
% 
% 17 96  
% 33 b 33*b(nx-1)


iZeros=[17 1 32 %Corner
    89 96 89+96*(ny-3) %x-dir on left
    81+96*(ny-2) 1 96*(ny-2)+96
    33 1 48
    96*(ny-2)+177 b 96*(ny-2)+177+b*(nx-3) %xy-term on left
    96*(ny-2)+161+b*(nx-2) 1 96*(ny-2)+176+b*(nx-2)
    96*(ny-2)+145+b*(nx-2) 1 96*(ny-2)+160+b*(nx-2)
    96*(ny-2)+b*(nx-2)+217 96 96*(ny-2)+b*(nx-2)+217+96*(ny-3)
    nnz_-47 1 nnz_-32
    97+96*(ny-2) 1 96*(ny-2)+112
    240*ny-207 b 240*ny-207+b*(nx-3)
    nnz_-31 1 nnz_-16
    49 1 64 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    137 96 137+96*(ny-3) %x-dir on left
    113+96*(ny-2) 1 96*(ny-2)+128
    96*(ny-2)+201 b 96*(ny-2)+201+b*(nx-3)
    96*(ny-2)+177+b*(nx-2) 1 96*(ny-2)+192+b*(nx-2)
    96*(ny-2)+b*(nx-2)+265 96 96*(ny-2)+b*(nx-2)+265+96*(ny-3)
    nnz_-15 1 nnz_
    240*ny-183 b 240*ny-183+b*(nx-3)];


% iZeros=...
%     [17 1 64 %Corner
%     89 96 89+96*(ny-3) %x-dir on left
%     137 96 137+96*(ny-3) %xy-term on left
%     81+96*(ny-2) 1 96*(ny-2)+128 %corner, top left
%     96*(ny-2)+177 b 96*(ny-2)+177+b*(nx-3) %bottom, xdir
%     240*ny-207 b 240*ny-207+b*(nx-3) %bottom, xy
%     96*(ny-2)+145+b*(nx-2) 1 96*(ny-2)+192+b*(nx-2) %corner, bottom right
%     96*(ny-2)+b*(nx-2)+217 96 96*(ny-2)+b*(nx-2)+217+96*(ny-3) % right, x-dir
%     96*(ny-2)+b*(nx-2)+265 96 96*(ny-2)+b*(nx-2)+265+96*(ny-3) %right, xy-dir
%     nnz_-47 1 nnz_]; %last corner

% ig=0;

% iOnes=0;
% return
[iOnes,ig]=bruteForceiOnesFinder(iZeros,irow,icol,nx,ny,...
    hl,hb,hr,ht,gl,gb,gr,gt,nnz_);

end

function [iOnes,ig]=bruteForceiOnesFinder(iZeros,irow,icol,nx,ny,...
    hl,hb,hr,ht,gl,gb,gr,gt,nnz_)
%This algorithm is ugly

%========================
%==========================
% VERY UGLY ALGORITHM HERE
ig=zeros(4*(nx+ny-1),2);

indexer=zeros(4*(nx+ny-1),3);
index=0;
n=nx*ny;
nfour=4*n;
il=0;
ib=0;
ir=0;
it=0;
for gbf1=1:1:nfour
    [node,class]=classify_gbf(gbf1);
    [inx,iny]=inxiny(node,ny);
    if inx==1
        if class==2
            index=index+1;
            il=il+1;
            indexer(index,:)=[gbf1,gbf1-1,-gl(il)];
        end
    end
    if inx==nx
        if class==2
            index=index+1;
            ir=ir+1;
            indexer(index,:)=[gbf1,gbf1-1,-gr(ir)];
        end
    end
    if iny==1
        if class==3
            index=index+1;
            ib=ib+1;
            indexer(index,:)=[gbf1,gbf1-2,-gb(ib)];
        end
    end
    if iny==ny
        if class==3
            index=index+1;
            it=it+1;
            indexer(index,:)=[gbf1,gbf1-2,-gt(it)];
        end
    end
end

%==================
il=0;
ib=1;
ir=0;
it=1;
for gbf1=1:1:nfour
    [node,class]=classify_gbf(gbf1);
    [inx,iny]=inxiny(node,ny);
    if inx==1
        if class==4
            index=index+1;
            il=il+1;
            indexer(index,:)=[gbf1,gbf1-1,-gl(il)];
        end
    end
    if inx==nx
        if class==4
            index=index+1;
            ir=ir+1;
            indexer(index,:)=[gbf1,gbf1-1,-gr(ir)];
        end
    end
    if iny==1 && inx>1 && inx<nx
        if class==4
            index=index+1;
            ib=ib+1;
            indexer(index,:)=[gbf1,gbf1-2,-gb(ib)];
        end
    end
    if iny==ny && inx>1 && inx<nx
        if class==4
            index=index+1;
            it=it+1;
            indexer(index,:)=[gbf1,gbf1-2,-gt(it)];
        end
    end
end
%==============

for i=1:1:4*(nx+ny-1)
    for j=1:1:nnz_
        if irow(j)==indexer(i,1) && icol(j)==indexer(i,2)
            ig(i,1)=j;
            ig(i,2)=indexer(i,3);
            break
        end
    end
end


%=========================
%===========================

iOnes=zeros(1,4*(nx+ny-1));

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
            break
        end
    end
end
for i=iZeros(3,1):iZeros(3,2):iZeros(3,3)
    if irow(i)==icol(i)
        index=index+1;
        iOnes(index)=i;
    end
end
for i=iZeros(4,1):iZeros(4,2):iZeros(4,3)
    if irow(i)==icol(i)
        index=index+1;
        iOnes(index)=i;
    end
end
for i=iZeros(5,1):iZeros(5,2):iZeros(5,3)
    for j=0:1:23
        if irow(i+j)==icol(i+j)
            index=index+1;
            iOnes(index)=i+j;
            break
        end
    end
end
for i=iZeros(6,1):iZeros(6,2):iZeros(6,3)
    if irow(i)==icol(i)
        index=index+1;
        iOnes(index)=i;
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
            break
        end
    end
end
for i=iZeros(9,1):iZeros(9,2):iZeros(9,3)
    if irow(i)==icol(i)
        index=index+1;
        iOnes(index)=i;
    end
end
for i=iZeros(10,1):iZeros(10,2):iZeros(10,3)
    if irow(i)==icol(i)
        index=index+1;
        iOnes(index)=i;
    end
end
for i=iZeros(11,1):iZeros(11,2):iZeros(11,3)
    for j=0:1:23
        if irow(i+j)==icol(i+j)
            index=index+1;
            iOnes(index)=i+j;
            break
        end
    end
end
for i=iZeros(12,1):iZeros(10,2):iZeros(12,3)
    if irow(i)==icol(i)
        index=index+1;
        iOnes(index)=i;
    end
end
%==================================================
%==================================================
%==================================================
%==================================================
%==================================================
%==================================================
for i=iZeros(13,1):iZeros(13,2):iZeros(13,3)
    if irow(i)==icol(i)
        index=index+1;
        iOnes(index)=i;
    end
end
for i=iZeros(14,1):iZeros(14,2):iZeros(14,3)
    for j=0:1:23
        if irow(i+j)==icol(i+j)
            index=index+1;
            iOnes(index)=i+j;
            break
        end
    end
end
for i=iZeros(15,1):iZeros(15,2):iZeros(15,3)
    if irow(i)==icol(i)
        index=index+1;
        iOnes(index)=i;
    end
end
for i=iZeros(16,1):iZeros(16,2):iZeros(16,3)
    for j=0:1:23
        if irow(i+j)==icol(i+j)
            index=index+1;
            iOnes(index)=i+j;
            break
        end
    end
end
for i=iZeros(17,1):iZeros(17,2):iZeros(17,3)
    if irow(i)==icol(i)
        index=index+1;
        iOnes(index)=i;
    end
end
for i=iZeros(18,1):iZeros(18,2):iZeros(18,3)
    for j=0:1:23
        if irow(i+j)==icol(i+j)
            index=index+1;
            iOnes(index)=i+j;
            break
        end
    end
end
for i=iZeros(19,1):iZeros(19,2):iZeros(19,3)
    if irow(i)==icol(i)
        index=index+1;
        iOnes(index)=i;
    end
end
for i=iZeros(20,1):iZeros(20,2):iZeros(20,3)
    for j=0:1:23
        if irow(i+j)==icol(i+j)
            index=index+1;
            iOnes(index)=i+j;
            break
        end
    end
end

end

function [inx,iny]=inxiny(node,ny)
iny=mod(node,ny);
if iny==0
    iny=ny;
end
inx=(node-iny)/ny+1;
end