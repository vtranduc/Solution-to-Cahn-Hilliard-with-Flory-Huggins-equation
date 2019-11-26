function sj=sparseBC(sj_,nx,ny,nnz_)
sj=sj_;
for i=17:1:64
    sj(i)=0.0;
end
for i=89:96:89+96*(ny-3)
    for j=0:1:23
        sj(i+j)=0.0;
    end
end
for i=137:96:137+96*(ny-3)
    for j=0:1:23
        sj(i+j)=0.0;
    end
end
for i=81+96*(ny-2):1:96*(ny-2)+128
    sj(i)=0.0;
end
b=144*(ny-2)+192;
for i=96*(ny-2)+177:b:96*(ny-2)+177+b*(nx-3)
    for j=0:1:27
        sj(i+j)=0.0;
    end
end
for i=96*(ny-2)+81:b:96*(ny-2)+81+b*(ny-3)
    for j=0:1:23
        sj(i+j)=0.0;
    end
end
for i=96*(ny-2)+145+b*(nx-2):1:96*(ny-2)+192+b*(nx-2)
    sj(i)=0.0;
end
for i=96*(ny-2)+b*(nx-2)+217:96:96*(ny-2)+b*(nx-2)+217+96*(ny-3)
    for j=0:1:23
        sj(i+j)=0.0;
    end
end
for i=96*(ny-2)+b*(nx-2)+241:96:96*(ny-2)+b*(nx-2)+241+96*(ny-3)
    for j=0:1:23
        sj(i+j)=0.0;
    end
end
for i=nnz_-47:1:nnz_
    sj(i)=0.0;
end
end