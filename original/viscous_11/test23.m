function sj=test23(sj_,nx,ny,nnz_,irow,icol,correction)
sj=sj_;
counter=0;
for i=17:1:64
    if irow(i)==icol(i)
        sj(i)=1.0;
    else
        sj(i)=0.0;
    end
    
    counter=counter+1;index=i;
    row=irow(index);
    col=icol(index);
    if correction(row,col)~=0.0 && correction(row,col)~=1.0
        row
        col
        correction(row,col)
        error('wait')
    end
    
end
for i=89:96:89+96*(ny-3)
    for j=0:1:23
        
        if sj(i+j)==0.0
            error('Overlap')
        end
        
        if irow(i+j)==icol(i+j)
            sj(i+j)=1.0;
        else
            sj(i+j)=0.0;
        end
        
        counter=counter+1;index=i+j;
        row=irow(index);
        col=icol(index);
        if correction(row,col)~=0.0 && correction(row,col)~=1.0
            row
            col
            correction(row,col)
            error('wait')
        end
        
        
    end
end
for i=137:96:137+96*(ny-3)
    for j=0:1:23
        
        if sj(i+j)==0.0
            error('Overlap')
        end
        
        if irow(i+j)==icol(i+j)
            sj(i+j)=1.0;
        else
            sj(i+j)=0.0;
        end
        
        
        
        counter=counter+1;index=i+j;
        row=irow(index);
        col=icol(index);
        if correction(row,col)~=0.0 && correction(row,col)~=1.0
            row
            col
            correction(row,col)
            error('wait')
        end
        
    end
end
for i=81+96*(ny-2):1:96*(ny-2)+128
    
    if sj(i)==0.0
        error('Overlap')
    end
    
    if irow(i)==icol(i)
        sj(i)=1.0;
    else
        sj(i)=0.0;
    end
    
    
    
    counter=counter+1;index=i;
    row=irow(index);
    col=icol(index);
    if correction(row,col)~=0.0 && correction(row,col)~=1.0
        row
        col
        correction(row,col)
        error('wait')
    end
    
end
b=144*(ny-2)+192;
for i=96*(ny-2)+177:b:96*(ny-2)+177+b*(nx-3)
    for j=0:1:47
        
        if sj(i+j)==0.0
            error('Overlap')
        end
        
        if irow(i+j)==icol(i+j)
            sj(i+j)=1.0;
        else
            sj(i+j)=0.0;
        end
        
        counter=counter+1;index=i+j;
        row=irow(index);
        col=icol(index);
        if correction(row,col)~=0.0 && correction(row,col)~=1.0
            row
            col
            correction(row,col)
            error('wait')
        end
        
    end
end
%Top middle
counter12312=0;
for i=240*ny-207:b:240*ny-207+b*(nx-3)
    counter12312=counter12312+1;
    for j=0:1:47
        
        irow(i+j)
        icol(i+j)
        
        if sj(i+j)==0.0
            error('Overlap')
        end
        
        if irow(i+j)==icol(i+j)
            sj(i+j)=1.0;
        else
            sj(i+j)=0.0;
        end
        
        irow(i+j)
        icol(i+j)
        warning('yoiko')
        
%         counter=counter+1;index=i+j;
%         row=irow(index);
%         col=icol(index);
%         if correction(row,col)~=0.0 && correction(row,col)~=1.0
%             row
%             col
%             correction(row,col)
%             error('wait')
%         end
        
    end
end
counter12312
nx-2
warning('adfasdf')
%Bottom right
for i=96*(ny-2)+145+b*(nx-2):1:96*(ny-2)+192+b*(nx-2)
    
	if sj(i+j)==0.0
        error('Overlap')
	end
    
    if irow(i)==icol(i)
        sj(i)=1.0;
    else
        sj(i)=0.0;
    end
    
    counter=counter+1;index=i;
    row=irow(index);
    col=icol(index);
    if correction(row,col)~=0.0 && correction(row,col)~=1.0
        row
        col
        correction(row,col)
        error('wait')
    end
    
end
for i=96*(ny-2)+b*(nx-2)+217:96:96*(ny-2)+b*(nx-2)+217+96*(ny-3)
    for j=0:1:23
        
        if sj(i+j)==0.0
            error('Overlap')
        end
        
        if irow(i+j)==icol(i+j)
            sj(i+j)=1.0;
        else
            sj(i+j)=0.0;
        end
        
        counter=counter+1;index=i+j;
        row=irow(index);
        col=icol(index);
        if correction(row,col)~=0.0 && correction(row,col)~=1.0
            row
            col
            correction(row,col)
            error('wait')
        end
        
    end
end
%ERROR DETECTED HERE==========
for i=96*(ny-2)+b*(nx-2)+265:96:96*(ny-2)+b*(nx-2)+265+96*(ny-3)
    for j=0:1:23
        
        if sj(i+j)==0.0
            error('Overlap')
        end
        
        if irow(i+j)==icol(i+j)
            sj(i+j)=1.0;
        else
            sj(i+j)=0.0;
        end
        
        counter=counter+1;index=i+j;
        row=irow(index);
        col=icol(index);
        if correction(row,col)~=0.0 && correction(row,col)~=1.0
            
            row
            col
            correction(row,col)
            error('wait suspiciion')
        end
        
    end
end
for i=nnz_-47:1:nnz_
    
    if sj(i)==0.0
        error('Overlap')
    end
    
    if irow(i)==icol(i)
        sj(i)=1.0;
    else
        sj(i)=0.0;
    end
    
    counter=counter+1;index=i;
    row=irow(index);
    col=icol(index);
    if correction(row,col)~=0.0 && correction(row,col)~=1.0
        row
        col
        correction(row,col)
        error('wait')
    end
    
end
counter
z=4*(nx+ny)
warning('just stio')
end