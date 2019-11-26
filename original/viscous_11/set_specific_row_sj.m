function sj=set_specific_row_sj(sj,value,row,nx,ny)

element=(row-mod(row-1,4)-1)/4+1;

if element<=ny %Left
    if element == 1 %left bottom
        middle=1;
        right=middle+4*ny;
        for i=0:1:7
            sj(row,middle+i)=value;
            sj(row,right+i)=value;
        end
    elseif element == ny %left top
        middle=4*(element-2)+1;
        right=middle+4*ny;
        for i=0:1:7
            sj(row,middle+i)=value;
            sj(row,right+i)=value;
        end
    else %left middle
        middle=4*(element-2)+1;
        right=middle+4*ny;
        for i=0:1:11
            sj(row,middle+i)=value;
            sj(row,right+i)=value;
        end
    end
elseif mod(element,ny)==0 %Top
    if element==nx*ny %Top right
        left=4*(element-ny-2)+1;
        middle=left+4*ny;
        for i=0:1:7
            sj(row,left+i)=value;
            sj(row,middle+i)=value;
        end
    else %Top middle
        left=4*(element-ny-2)+1;
        middle=left+4*ny;
        right=middle+4*ny;
        for i=0:1:7
            sj(row,left+i)=value;
            sj(row,middle+i)=value;
            sj(row,right+i)=value;
        end
    end
elseif mod(element-1,ny)==0 %Bottom
    if element==ny*(nx-1)+1 %Bottom right
        left=4*(element-ny-1)+1;
        middle=left+4*ny;
        for i=0:1:7
            sj(row,left+i)=value;
            sj(row,middle+i)=value;
        end
    else % Bottom middle
        left=4*(element-ny-1)+1;
        middle=left+4*ny;
        right=middle+4*ny;
        for i=0:1:7
            sj(row,left+i)=value;
            sj(row,middle+i)=value;
            sj(row,right+i)=value;
        end
    end
elseif element > ny*(nx-1)+1 %Right middle
    left=4*(element-ny-2)+1;
    middle=left+4*ny;
    for i=0:1:11
        sj(row,left+i)=value;
        sj(row,middle+i)=value;
    end
else % Middle
    left=4*(element-ny-2)+1;
    middle=left+4*ny;
    right=middle+4*ny;
    for i=0:1:11
        sj(row,left+i)=value;
        sj(row,middle+i)=value;
        sj(row,right+i)=value;
    end
end

end