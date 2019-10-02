function solution=basis(vals,orientation,type,order)

%THIS IS LOCAL
%ALL ELEMENTS IN VALS MUST BE BETWEEN 0 AND 1

if order==0
    if type==0
        if orientation==0
            solution=1-3*vals.^2+2*vals.^3;
        elseif orientation==1
            solution=3*vals.^2-2*vals.^3;
        else
            error('Orientation must be 0 or 1!')
        end
    elseif type==1
        if orientation==0
            solution=vals-2*vals.^2+vals.^3;
        elseif orientation==1
            solution=-vals.^2+vals.^3;
        else
            error('Orientation must be 0 or 1!')
        end
    else
        error('Type must be 0 or 1!')
    end
elseif order==1
    if type==0
        if orientation==0
            solution=-6*vals+6*vals.^2;
        elseif orientation==1
            solution=6*vals-6*vals.^2;
        else
            error('Orientation must be 0 or 1!')
        end
    elseif type==1
        if orientation==0
            [nrows,ncols]=size(vals);
            solution=ones(nrows,ncols)-4*vals+3*vals.^2;
        elseif orientation==1
            solution=-2*vals+3*vals.^2;
        else
            error('Orientation must be 0 or 1!')
        end
    else
        error('Type must be 0 or 1!')
    end
elseif order==2
    if type==0
        if orientation==0
            [nrows,ncols]=size(vals);
            solution=-6*ones(nrows,ncols)+12*vals;
        elseif orientation==1
            [nrows,ncols]=size(vals);
            solution=6*ones(nrows,ncols)-12*vals;
        else
            error('Orientation must be 0 or 1!')
        end
    elseif type==1
        if orientation==0
            [nrows,ncols]=size(vals);
            solution=-4*ones(nrows,ncols)+6*vals;
        elseif orientation==1
            [nrows,ncols]=size(vals);
            solution=-2*ones(nrows,ncols)+6*vals;
        else
            error('Orientation must be 0 or 1!')
        end
    else
        error('Type must be 0 or 1!')
    end
    %==================
elseif order==3
    if type==0
        if orientation==0
            [nrows,ncols]=size(vals);
            solution=12*ones(nrows,ncols);
        elseif orientation==1
            [nrows,ncols]=size(vals);
            solution=-12*ones(nrows,ncols);
        else
            error('Orientation must be 0 or 1!')
        end
    elseif type==1
        if orientation==0
            [nrows,ncols]=size(vals);
            solution=6*ones(nrows,ncols);
        elseif orientation==1
            [nrows,ncols]=size(vals);
            solution=6*ones(nrows,ncols);
        else
            error('Orientation must be 0 or 1!')
        end
    else
        error('Type must be 0 or 1!')
    end
    
else
    error('Only derivatives of order 0, 1, 2, and 3 are available!')
    %=================================
end

end