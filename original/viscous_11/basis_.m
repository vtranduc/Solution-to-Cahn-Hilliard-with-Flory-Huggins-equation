function solution=basis_(vals,orientation,type,order,alpha)

%THIS IS LOCAL
%ALL ELEMENTS IN VALS MUST BE BETWEEN 0 AND 1

if order==0
    if type==0
        if orientation==0
            [nrows,ncols]=size(vals);
            solution=ones(nrows,ncols)-3*vals.^2+2*vals.^3;
        elseif orientation==1
            solution=3*vals.^2-2*vals.^3;
        else
            error('Orientation must be 0 or 1!')
        end
    elseif type==1
        if orientation==0
            solution=(vals-2*vals.^2+vals.^3).*alpha;
        elseif orientation==1
            solution=(-vals.^2+vals.^3).*alpha;
        else
            error('Orientation must be 0 or 1!')
        end
    else
        error('Type must be 0 or 1!')
    end
elseif order==1
    if type==0
        if orientation==0
            solution=(-6*vals+6*vals.^2)./alpha;
        elseif orientation==1
            solution=(6*vals-6*vals.^2)./alpha;
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
            solution=(-6*ones(nrows,ncols)+12*vals)./(alpha^2);
        elseif orientation==1
            [nrows,ncols]=size(vals);
            solution=(6*ones(nrows,ncols)-12*vals)./(alpha^2);
        else
            error('Orientation must be 0 or 1!')
        end
    elseif type==1
        if orientation==0
            [nrows,ncols]=size(vals);
            solution=(-4*ones(nrows,ncols)+6*vals)./alpha;
        elseif orientation==1
            [nrows,ncols]=size(vals);
            solution=(-2*ones(nrows,ncols)+6*vals)./alpha;
        else
            error('Orientation must be 0 or 1!')
        end
    else
        error('Type must be 0 or 1!')
    end
else
    error('Only derivatives of order 0, 1, and 2 are available!')
end

end