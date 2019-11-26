function sol=isoBasis(val,orientation,order)
%val must be between 0 and 1
if order==0
    if orientation==0
        sol=1-val;
    elseif orientation==1
        sol=val;
    end
elseif order==1
    if orientation==0
        sol=-1;
    elseif orientation==1
        sol=1;
    end
end
end