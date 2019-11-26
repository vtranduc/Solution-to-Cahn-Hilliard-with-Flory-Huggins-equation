function sol=mirror(val,dir)
if dir==1
    if val==1
        sol=2;
    elseif val==3
        sol=4;
    elseif val==5
        sol=6;
    elseif val==7
        sol=8;
    end
elseif dir==2
    if val==2
        sol=1;
    elseif val==4
        sol=3;
    elseif val==6
        sol=5;
    elseif val==8
        sol=7;
    end
elseif dir==3
    if val==1
        sol=3;
    elseif val==2
        sol=4;
    elseif val==5
        sol=7;
    elseif val==6
        sol=8;
    end
elseif dir==4
    if val==3
        sol=1;
    elseif val==4
        sol=2;
    elseif val==7
        sol=5;
    elseif val==8
        sol=6;
    end
elseif dir==5
    if val==1
        sol=5;
    elseif val==2
        sol=6;
    elseif val==3
        sol=7;
    elseif val==4
        sol=8;
    end
elseif dir==6
    if val==5
        sol=1;
    elseif val==6
        sol=2;
    elseif val==7
        sol=3;
    elseif val==8
        sol=4;
    end
end
end