function sol=get_extruding_element_sph(e,nex,ney,nez,nexney,n_eSurf)
sol=zeros(1,3);
[xth,yth,zth]=get_n1xyzth_3d(e,nex,nexney);
if xth==1
    if yth==1
        if zth==1
            sol(1)=get_esurf_minus_x(ney,yth,zth);
            sol(2)=get_esurf_minus_y(n_eSurf,nex,xth,zth);
            sol(3)=get_esurf_minus_z(n_eSurf,nex,xth,yth);
        elseif zth==nez
            sol(1)=get_esurf_minus_x(ney,yth,zth);
            sol(2)=get_esurf_minus_y(n_eSurf,nex,xth,zth);
            sol(3)=get_esurf_plus_z(n_eSurf,nex,xth,yth);
        else
            sol(1)=get_esurf_minus_x(ney,yth,zth);
            sol(2)=get_esurf_minus_y(n_eSurf,nex,xth,zth);
        end
    elseif yth==ney
        if zth==1
            sol(1)=get_esurf_minus_x(ney,yth,zth);
            sol(2)=get_esurf_plus_y(n_eSurf,nex,xth,zth);
            sol(3)=get_esurf_minus_z(n_eSurf,nex,xth,yth);
        elseif zth==nez
            sol(1)=get_esurf_minus_x(ney,yth,zth);
            sol(2)=get_esurf_plus_y(n_eSurf,nex,xth,zth);
            sol(3)=get_esurf_plus_z(n_eSurf,nex,xth,yth);
        else
            sol(1)=get_esurf_minus_x(ney,yth,zth);
            sol(2)=get_esurf_plus_y(n_eSurf,nex,xth,zth);
        end
    else
        if zth==1
            sol(1)=get_esurf_minus_x(ney,yth,zth);
            sol(3)=get_esurf_minus_z(n_eSurf,nex,xth,yth);
        elseif zth==nez
            sol(1)=get_esurf_minus_x(ney,yth,zth);
            sol(3)=get_esurf_plus_z(n_eSurf,nex,xth,yth);
        else
            sol(1)=get_esurf_minus_x(ney,yth,zth);
        end
    end
elseif xth==nex
    if yth==1
        if zth==1
            sol(1)=get_esurf_plus_x(n_eSurf,ney,yth,zth);
            sol(2)=get_esurf_minus_y(n_eSurf,nex,xth,zth);
            sol(3)=get_esurf_minus_z(n_eSurf,nex,xth,yth);
        elseif zth==nez
            sol(1)=get_esurf_plus_x(n_eSurf,ney,yth,zth);
            sol(2)=get_esurf_minus_y(n_eSurf,nex,xth,zth);
            sol(3)=get_esurf_plus_z(n_eSurf,nex,xth,yth);
        else
            sol(1)=get_esurf_plus_x(n_eSurf,ney,yth,zth);
            sol(2)=get_esurf_minus_y(n_eSurf,nex,xth,zth);
        end
    elseif yth==ney
        if zth==1
            sol(1)=get_esurf_plus_x(n_eSurf,ney,yth,zth);
            sol(2)=get_esurf_plus_y(n_eSurf,nex,xth,zth);
            sol(3)=get_esurf_minus_z(n_eSurf,nex,xth,yth);
        elseif zth==nez
            sol(1)=get_esurf_plus_x(n_eSurf,ney,yth,zth);
            sol(2)=get_esurf_plus_y(n_eSurf,nex,xth,zth);
            sol(3)=get_esurf_plus_z(n_eSurf,nex,xth,yth);
        else
            sol(1)=get_esurf_plus_x(n_eSurf,ney,yth,zth);
            sol(2)=get_esurf_plus_y(n_eSurf,nex,xth,zth);
        end
    else
        if zth==1
            sol(1)=get_esurf_plus_x(n_eSurf,ney,yth,zth);
            sol(3)=get_esurf_minus_z(n_eSurf,nex,xth,yth);
        elseif zth==nez
            sol(1)=get_esurf_plus_x(n_eSurf,ney,yth,zth);
            sol(3)=get_esurf_plus_z(n_eSurf,nex,xth,yth);
        else
            sol(1)=get_esurf_plus_x(n_eSurf,ney,yth,zth);
        end
    end
else
    if yth==1
        if zth==1
            sol(2)=get_esurf_minus_y(n_eSurf,nex,xth,zth);
            sol(3)=get_esurf_minus_z(n_eSurf,nex,xth,yth);
        elseif zth==nez
            sol(2)=get_esurf_minus_y(n_eSurf,nex,xth,zth);
            sol(3)=get_esurf_plus_z(n_eSurf,nex,xth,yth);
        else
            sol(2)=get_esurf_minus_y(n_eSurf,nex,xth,zth);
        end
    elseif yth==ney
        if zth==1
            sol(2)=get_esurf_plus_y(n_eSurf,nex,xth,zth);
            sol(3)=get_esurf_minus_z(n_eSurf,nex,xth,yth);
        elseif zth==nez
            sol(2)=get_esurf_plus_y(n_eSurf,nex,xth,zth);
            sol(3)=get_esurf_plus_z(n_eSurf,nex,xth,yth);
        else
            sol(2)=get_esurf_plus_y(n_eSurf,nex,xth,zth);
        end
    else
        if zth==1
            sol(3)=get_esurf_minus_z(n_eSurf,nex,xth,yth);
        elseif zth==nez
            sol(3)=get_esurf_plus_z(n_eSurf,nex,xth,yth);
        else
            error('The node is not on box surface')
        end
    end
end
end

function sol=get_esurf_minus_x(ney,yth,zth)
sol=(zth-1)*ney+yth;
end

function sol=get_esurf_plus_x(n_eSurf,ney,yth,zth)
sol=n_eSurf(1)+(zth-1)*ney+yth;
end

function sol=get_esurf_minus_y(n_eSurf,nex,xth,zth)
sol=n_eSurf(2)+(zth-1)*nex+xth;
end

function sol=get_esurf_plus_y(n_eSurf,nex,xth,zth)
sol=n_eSurf(3)+(zth-1)*nex+xth;
end

function sol=get_esurf_minus_z(n_eSurf,nex,xth,yth)
sol=n_eSurf(4)+(yth-1)*nex+xth;
end

function sol=get_esurf_plus_z(n_eSurf,nex,xth,yth)
sol=n_eSurf(5)+(yth-1)*nex+xth;
end