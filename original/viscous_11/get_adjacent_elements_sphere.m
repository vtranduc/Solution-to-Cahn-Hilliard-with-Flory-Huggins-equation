function [elements,eExtruding_x,eExtruding_y,eExtruding_z,...
    xOrientations,yOrientations,zOrientations]=...
    get_adjacent_elements_sphere(...
    node,n,nx,ny,nz,ne,nex,ney,n_eSurf,n_nSurf,nTotal,nx_)

if node<=n
    [xth,yth,zth]=get_xyzth_3d(node,nx,ny);
    elements=get_elements_3d(xth,yth,zth,nex,ney,nx,ny,nz);
    [eExtruding_x,eExtruding_y,eExtruding_z,xOrientations,yOrientations,zOrientations]=...
        get_adjacent_extruding_elements(...
        xth,yth,zth,nex,ney,nx,ny,nz,ne,n_eSurf,0);
elseif node<=nTotal
    elements=NaN;
    [xth,yth,zth]=get_root_of_extrusion_sph(node,n_nSurf,n,nx,ny,nz,nx_);
    [eExtruding_x,eExtruding_y,eExtruding_z,xOrientations,yOrientations,zOrientations]=...
        get_adjacent_extruding_elements(...
        xth,yth,zth,nex,ney,nx,ny,nz,ne,n_eSurf,1);
end

end

function [eExtruding_x,eExtruding_y,eExtruding_z,xOrientations,yOrientations,zOrientations]=...
    get_adjacent_extruding_elements(...
    xth,yth,zth,nex,ney,nx,ny,nz,ne,n_eSurf,node_is_on_surface)

if node_is_on_surface==0
    if xth==1
        xOrientations=[8 6 4 2];
        eExtruding_x=zeros(1,4);
        if yth~=1
            if zth~=1
                eExtruding_x(1)=get_esurf_minus_x(ney,yth-1,zth-1)+ne;
            end
            if zth~=nz
                eExtruding_x(3)=get_esurf_minus_x(ney,yth-1,zth)+ne;
            end
        end
        if yth~=ny
            if zth~=1
                eExtruding_x(2)=get_esurf_minus_x(ney,yth,zth-1)+ne;
            end
            if zth~=nz
                eExtruding_x(4)=get_esurf_minus_x(ney,yth,zth)+ne;
            end
        end
    elseif xth==nx
        xOrientations=[7 5 3 1];
        eExtruding_x=zeros(1,4);
        if yth~=1
            if zth~=1
                eExtruding_x(1)=get_esurf_plus_x(n_eSurf,ney,yth-1,zth-1)+ne;
            end
            if zth~=nz
                eExtruding_x(3)=get_esurf_plus_x(n_eSurf,ney,yth-1,zth)+ne;
            end
        end
        if yth~=ny
            if zth~=1
                eExtruding_x(2)=get_esurf_plus_x(n_eSurf,ney,yth,zth-1)+ne;
            end
            if zth~=nz
                eExtruding_x(4)=get_esurf_plus_x(n_eSurf,ney,yth,zth)+ne;
            end
        end
    else
        xOrientations=NaN;
        eExtruding_x=NaN;
    end

    if yth==1
        yOrientations=[8 7 4 3];
        eExtruding_y=zeros(1,4);
        if xth~=1
            if zth~=1
                eExtruding_y(1)=get_esurf_minus_y(n_eSurf,nex,xth-1,zth-1)+ne;
            end
            if zth~=nz
                eExtruding_y(3)=get_esurf_minus_y(n_eSurf,nex,xth-1,zth)+ne;
            end
        end
        if xth~=nx
            if zth~=1
                eExtruding_y(2)=get_esurf_minus_y(n_eSurf,nex,xth,zth-1)+ne;
            end
            if zth~=nz
                eExtruding_y(4)=get_esurf_minus_y(n_eSurf,nex,xth,zth)+ne;
            end
        end
    elseif yth==ny
        yOrientations=[6 5 2 1];
        eExtruding_y=zeros(1,4);
        if xth~=1
            if zth~=1
                eExtruding_y(1)=get_esurf_plus_y(n_eSurf,nex,xth-1,zth-1)+ne;
            end
            if zth~=nz
                eExtruding_y(3)=get_esurf_plus_y(n_eSurf,nex,xth-1,zth)+ne;
            end
        end
        if xth~=nx
            if zth~=1
                eExtruding_y(2)=get_esurf_plus_y(n_eSurf,nex,xth,zth-1)+ne;
            end
            if zth~=nz
                eExtruding_y(4)=get_esurf_plus_y(n_eSurf,nex,xth,zth)+ne;
            end
        end
    else
        yOrientations=NaN;
        eExtruding_y=NaN;
    end

    if zth==1
        zOrientations=[8 7 6 5];
        eExtruding_z=zeros(1,4);
        if xth~=1
            if yth~=1
                eExtruding_z(1)=get_esurf_minus_z(n_eSurf,nex,xth-1,yth-1)+ne;
            end
            if yth~=ny
                eExtruding_z(3)=get_esurf_minus_z(n_eSurf,nex,xth-1,yth)+ne;
            end
        end
        if xth~=nx
            if yth~=1
                eExtruding_z(2)=get_esurf_minus_z(n_eSurf,nex,xth,yth-1)+ne;
            end
            if yth~=ny
                eExtruding_z(4)=get_esurf_minus_z(n_eSurf,nex,xth,yth)+ne;
            end
        end
    elseif zth==nz
        zOrientations=[4 3 2 1];
        eExtruding_z=zeros(1,4);
        if xth~=1
            if yth~=1
                eExtruding_z(1)=get_esurf_plus_z(n_eSurf,nex,xth-1,yth-1)+ne;
            end
            if yth~=ny
                eExtruding_z(3)=get_esurf_plus_z(n_eSurf,nex,xth-1,yth)+ne;
            end
        end
        if xth~=nx
            if yth~=1
                eExtruding_z(2)=get_esurf_plus_z(n_eSurf,nex,xth,yth-1)+ne;
            end
            if yth~=ny
                eExtruding_z(4)=get_esurf_plus_z(n_eSurf,nex,xth,yth)+ne;
            end
        end
    else
        zOrientations=NaN;
        eExtruding_z=NaN;
    end
    %-------------------------------------------------
elseif node_is_on_surface==1
    if xth==1
        xOrientations=[7 5 3 1];
        eExtruding_x=zeros(1,4);
        if yth~=1
            if zth~=1
                eExtruding_x(1)=get_esurf_minus_x(ney,yth-1,zth-1)+ne;
            end
            if zth~=nz
                eExtruding_x(3)=get_esurf_minus_x(ney,yth-1,zth)+ne;
            end
        end
        if yth~=ny
            if zth~=1
                eExtruding_x(2)=get_esurf_minus_x(ney,yth,zth-1)+ne;
            end
            if zth~=nz
                eExtruding_x(4)=get_esurf_minus_x(ney,yth,zth)+ne;
            end
        end
    elseif xth==nx
        xOrientations=[8 6 4 2];
        eExtruding_x=zeros(1,4);
        if yth~=1
            if zth~=1
                eExtruding_x(1)=get_esurf_plus_x(n_eSurf,ney,yth-1,zth-1)+ne;
            end
            if zth~=nz
                eExtruding_x(3)=get_esurf_plus_x(n_eSurf,ney,yth-1,zth)+ne;
            end
        end
        if yth~=ny
            if zth~=1
                eExtruding_x(2)=get_esurf_plus_x(n_eSurf,ney,yth,zth-1)+ne;
            end
            if zth~=nz
                eExtruding_x(4)=get_esurf_plus_x(n_eSurf,ney,yth,zth)+ne;
            end
        end
    else
        xOrientations=NaN;
        eExtruding_x=NaN;
    end

    if yth==1
        yOrientations=[6 5 2 1];
        eExtruding_y=zeros(1,4);
        if xth~=1
            if zth~=1
                eExtruding_y(1)=get_esurf_minus_y(n_eSurf,nex,xth-1,zth-1)+ne;
            end
            if zth~=nz
                eExtruding_y(3)=get_esurf_minus_y(n_eSurf,nex,xth-1,zth)+ne;
            end
        end
        if xth~=nx
            if zth~=1
                eExtruding_y(2)=get_esurf_minus_y(n_eSurf,nex,xth,zth-1)+ne;
            end
            if zth~=nz
                eExtruding_y(4)=get_esurf_minus_y(n_eSurf,nex,xth,zth)+ne;
            end
        end
    elseif yth==ny
        yOrientations=[8 7 4 3];
        eExtruding_y=zeros(1,4);
        if xth~=1
            if zth~=1
                eExtruding_y(1)=get_esurf_plus_y(n_eSurf,nex,xth-1,zth-1)+ne;
            end
            if zth~=nz
                eExtruding_y(3)=get_esurf_plus_y(n_eSurf,nex,xth-1,zth)+ne;
            end
        end
        if xth~=nx
            if zth~=1
                eExtruding_y(2)=get_esurf_plus_y(n_eSurf,nex,xth,zth-1)+ne;
            end
            if zth~=nz
                eExtruding_y(4)=get_esurf_plus_y(n_eSurf,nex,xth,zth)+ne;
            end
        end
    else
        yOrientations=NaN;
        eExtruding_y=NaN;
    end

    if zth==1
        zOrientations=[4 3 2 1];
        eExtruding_z=zeros(1,4);
        if xth~=1
            if yth~=1
                eExtruding_z(1)=get_esurf_minus_z(n_eSurf,nex,xth-1,yth-1)+ne;
            end
            if yth~=ny
                eExtruding_z(3)=get_esurf_minus_z(n_eSurf,nex,xth-1,yth)+ne;
            end
        end
        if xth~=nx
            if yth~=1
                eExtruding_z(2)=get_esurf_minus_z(n_eSurf,nex,xth,yth-1)+ne;
            end
            if yth~=ny
                eExtruding_z(4)=get_esurf_minus_z(n_eSurf,nex,xth,yth)+ne;
            end
        end
    elseif zth==nz
        zOrientations=[8 7 6 5];
        eExtruding_z=zeros(1,4);
        if xth~=1
            if yth~=1
                eExtruding_z(1)=get_esurf_plus_z(n_eSurf,nex,xth-1,yth-1)+ne;
            end
            if yth~=ny
                eExtruding_z(3)=get_esurf_plus_z(n_eSurf,nex,xth-1,yth)+ne;
            end
        end
        if xth~=nx
            if yth~=1
                eExtruding_z(2)=get_esurf_plus_z(n_eSurf,nex,xth,yth-1)+ne;
            end
            if yth~=ny
                eExtruding_z(4)=get_esurf_plus_z(n_eSurf,nex,xth,yth)+ne;
            end
        end
    else
        zOrientations=NaN;
        eExtruding_z=NaN;
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

