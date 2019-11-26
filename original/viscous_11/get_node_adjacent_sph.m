function [nodes_adjacent,extruding]=get_node_adjacent_sph(node,n,nx,ny,nz,nxny,n_nSurf,nx_,nTotal)
%If loop can be shortened, so that less statements will be checked.
%This can be fixed later if time permitted


if node<=n
    
    nodes_adjacent=zeros(1,27);
    
    [xth,yth,zth]=get_xyzth_3d(node,nx,ny);
    if zth~=1
        node_minus_nxny=node-nx*ny;
        if xth~=1 && yth~=1
            nodes_adjacent(1)=node_minus_nxny-nx-1;
        end
        if yth~=1
            nodes_adjacent(2)=node_minus_nxny-nx;
        end
        if xth~=nx && yth~=1
            nodes_adjacent(3)=node_minus_nxny-nx+1;
        end
        if xth~=1
            nodes_adjacent(4)=node_minus_nxny-1;
        end

        nodes_adjacent(5)=node_minus_nxny;

        if xth~=nx
            nodes_adjacent(6)=node_minus_nxny+1;
        end
        if xth~=1 && yth~=ny
            nodes_adjacent(7)=node_minus_nxny+nx-1;
        end
        if yth~=ny
            nodes_adjacent(8)=node_minus_nxny+nx;
        end
        if xth~=nx && yth~=ny
            nodes_adjacent(9)=node_minus_nxny+nx+1;
        end
    end
    if xth~=1 && yth~=1
        nodes_adjacent(10)=node-nx-1;
    end
    if yth~=1
        nodes_adjacent(11)=node-nx;
    end
    if xth~=nx && yth~=1
        nodes_adjacent(12)=node-nx+1;
    end
    if xth~=1
        nodes_adjacent(13)=node-1;
    end

    nodes_adjacent(14)=node;

    if xth~=nx
        nodes_adjacent(15)=node+1;
    end
    if xth~=1 && yth~=ny
        nodes_adjacent(16)=node+nx-1;
    end
    if yth~=ny
        nodes_adjacent(17)=node+nx;
    end
    if xth~=nx && yth~=ny
        nodes_adjacent(18)=node+nx+1;
    end
    if zth~=nz
        node_plus_nxny=node+nx*ny;
        if xth~=1 && yth~=1
            nodes_adjacent(19)=node_plus_nxny-nx-1;
        end
        if yth~=1
            nodes_adjacent(20)=node_plus_nxny-nx;
        end
        if xth~=nx && yth~=1
            nodes_adjacent(21)=node_plus_nxny-nx+1;
        end
        if xth~=1
            nodes_adjacent(22)=node_plus_nxny-1;
        end

        nodes_adjacent(23)=node_plus_nxny;

        if xth~=nx
            nodes_adjacent(24)=node_plus_nxny+1;
        end
        if xth~=1 && yth~=ny
            nodes_adjacent(25)=node_plus_nxny+nx-1;
        end
        if yth~=ny
            nodes_adjacent(26)=node_plus_nxny+nx;
        end
        if xth~=nx && yth~=ny
            nodes_adjacent(27)=node_plus_nxny+nx+1;
        end
    end
    % --- Compute extrusion ---------------------------------------------
    if xth==1
        if yth==1
            if zth==1
                csurf=zeros(1,7);
                csurf(1)=1;
                csurf(2)=2;
                csurf(3)=nx+1;
                csurf(4)=nx+2;
                csurf(5)=nxny+1;
                csurf(6)=nxny+2;
                csurf(7)=nxny+nx+1;
                extruding=get_extrusion_sph(csurf,n,nx,ny,nz,n_nSurf,nx_);
            elseif zth==nz
                csurf=zeros(1,7);
                z1=(nz-2)*nxny+1;
                csurf(1)=z1;
                csurf(2)=z1+1;
                csurf(3)=z1+nx;
                z2=(nz-1)*nxny+1;
                csurf(4)=z2;
                csurf(5)=z2+1;
                csurf(6)=z2+nx;
                csurf(7)=csurf(6)+1;
                extruding=get_extrusion_sph(csurf,n,nx,ny,nz,n_nSurf,nx_);
            else
                csurf=zeros(1,9);
                z2=get_node_3d(xth,yth,zth,nx,ny);
                z1=z2-nxny;
                z3=z2+nxny;
                csurf(1)=z1;
                csurf(2)=z1+1;
                csurf(3)=z1+nx;
                csurf(4)=z2;
                csurf(5)=z2+1;
                csurf(6)=z2+nx;
                csurf(7)=z3;
                csurf(8)=z3+1;
                csurf(9)=z3+nx;
                extruding=get_extrusion_sph(csurf,n,nx,ny,nz,n_nSurf,nx_);
            end
        elseif yth==ny
            if zth==1
                csurf=zeros(1,7);
                z1=nxny-2*nx+1;
                z2=2*nxny-nx+1;
                csurf(1)=z1;
                csurf(2)=z1+1;
                csurf(3)=z1+nx;
                csurf(4)=csurf(3)+1;
                csurf(5)=z2-nx;
                csurf(6)=z2;
                csurf(7)=z2+1;
                extruding=get_extrusion_sph(csurf,n,nx,ny,nz,n_nSurf,nx_);
            elseif zth==nz
                csurf=zeros(1,7);
                z2=get_node_3d(xth,yth,zth,nx,ny);
                z1=z2-nxny;
                csurf(1)=z1-nx;
                csurf(2)=z1;
                csurf(3)=z1+1;
                csurf(4)=z2-nx;
                csurf(5)=csurf(4)+1;
                csurf(6)=z2;
                csurf(7)=z2+1;
                extruding=get_extrusion_sph(csurf,n,nx,ny,nz,n_nSurf,nx_);
            else
                csurf=zeros(1,9);
                z2=get_node_3d(xth,yth,zth,nx,ny);
                z1=z2-nxny;
                z3=z2+nxny;
                csurf(1)=z1-nx;
                csurf(2)=z1;
                csurf(3)=z1+1;
                csurf(4)=z2-nx;
                csurf(5)=z2;
                csurf(6)=z2+1;
                csurf(7)=z3-nx;
                csurf(8)=z3;
                csurf(9)=z3+1;
                extruding=get_extrusion_sph(csurf,n,nx,ny,nz,n_nSurf,nx_);
            end
        else
            if zth==1
                csurf=zeros(1,9);
                csurf(1)=node-nx;
                csurf(2)=csurf(1)+1;
                csurf(3)=node;
                csurf(4)=node+1;
                csurf(5)=node+nx;
                csurf(6)=csurf(5)+1;
                csurf(8)=node+nxny;
                csurf(7)=csurf(8)-nx;
                csurf(9)=csurf(8)+nx;
                extruding=get_extrusion_sph(csurf,n,nx,ny,nz,n_nSurf,nx_);
            elseif zth==nz
                csurf=zeros(1,9);
                z1=node-nxny;
                csurf(1)=z1-nx;
                csurf(2)=z1;
                csurf(3)=z1+nx;
                csurf(4)=node-nx;
                csurf(5)=node;
                csurf(6)=node+nx;
                z1=node+1;
                csurf(7)=z1-nx;
                csurf(8)=z1;
                csurf(9)=z1+nx;
                extruding=get_extrusion_sph(csurf,n,nx,ny,nz,n_nSurf,nx_);
            else
                csurf=zeros(1,9);
                z2=node;
                z1=node-nxny;
                z3=node+nxny;
                csurf(1)=z1-nx;
                csurf(2)=z1;
                csurf(3)=z1+nx;
                csurf(4)=z2-nx;
                csurf(5)=z2;
                csurf(6)=z2+nx;
                csurf(7)=z3-nx;
                csurf(8)=z3;
                csurf(9)=z3+nx;
                extruding=get_extrusion_sph(csurf,n,nx,ny,nz,n_nSurf,nx_);
            end
        end
    elseif xth==nx
        if yth==1
            if zth==1
                csurf=zeros(1,7);
                csurf(1)=node+nx-1;
                csurf(2)=node-1;
                csurf(3)=node;
                csurf(4)=node+nx;
                csurf(6)=node+nxny;
                csurf(5)=csurf(6)-1;
                csurf(7)=csurf(6)+nx;
                extruding=get_extrusion_sph(csurf,n,nx,ny,nz,n_nSurf,nx_);
            elseif zth==nz
                csurf=zeros(1,7);
                csurf(2)=node-nxny;
                csurf(1)=csurf(2)-1;
                csurf(3)=csurf(2)+nx;
                csurf(4)=node-1;
                csurf(5)=node;
                csurf(7)=node+nx;
                csurf(6)=csurf(7)-1;
                extruding=get_extrusion_sph(csurf,n,nx,ny,nz,n_nSurf,nx_);
            else
                csurf=zeros(1,9);
                z1=node-nxny;
                z2=node;
                z3=node+nxny;
                csurf(1)=z1-1;
                csurf(2)=z1;
                csurf(3)=z1+nx;
                csurf(4)=z2-1;
                csurf(5)=z2;
                csurf(6)=z2+nx;
                csurf(7)=z3-1;
                csurf(8)=z3;
                csurf(9)=z3+nx;
                extruding=get_extrusion_sph(csurf,n,nx,ny,nz,n_nSurf,nx_);
            end
        elseif yth==ny
            if zth==1
                csurf=zeros(1,7);
                csurf(2)=node-nx;
                csurf(1)=csurf(2)-1;
                csurf(3)=node-1;
                csurf(4)=node;
                csurf(7)=node+nxny;
                csurf(5)=csurf(7)-nx;
                csurf(6)=csurf(7)-1;
                extruding=get_extrusion_sph(csurf,n,nx,ny,nz,n_nSurf,nx_);
            elseif zth==nz
                csurf=zeros(1,7);
                csurf(3)=node-nxny;
                csurf(1)=csurf(3)-nx;
                csurf(2)=csurf(3)-1;
                csurf(5)=node-nx;
                csurf(4)=csurf(5)-1;
                csurf(6)=node-1;
                csurf(7)=node;
                extruding=get_extrusion_sph(csurf,n,nx,ny,nz,n_nSurf,nx_);
            else
                csurf=zeros(1,9);
                z1=node-nxny;
                z2=node;
                z3=node+nxny;
                csurf(1)=z1-nx;
                csurf(2)=z1-1;
                csurf(3)=z1;
                csurf(4)=z2-nx;
                csurf(5)=z2-1;
                csurf(6)=z2;
                csurf(7)=z3-nx;
                csurf(8)=z3-1;
                csurf(9)=z3;
                extruding=get_extrusion_sph(csurf,n,nx,ny,nz,n_nSurf,nx_);
            end
        else
            if zth==1
                csurf=zeros(1,9);
                csurf(2)=node-nx;
                csurf(1)=csurf(2)-1;
                csurf(3)=node-1;
                csurf(4)=node;
                csurf(6)=node+nx;
                csurf(5)=csurf(6)-1;
                csurf(8)=node+nxny;
                csurf(7)=csurf(8)-nx;
                csurf(9)=csurf(8)+nx;
                extruding=get_extrusion_sph(csurf,n,nx,ny,nz,n_nSurf,nx_);
            elseif zth==nz
                csurf=zeros(1,9);
                csurf(2)=node-nxny;
                csurf(1)=csurf(2)-nx;
                csurf(3)=csurf(2)+nx;
                csurf(5)=node-nx;
                csurf(4)=csurf(5)-1;
                csurf(6)=node-1;
                csurf(7)=node;
                csurf(9)=node+nx;
                csurf(8)=csurf(9)-1;
                extruding=get_extrusion_sph(csurf,n,nx,ny,nz,n_nSurf,nx_);
            else
                z1=node-nxny;
                z2=node;
                z3=node+nxny;
                csurf(1)=z1-nx;
                csurf(2)=z1;
                csurf(3)=z1+nx;
                csurf(4)=z2-nx;
                csurf(5)=z2;
                csurf(6)=z2+nx;
                csurf(7)=z3-nx;
                csurf(8)=z3;
                csurf(9)=z3+nx;
                extruding=get_extrusion_sph(csurf,n,nx,ny,nz,n_nSurf,nx_);
            end
        end
    else
        if yth==1
            if zth==1
                csurf=zeros(1,9);
                csurf(1)=node-1;
                csurf(2)=node;
                csurf(3)=node+1;
                z1=node+nx;
                csurf(4)=z1-1;
                csurf(5)=z1;
                csurf(6)=z1+1;
                z2=node+nxny;
                csurf(7)=z2-1;
                csurf(8)=z2;
                csurf(9)=z2+1;
                extruding=get_extrusion_sph(csurf,n,nx,ny,nz,n_nSurf,nx_);
            elseif zth==nz
                csurf=zeros(1,9);
                z1=node-nxny;
                z2=node;
                z3=node+nx;
                csurf(1)=z1-1;
                csurf(2)=z1;
                csurf(3)=z1+1;
                csurf(4)=z2-1;
                csurf(5)=z2;
                csurf(6)=z2+1;
                csurf(7)=z3-1;
                csurf(8)=z3;
                csurf(9)=z3+1;
                extruding=get_extrusion_sph(csurf,n,nx,ny,nz,n_nSurf,nx_);
            else
                csurf=zeros(1,9);
                z1=node-nxny;
                z2=node;
                z3=node+nxny;
                csurf(1)=z1-1;
                csurf(2)=z1;
                csurf(3)=z1+1;
                csurf(4)=z2-1;
                csurf(5)=z2;
                csurf(6)=z2+1;
                csurf(7)=z3-1;
                csurf(8)=z3;
                csurf(9)=z3+1;
                extruding=get_extrusion_sph(csurf,n,nx,ny,nz,n_nSurf,nx_);
            end
        elseif yth==ny
            if zth==1
                csurf=zeros(1,9);
                z1=node-nx;
                z2=node;
                z3=node+nxny;
                csurf(1)=z1-1;
                csurf(2)=z1;
                csurf(3)=z1+1;
                csurf(4)=z2-1;
                csurf(5)=z2;
                csurf(6)=z2+1;
                csurf(7)=z3-1;
                csurf(8)=z3;
                csurf(9)=z3+1;
                extruding=get_extrusion_sph(csurf,n,nx,ny,nz,n_nSurf,nx_);
            elseif zth==nz
                csurf=zeros(1,9);
                z1=node-nxny;
                z2=node-nx;
                z3=node;
                csurf(1)=z1-1;
                csurf(2)=z1;
                csurf(3)=z1+1;
                csurf(4)=z2-1;
                csurf(5)=z2;
                csurf(6)=z2+1;
                csurf(7)=z3-1;
                csurf(8)=z3;
                csurf(9)=z3+1;
                extruding=get_extrusion_sph(csurf,n,nx,ny,nz,n_nSurf,nx_);
            else
                csurf=zeros(1,9);
                z1=node-nxny;
                z2=node;
                z3=node+nxny;
                csurf(1)=z1-1;
                csurf(2)=z1;
                csurf(3)=z1+1;
                csurf(4)=z2-1;
                csurf(5)=z2;
                csurf(6)=z2+1;
                csurf(7)=z3-1;
                csurf(8)=z3;
                csurf(9)=z3+1;
                extruding=get_extrusion_sph(csurf,n,nx,ny,nz,n_nSurf,nx_);
            end
        else
            if zth==1 || zth==nz
                csurf=zeros(1,9);
                z1=node-nx;
                z2=node;
                z3=node+nx;
                csurf(1)=z1-1;
                csurf(2)=z1;
                csurf(3)=z1+1;
                csurf(4)=z2-1;
                csurf(5)=z2;
                csurf(6)=z2+1;
                csurf(7)=z3-1;
                csurf(8)=z3;
                csurf(9)=z3+1;
                extruding=get_extrusion_sph(csurf,n,nx,ny,nz,n_nSurf,nx_);
            else
                extruding=NaN;
            end
        end
    end
    % --- Extruding nodes ----------------------------------------------
    % --- Extruding nodes ----------------------------------------------
    % --- Extruding nodes ----------------------------------------------
    % --- Extruding nodes ----------------------------------------------
    % --- Extruding nodes ----------------------------------------------
    % --- Extruding nodes ----------------------------------------------
    % --- Extruding nodes ----------------------------------------------
    % --- Extruding nodes ----------------------------------------------
    % --- Extruding nodes ----------------------------------------------
    % --- Extruding nodes ----------------------------------------------
    % --- Extruding nodes ----------------------------------------------
    % --- Extruding nodes ----------------------------------------------
    % --- Extruding nodes ----------------------------------------------
    % --- Extruding nodes ----------------------------------------------
    % --- Extruding nodes ----------------------------------------------
    % --- Extruding nodes ----------------------------------------------
    % --- Extruding nodes ----------------------------------------------
    % --- Extruding nodes ----------------------------------------------
    % --- Extruding nodes ----------------------------------------------
    % --- Extruding nodes ----------------------------------------------
    % --- Extruding nodes ----------------------------------------------
    % --- Extruding nodes ----------------------------------------------
    % --- Extruding nodes ----------------------------------------------
    % --- Extruding nodes ----------------------------------------------
    % --- Extruding nodes ----------------------------------------------
    % --- Extruding nodes ----------------------------------------------
    % --- Extruding nodes ----------------------------------------------
    % --- Extruding nodes ----------------------------------------------
    % --- Extruding nodes ----------------------------------------------
    % --- Extruding nodes ----------------------------------------------
    % --- Extruding nodes ----------------------------------------------
    % --- Extruding nodes ----------------------------------------------
    % --- Extruding nodes ----------------------------------------------
    % --- Extruding nodes ----------------------------------------------
    % --- Extruding nodes ----------------------------------------------
    % --- Extruding nodes ----------------------------------------------
    % --- Extruding nodes ----------------------------------------------
    % --- Extruding nodes ----------------------------------------------
    % --- Extruding nodes ----------------------------------------------
    % --- Extruding nodes ----------------------------------------------
    % --- Extruding nodes ----------------------------------------------
    % --- Extruding nodes ----------------------------------------------
    % --- Extruding nodes ----------------------------------------------
    % --- Extruding nodes ----------------------------------------------
    % --- Extruding nodes ----------------------------------------------
    % --- Extruding nodes ----------------------------------------------
    % --- Extruding nodes ----------------------------------------------
    % --- Extruding nodes ----------------------------------------------
    % --- Extruding nodes ----------------------------------------------
    % --- Extruding nodes ----------------------------------------------
    % --- Extruding nodes ----------------------------------------------
    % --- Extruding nodes ----------------------------------------------
    % --- Extruding nodes ----------------------------------------------
    % --- Extruding nodes ----------------------------------------------
    % --- Extruding nodes ----------------------------------------------
    % --- Extruding nodes ----------------------------------------------
    % --- Extruding nodes ----------------------------------------------
    % --- Extruding nodes ----------------------------------------------
    % --- Extruding nodes ----------------------------------------------
    % --- Extruding nodes ----------------------------------------------
    
elseif node<=nTotal
    
    nodes_adjacent=NaN;
    
    [xth,yth,zth]=get_root_of_extrusion_sph(node,n_nSurf,n,nx,ny,nz,nx_);
    node_=get_node_3d(xth,yth,zth,nx,ny);
    if xth==1
        if yth==1
            if zth==1
                csurf=zeros(1,7);
                csurf(1)=1;
                csurf(2)=2;
                csurf(3)=nx+1;
                csurf(4)=nx+2;
                csurf(5)=nxny+1;
                csurf(6)=nxny+2;
                csurf(7)=nxny+nx+1;
                extruding=get_extrusion_sph(csurf,n,nx,ny,nz,n_nSurf,nx_);
                extruding=[csurf,extruding];
            elseif zth==nz
                csurf=zeros(1,7);
                z1=(nz-2)*nxny+1;
                csurf(1)=z1;
                csurf(2)=z1+1;
                csurf(3)=z1+nx;
                z2=(nz-1)*nxny+1;
                csurf(4)=z2;
                csurf(5)=z2+1;
                csurf(6)=z2+nx;
                csurf(7)=csurf(6)+1;
                extruding=get_extrusion_sph(csurf,n,nx,ny,nz,n_nSurf,nx_);
                extruding=[csurf,extruding];
            else
                csurf=zeros(1,9);
                z2=get_node_3d(xth,yth,zth,nx,ny);
                z1=z2-nxny;
                z3=z2+nxny;
                csurf(1)=z1;
                csurf(2)=z1+1;
                csurf(3)=z1+nx;
                csurf(4)=z2;
                csurf(5)=z2+1;
                csurf(6)=z2+nx;
                csurf(7)=z3;
                csurf(8)=z3+1;
                csurf(9)=z3+nx;
                extruding=get_extrusion_sph(csurf,n,nx,ny,nz,n_nSurf,nx_);
                extruding=[csurf,extruding];
            end
        elseif yth==ny
            if zth==1
                csurf=zeros(1,7);
                z1=nxny-2*nx+1;
                z2=2*nxny-nx+1;
                csurf(1)=z1;
                csurf(2)=z1+1;
                csurf(3)=z1+nx;
                csurf(4)=csurf(3)+1;
                csurf(5)=z2-nx;
                csurf(6)=z2;
                csurf(7)=z2+1;
                extruding=get_extrusion_sph(csurf,n,nx,ny,nz,n_nSurf,nx_);
                extruding=[csurf,extruding];
            elseif zth==nz
                csurf=zeros(1,7);
                z2=get_node_3d(xth,yth,zth,nx,ny);
                z1=z2-nxny;
                csurf(1)=z1-nx;
                csurf(2)=z1;
                csurf(3)=z1+1;
                csurf(4)=z2-nx;
                csurf(5)=csurf(4)+1;
                csurf(6)=z2;
                csurf(7)=z2+1;
                extruding=get_extrusion_sph(csurf,n,nx,ny,nz,n_nSurf,nx_);
                extruding=[csurf,extruding];
            else
                csurf=zeros(1,9);
                z2=get_node_3d(xth,yth,zth,nx,ny);
                z1=z2-nxny;
                z3=z2+nxny;
                csurf(1)=z1-nx;
                csurf(2)=z1;
                csurf(3)=z1+1;
                csurf(4)=z2-nx;
                csurf(5)=z2;
                csurf(6)=z2+1;
                csurf(7)=z3-nx;
                csurf(8)=z3;
                csurf(9)=z3+1;
                extruding=get_extrusion_sph(csurf,n,nx,ny,nz,n_nSurf,nx_);
                extruding=[csurf,extruding];
            end
        else
            if zth==1
                csurf=zeros(1,9);
                csurf(1)=node_-nx;
                csurf(2)=csurf(1)+1;
                csurf(3)=node_;
                csurf(4)=node_+1;
                csurf(5)=node_+nx;
                csurf(6)=csurf(5)+1;
                csurf(8)=node_+nxny;
                csurf(7)=csurf(8)-nx;
                csurf(9)=csurf(8)+nx;
                extruding=get_extrusion_sph(csurf,n,nx,ny,nz,n_nSurf,nx_);
                extruding=[csurf,extruding];
            elseif zth==nz
                csurf=zeros(1,9);
                z1=node_-nxny;
                csurf(1)=z1-nx;
                csurf(2)=z1;
                csurf(3)=z1+nx;
                csurf(4)=node_-nx;
                csurf(5)=node_;
                csurf(6)=node_+nx;
                z1=node_+1;
                csurf(7)=z1-nx;
                csurf(8)=z1;
                csurf(9)=z1+nx;
                extruding=get_extrusion_sph(csurf,n,nx,ny,nz,n_nSurf,nx_);
                extruding=[csurf,extruding];
            else
                csurf=zeros(1,9);
                z2=node_;
                z1=node_-nxny;
                z3=node_+nxny;
                csurf(1)=z1-nx;
                csurf(2)=z1;
                csurf(3)=z1+nx;
                csurf(4)=z2-nx;
                csurf(5)=z2;
                csurf(6)=z2+nx;
                csurf(7)=z3-nx;
                csurf(8)=z3;
                csurf(9)=z3+nx;
                extruding=get_extrusion_sph(csurf,n,nx,ny,nz,n_nSurf,nx_);
                extruding=[csurf,extruding];
            end
        end
    elseif xth==nx
        if yth==1
            if zth==1
                csurf=zeros(1,7);
                csurf(1)=node_+nx-1;
                csurf(2)=node_-1;
                csurf(3)=node_;
                csurf(4)=node_+nx;
                csurf(6)=node_+nxny;
                csurf(5)=csurf(6)-1;
                csurf(7)=csurf(6)+nx;
                extruding=get_extrusion_sph(csurf,n,nx,ny,nz,n_nSurf,nx_);
                extruding=[csurf,extruding];
            elseif zth==nz
                csurf=zeros(1,7);
                csurf(2)=node_-nxny;
                csurf(1)=csurf(2)-1;
                csurf(3)=csurf(2)+nx;
                csurf(4)=node_-1;
                csurf(5)=node_;
                csurf(7)=node_+nx;
                csurf(6)=csurf(7)-1;
                extruding=get_extrusion_sph(csurf,n,nx,ny,nz,n_nSurf,nx_);
                extruding=[csurf,extruding];
            else
                csurf=zeros(1,9);
                z1=node_-nxny;
                z2=node_;
                z3=node_+nxny;
                csurf(1)=z1-1;
                csurf(2)=z1;
                csurf(3)=z1+nx;
                csurf(4)=z2-1;
                csurf(5)=z2;
                csurf(6)=z2+nx;
                csurf(7)=z3-1;
                csurf(8)=z3;
                csurf(9)=z3+nx;
                extruding=get_extrusion_sph(csurf,n,nx,ny,nz,n_nSurf,nx_);
                extruding=[csurf,extruding];
            end
        elseif yth==ny
            if zth==1
                csurf=zeros(1,7);
                csurf(2)=node_-nx;
                csurf(1)=csurf(2)-1;
                csurf(3)=node_-1;
                csurf(4)=node_;
                csurf(7)=node_+nxny;
                csurf(5)=csurf(7)-nx;
                csurf(6)=csurf(7)-1;
                extruding=get_extrusion_sph(csurf,n,nx,ny,nz,n_nSurf,nx_);
                extruding=[csurf,extruding];
            elseif zth==nz
                csurf=zeros(1,7);
                csurf(3)=node_-nxny;
                csurf(1)=csurf(3)-nx;
                csurf(2)=csurf(3)-1;
                csurf(5)=node_-nx;
                csurf(4)=csurf(5)-1;
                csurf(6)=node_-1;
                csurf(7)=node_;
                extruding=get_extrusion_sph(csurf,n,nx,ny,nz,n_nSurf,nx_);
                extruding=[csurf,extruding];
            else
                csurf=zeros(1,9);
                z1=node_-nxny;
                z2=node_;
                z3=node_+nxny;
                csurf(1)=z1-nx;
                csurf(2)=z1-1;
                csurf(3)=z1;
                csurf(4)=z2-nx;
                csurf(5)=z2-1;
                csurf(6)=z2;
                csurf(7)=z3-nx;
                csurf(8)=z3-1;
                csurf(9)=z3;
                extruding=get_extrusion_sph(csurf,n,nx,ny,nz,n_nSurf,nx_);
                extruding=[csurf,extruding];
            end
        else
            if zth==1
                csurf=zeros(1,9);
                csurf(2)=node_-nx;
                csurf(1)=csurf(2)-1;
                csurf(3)=node_-1;
                csurf(4)=node_;
                csurf(6)=node_+nx;
                csurf(5)=csurf(6)-1;
                csurf(8)=node_+nxny;
                csurf(7)=csurf(8)-nx;
                csurf(9)=csurf(8)+nx;
                extruding=get_extrusion_sph(csurf,n,nx,ny,nz,n_nSurf,nx_);
                extruding=[csurf,extruding];
            elseif zth==nz
                csurf=zeros(1,9);
                csurf(2)=node_-nxny;
                csurf(1)=csurf(2)-nx;
                csurf(3)=csurf(2)+nx;
                csurf(5)=node_-nx;
                csurf(4)=csurf(5)-1;
                csurf(6)=node_-1;
                csurf(7)=node_;
                csurf(9)=node_+nx;
                csurf(8)=csurf(9)-1;
                extruding=get_extrusion_sph(csurf,n,nx,ny,nz,n_nSurf,nx_);
                extruding=[csurf,extruding];
            else
                z1=node_-nxny;
                z2=node_;
                z3=node_+nxny;
                csurf(1)=z1-nx;
                csurf(2)=z1;
                csurf(3)=z1+nx;
                csurf(4)=z2-nx;
                csurf(5)=z2;
                csurf(6)=z2+nx;
                csurf(7)=z3-nx;
                csurf(8)=z3;
                csurf(9)=z3+nx;
                extruding=get_extrusion_sph(csurf,n,nx,ny,nz,n_nSurf,nx_);
                extruding=[csurf,extruding];
            end
        end
    else
        if yth==1
            if zth==1
                csurf=zeros(1,9);
                csurf(1)=node_-1;
                csurf(2)=node_;
                csurf(3)=node_+1;
                z1=node_+nx;
                csurf(4)=z1-1;
                csurf(5)=z1;
                csurf(6)=z1+1;
                z2=node_+nxny;
                csurf(7)=z2-1;
                csurf(8)=z2;
                csurf(9)=z2+1;
                extruding=get_extrusion_sph(csurf,n,nx,ny,nz,n_nSurf,nx_);
                extruding=[csurf,extruding];
            elseif zth==nz
                csurf=zeros(1,9);
                z1=node_-nxny;
                z2=node_;
                z3=node_+nx;
                csurf(1)=z1-1;
                csurf(2)=z1;
                csurf(3)=z1+1;
                csurf(4)=z2-1;
                csurf(5)=z2;
                csurf(6)=z2+1;
                csurf(7)=z3-1;
                csurf(8)=z3;
                csurf(9)=z3+1;
                extruding=get_extrusion_sph(csurf,n,nx,ny,nz,n_nSurf,nx_);
                extruding=[csurf,extruding];
            else
                csurf=zeros(1,9);
                z1=node_-nxny;
                z2=node_;
                z3=node_+nxny;
                csurf(1)=z1-1;
                csurf(2)=z1;
                csurf(3)=z1+1;
                csurf(4)=z2-1;
                csurf(5)=z2;
                csurf(6)=z2+1;
                csurf(7)=z3-1;
                csurf(8)=z3;
                csurf(9)=z3+1;
                extruding=get_extrusion_sph(csurf,n,nx,ny,nz,n_nSurf,nx_);
                extruding=[csurf,extruding];
            end
        elseif yth==ny
            if zth==1
                csurf=zeros(1,9);
                z1=node_-nx;
                z2=node_;
                z3=node_+nxny;
                csurf(1)=z1-1;
                csurf(2)=z1;
                csurf(3)=z1+1;
                csurf(4)=z2-1;
                csurf(5)=z2;
                csurf(6)=z2+1;
                csurf(7)=z3-1;
                csurf(8)=z3;
                csurf(9)=z3+1;
                extruding=get_extrusion_sph(csurf,n,nx,ny,nz,n_nSurf,nx_);
                extruding=[csurf,extruding];
            elseif zth==nz
                csurf=zeros(1,9);
                z1=node_-nxny;
                z2=node_-nx;
                z3=node_;
                csurf(1)=z1-1;
                csurf(2)=z1;
                csurf(3)=z1+1;
                csurf(4)=z2-1;
                csurf(5)=z2;
                csurf(6)=z2+1;
                csurf(7)=z3-1;
                csurf(8)=z3;
                csurf(9)=z3+1;
                extruding=get_extrusion_sph(csurf,n,nx,ny,nz,n_nSurf,nx_);
                extruding=[csurf,extruding];
            else
                csurf=zeros(1,9);
                z1=node_-nxny;
                z2=node_;
                z3=node_+nxny;
                csurf(1)=z1-1;
                csurf(2)=z1;
                csurf(3)=z1+1;
                csurf(4)=z2-1;
                csurf(5)=z2;
                csurf(6)=z2+1;
                csurf(7)=z3-1;
                csurf(8)=z3;
                csurf(9)=z3+1;
                extruding=get_extrusion_sph(csurf,n,nx,ny,nz,n_nSurf,nx_);
                extruding=[csurf,extruding];
            end
        else
            csurf=zeros(1,9);
            z1=node_-nx;
            z2=node_;
            z3=node_+nx;
            csurf(1)=z1-1;
            csurf(2)=z1;
            csurf(3)=z1+1;
            csurf(4)=z2-1;
            csurf(5)=z2;
            csurf(6)=z2+1;
            csurf(7)=z3-1;
            csurf(8)=z3;
            csurf(9)=z3+1;
            extruding=get_extrusion_sph(csurf,n,nx,ny,nz,n_nSurf,nx_);
            extruding=[csurf,extruding];
        end
    end
end
end