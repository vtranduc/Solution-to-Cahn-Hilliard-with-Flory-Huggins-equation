function sols=get_extrusion_sph(nodes,n,nx,ny,nz,n_nSurf,nx_)
len=length(nodes);
sols=zeros(1,len);
for i=1:1:len
    node=nodes(i);
    [xth,yth,zth]=get_xyzth_3d(node,nx,ny);
    if xth==1
        sol=ny*(zth-1)+yth;
    elseif xth==nx
        sol=n_nSurf(1)+ny*(zth-1)+yth;
    else
        if yth==1
            sol=n_nSurf(2)+nx_*(zth-1)+xth-1;
        elseif yth==ny
            sol=n_nSurf(3)+nx_*(zth-1)+xth-1;
        else
            if zth==1
                sol=n_nSurf(4)+nx_*(yth-2)+xth-1;
            elseif zth==nz
                sol=n_nSurf(5)+nx_*(yth-2)+xth-1;
            end
        end
    end
    sols(i)=sol+n;
end
end