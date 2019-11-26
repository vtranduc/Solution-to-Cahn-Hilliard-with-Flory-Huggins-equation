function nodes=get_nodes_of_element_sph(e,ne,nex,ney,nexney,nx,ny,nz,n,nxny,n_eSurf,n_nSurf,nx_)

if e<=ne
    nodes=zeros(1,8);
    [n1xth,n1yth,n1zth]=get_n1xyzth_3d(e,nex,nexney);
    nodes(1)=(n1zth-1)*nxny+nx*(n1yth-1)+n1xth;
    nodes(2)=nodes(1)+1;
    nodes(3)=nodes(1)+nx;
    nodes(4)=nodes(3)+1;
    nodes(5)=nodes(1)+nxny;
    nodes(6)=nodes(2)+nxny;
    nodes(7)=nodes(3)+nxny;
    nodes(8)=nodes(4)+nxny;
else
    esurf=e-ne;
    if esurf<=n_eSurf(1)
        nodes=get_nodes_es1(esurf,ney,n,nx,ny,nxny);
    elseif esurf<=n_eSurf(2)
        nodes=get_nodes_es2(esurf,n_eSurf,n_nSurf,ney,n,nx,ny,nxny);
    elseif esurf<=n_eSurf(3)
        nodes=get_nodes_es3(esurf,n_eSurf,n_nSurf,nex,n,nx,ny,nxny,nx_);
    elseif esurf<=n_eSurf(4)
        nodes=get_nodes_es4(esurf,n_eSurf,n_nSurf,nex,n,nx,ny,nxny,nx_);
    elseif esurf<=n_eSurf(5)
        nodes=get_nodes_es5(esurf,n_eSurf,n_nSurf,nex,ney,n,nx,ny,nx_);
    elseif esurf<=n_eSurf(6)
        nodes=get_nodes_es6(esurf,n_eSurf,n_nSurf,nex,ney,n,nx,ny,nz,nx_);
    end
end

end

function nodes=get_nodes_es1(esurf,ney,n,nx,ny,nxny)
nodes=zeros(1,8);
yth=mod_sp(esurf,ney);
zth=(esurf-yth)/ney+1;
nodes(1)=ny*(zth-1)+yth+n;
nodes(2)=get_node_3d(1,yth,zth,nx,ny);
nodes(3)=nodes(1)+1;
nodes(4)=nodes(2)+nx;
nodes(5)=nodes(1)+ny;
nodes(6)=nodes(2)+nxny;
nodes(7)=nodes(5)+1;
nodes(8)=nodes(4)+nxny;
end

function nodes=get_nodes_es2(esurf_,n_eSurf,n_nSurf,ney,n,nx,ny,nxny)
nodes=zeros(1,8);
esurf=esurf_-n_eSurf(1);
yth=mod_sp(esurf,ney);
zth=(esurf-yth)/ney+1;
nodes(1)=get_node_3d(nx,yth,zth,nx,ny);
nodes(2)=ny*(zth-1)+yth+n_nSurf(1)+n;
nodes(3)=nodes(1)+nx;
nodes(4)=nodes(2)+1;
nodes(5)=nodes(1)+nxny;
nodes(6)=nodes(2)+ny;
nodes(7)=nodes(3)+nxny;
nodes(8)=nodes(6)+1;
end

function nodes=get_nodes_es3(esurf_,n_eSurf,n_nSurf,nex,n,nx,ny,nxny,nx_)
nodes=zeros(1,8);
esurf=esurf_-n_eSurf(2);
xth=mod_sp(esurf,nex);
zth=(esurf-xth)/nex+1;
if xth==1
    nodes(1)=ny*(zth-1)+1+n;
    nodes(2)=n_nSurf(2)+nx_*(zth-1)+xth+n;
    nodes(5)=nodes(1)+ny;
    nodes(6)=nodes(2)+nx_;
elseif xth==nex
    nodes(1)=n_nSurf(2)+nx_*(zth-1)+xth-1+n;
    nodes(2)=n_nSurf(1)+ny*(zth-1)+1+n;
    nodes(5)=nodes(1)+nx_;
    nodes(6)=nodes(2)+ny;
else
    nodes(1)=n_nSurf(2)+nx_*(zth-1)+xth-1+n;
    nodes(2)=nodes(1)+1;
    nodes(5)=nodes(1)+nx_;
    nodes(6)=nodes(5)+1;
end
nodes(3)=get_node_3d(xth,1,zth,nx,ny);
nodes(4)=nodes(3)+1;
nodes(7)=nodes(3)+nxny;
nodes(8)=nodes(4)+nxny;
end

function nodes=get_nodes_es4(esurf_,n_eSurf,n_nSurf,nex,n,nx,ny,nxny,nx_)
nodes=zeros(1,8);
esurf=esurf_-n_eSurf(3);
xth=mod_sp(esurf,nex);
zth=(esurf-xth)/nex+1;
if xth==1
    nodes(3)=ny*zth+n;
    nodes(4)=n_nSurf(3)+nx_*(zth-1)+xth+n;
    nodes(7)=nodes(3)+ny;
    nodes(8)=nodes(4)+nx_;
elseif xth==nex
    nodes(3)=n_nSurf(3)+nx_*(zth-1)+nx_+n;
    nodes(4)=n_nSurf(1)+ny*zth+n;
    nodes(7)=nodes(3)+nx_;
    nodes(8)=nodes(4)+ny;
else
    nodes(3)=n_nSurf(3)+nx_*(zth-1)+xth-1+n;
    nodes(4)=nodes(3)+1;
    nodes(7)=nodes(3)+nx_;
    nodes(8)=nodes(7)+1;
end
nodes(1)=get_node_3d(xth,ny,zth,nx,ny);
nodes(2)=nodes(1)+1;
nodes(5)=nodes(1)+nxny;
nodes(6)=nodes(2)+nxny;
end

function nodes=get_nodes_es5(esurf_,n_eSurf,n_nSurf,nex,ney,n,nx,ny,nx_)
nodes=zeros(1,8);
esurf=esurf_-n_eSurf(4);
xth=mod_sp(esurf,nex);
yth=(esurf-xth)/nex+1;
if xth==1
    nodes(1)=yth+n;
    nodes(3)=nodes(1)+1;
    if yth==1
        nodes(2)=n_nSurf(2)+1+n;
        nodes(4)=n_nSurf(4)+1+n;
    elseif yth==ney
        nodes(2)=n_nSurf(5)-nx_+1+n;
        nodes(4)=n_nSurf(3)+1+n;
    else
        nodes(2)=n_nSurf(4)+nx_*(yth-2)+1+n;
        nodes(4)=nodes(2)+nx_;
    end
elseif xth==nex
    nodes(2)=n_nSurf(1)+yth+n;
    nodes(4)=nodes(2)+1;
    if yth==1
        nodes(1)=n_nSurf(2)+nx_+n;
        nodes(3)=n_nSurf(4)+nx_+n;
    elseif yth==ney
        nodes(1)=n_nSurf(5)+n;
        nodes(3)=n_nSurf(3)+nx_+n;
    else
        nodes(1)=n_nSurf(4)+nx_*(yth-1)+n;
        nodes(3)=nodes(1)+nx_;
    end
else
    if yth==1
        nodes(1)=n_nSurf(2)+xth-1+n;
        nodes(2)=nodes(1)+1;
        nodes(3)=n_nSurf(4)+xth-1+n;
        nodes(4)=nodes(3)+1;
    elseif yth==ney
        nodes(1)=n_nSurf(5)-nx+xth+1+n;
        nodes(2)=nodes(1)+1;
        nodes(3)=n_nSurf(3)+xth-1+n;
        nodes(4)=nodes(3)+1;
    else
        nodes(1)=nx_*(yth-2)+xth-1+n_nSurf(4)+n;
        nodes(2)=nodes(1)+1;
        nodes(3)=nodes(1)+nx_;
        nodes(4)=nodes(3)+1;
    end
end
nodes(5)=get_node_3d(xth,yth,1,nx,ny);
nodes(6)=nodes(5)+1;
nodes(7)=nodes(5)+nx;
nodes(8)=nodes(7)+1;
end

function nodes=get_nodes_es6(esurf_,n_eSurf,n_nSurf,nex,ney,n,nx,ny,nz,nx_)
nodes=zeros(1,8);
esurf=esurf_-n_eSurf(5);
xth=mod_sp(esurf,nex);
yth=(esurf-xth)/nex+1;
nodes(1)=get_node_3d(xth,yth,nz,nx,ny);
nodes(2)=nodes(1)+1;
nodes(3)=nodes(1)+nx;
nodes(4)=nodes(3)+1;
if xth==1
    if yth==1
        nodes(5)=n_nSurf(1)-ney+n;
        nodes(7)=nodes(5)+1;
        nodes(6)=n_nSurf(3)-nx_+1+n;
        nodes(8)=n_nSurf(5)+1+n;
    elseif yth==ney
        nodes(5)=n_nSurf(1)-1+n;
        nodes(7)=nodes(5)+1;
        nodes(6)=n_nSurf(6)-nx_+1+n;
        nodes(8)=n_nSurf(4)-nx_+1+n;
    else
        nodes(5)=n_nSurf(1)-ny+yth+n;
        nodes(7)=nodes(5)+1;
        nodes(6)=n_nSurf(5)+nx_*(yth-2)+xth+n;
        nodes(8)=nodes(6)+nx_;
    end
elseif xth==nex
    if yth==1
        nodes(5)=n_nSurf(3)+n;
        nodes(6)=n_nSurf(2)-ney+n;
        nodes(7)=n_nSurf(5)+nx_+n;
        nodes(8)=nodes(6)+1;
    elseif yth==ney
        nodes(5)=n_nSurf(6)+n;
        nodes(6)=n_nSurf(2)-1+n;
        nodes(7)=n_nSurf(4)+n;
        nodes(8)=n_nSurf(2)+n;
    else
        nodes(5)=n_nSurf(5)+nx_+nx_*(yth-2)+n;
        nodes(6)=n_nSurf(2)-ny+yth+n;
        nodes(7)=nodes(5)+nx_;
        nodes(8)=nodes(6)+1;
    end
else
    if yth==1
        nodes(5)=n_nSurf(3)-nx_+xth-1+n;
        nodes(6)=nodes(5)+1;
        nodes(7)=n_nSurf(5)+xth-1+n;
        nodes(8)=nodes(7)+1;
    elseif yth==ney
        nodes(5)=n_nSurf(6)-nex+xth+n;
        nodes(6)=nodes(5)+1;
        nodes(7)=n_nSurf(4)-nex+xth+n;
        nodes(8)=nodes(7)+1;
    else
        nodes(5)=n_nSurf(5)+nx_*(yth-2)+xth-1+n;
        nodes(6)=nodes(5)+1;
        nodes(7)=nodes(5)+nx_;
        nodes(8)=nodes(7)+1;
    end
end
end