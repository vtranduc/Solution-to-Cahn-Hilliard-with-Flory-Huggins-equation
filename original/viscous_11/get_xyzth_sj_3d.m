function [xth,yth,zth,nodeth,relationth]=get_xyzth_sj_3d(sjth,nx,ny,nz,...
    nnz_)

n1stlastzth=1152*(nx-2)*(ny-2)+1536*(nx+ny-4)+2048;
nbtwzth=1728*(nx-2)*(ny-2)+2304*(nx+ny-4)+3072;
%=============================
if sjth<=n1stlastzth
    zth=1;
    xyth=sjth;
elseif sjth>=nnz_-n1stlastzth+1
    zth=nz;
    xyth=sjth-((nz-2)*nbtwzth+n1stlastzth);
else
    a=sjth-n1stlastzth;
    xyth=mod(a,nbtwzth);
    if xyth==0
        xyth=nbtwzth;
    end
    zth=(a-xyth)/nbtwzth+2;
end
%=============================
if zth==1 || zth==nz
    n1stlastyth=96*nx-64;
    nbtwyth=144*nx-96;
else
    n1stlastyth=144*nx-96;
    nbtwyth=216*nx-144;
end
n1stlastyth=n1stlastyth*8;
nbtwyth=nbtwyth*8;
%=============================
if xyth<=n1stlastyth
    yth=1;
    xxyth=xyth;
elseif xyth>=n1stlastyth+(ny-2)*nbtwyth+1
    yth=ny;
    xxyth=xyth-((ny-2)*nbtwyth+n1stlastyth);
else
    b=xyth-n1stlastyth;
    xxyth=mod(b,nbtwyth);
    if xxyth==0
        xxyth=nbtwyth;
    end
    yth=(b-xxyth)/nbtwyth+2;
end
%=============================
if zth==1 || zth==nz
    if yth==1 || yth==ny
        n1stlastxth=64;
        nbtwxth=96;
    else
        n1stlastxth=96;
        nbtwxth=144;
    end
else
    if yth==1 || yth==ny
        n1stlastxth=96;
        nbtwxth=144;
    else
        n1stlastxth=144;
        nbtwxth=216;
    end
end
n1stlastxth=n1stlastxth*8;
nbtwxth=nbtwxth*8;
%=============================
if xxyth<=n1stlastxth
    xth=1;
    relationshipth=xxyth;
elseif xxyth >= n1stlastxth+(nx-2)*nbtwxth+1
    xth=nx;
    relationshipth=xxyth-((nx-2)*nbtwxth+n1stlastxth);
else
    c=xxyth-n1stlastxth;
    relationshipth=mod(c,nbtwxth);
    if relationshipth==0
        relationshipth=nbtwxth;
    end
    xth=(c-relationshipth)/nbtwxth+2;
end
%=============================
edge_degree=0;
if xth==1 || xth==nx
    edge_degree=edge_degree+1;
end
if yth==1 || yth==ny
    edge_degree=edge_degree+1;
end
if zth==1 || zth==nz
    edge_degree=edge_degree+1;
end
if edge_degree==0
    nInteractions=216;
elseif edge_degree==1
    nInteractions=144;
elseif edge_degree==2
    nInteractions=96;
elseif edge_degree==3
    nInteractions=64;
end
%=============================
relationth=mod(relationshipth,nInteractions);
if relationth==0
    relationth=nInteractions;
end
nodeth=(relationshipth-relationth)/nInteractions+1;
%=============================

%=============================
end