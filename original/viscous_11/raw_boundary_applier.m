function [sf,sj]=raw_boundary_applier(sf_,sj_,irow,sf_len)
sf=sf_;
sj=sj_;
sf(irow)=0;
for i=1:1:sf_len
    sj(irow,i)=0;
end
sj(irow,irow)=1;
end