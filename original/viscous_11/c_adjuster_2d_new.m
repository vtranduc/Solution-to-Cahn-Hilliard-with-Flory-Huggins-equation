function [c_adjusted,range_check]=c_adjuster_2d_new(...
    ne,ney,dxdy,ny,c,weights,ci_ave,xlen,ylen,nfour)
w=[5/18 4/9 5/18];
w_set=zeros(3,3);
for m=1:1:3
    for n=1:1:3
        w_set(m,n)=w(m)*w(n);
    end
end
conc_abs=get_conc_abs_new(ne,ney,ny,c,weights);
c_ave_current=0;
for e=1:1:ne
    for iny=1:1:3
        for inx=1:1:3
            c_ave_current=c_ave_current+w_set(inx,iny)*(conc_abs(inx,iny,e));
        end
    end
end
c_ave_current=dxdy*c_ave_current;
deviation=(c_ave_current-ci_ave)*xlen*ylen;
c_adjusted=c;
for i=1:4:nfour-3
    c_adjusted(i)=c(i)-deviation;
    if c_adjusted(i)<=0 || c_adjusted(i)>=1
        range_check=0;
        return
    end
end
range_check=1;
end

function conc_abs=get_conc_abs_new(ne,ney,ny,c,weights)
conc_abs=zeros(3,3,ne);
for element=1:1:ne
    [inx,iny]=inxiny_elemental(element,ney);
    gbfs=elemental_gbf(inx,iny,ny);
    cElemental=get_cElemental(gbfs,c);
    conc_abs(:,:,element)=conc(cElemental,weights,1);
end
end