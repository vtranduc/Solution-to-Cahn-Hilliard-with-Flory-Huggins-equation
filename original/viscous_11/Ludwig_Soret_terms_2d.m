function terms=Ludwig_Soret_terms_2d(c,ne,two_slope)

terms=zeros(ne,3,3,3);

%==========================
D_=0;
% D_=10^5;
%==========================

if D_==0
    return
end

slope=two_slope/2;


for e=1:1:ne
    for ix=1:1:3
        for iy=1:1:3
            sub1=1-2*c(ix,iy,1,e);
            sub2=slope*c(ix,iy,2,e);
            terms(e,ix,iy,1)=sub1*sub2;
            terms(e,ix,iy,2)=-2*sub2;
            terms(e,ix,iy,3)=sub1*slope;
        end
    end
end
terms=terms*D_;
end