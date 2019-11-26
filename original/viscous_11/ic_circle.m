function co=ic_circle(ci,include_fluc,ci_fluc,n,nx,ny,nfour,radius,nex,ney)

co=zeros(1,nfour);

if include_fluc==1
    index=-3;
    dx=1/nex;
    dy=1/ney;
    x=-dx;
    summation=0;
    for ix=1:1:nx
        x=x+dx;
        y=-dy;
        for iy=1:1:ny
            y=y+dy;
            index=index+4;
            dist=sqrt((x-0.5)^2+(y-0.5)^2);
            if dist<=radius
                co(index)=-ci_fluc;
                summation=summation+co(index);
%                 error('adfaga')
%             else
%                 warning('hey')
            end
        end
    end
else
    summation=0;
end

diff=ci-summation/n;

for i=1:4:nfour-3
%     if co(i)~=0
%         co(i);
% %         error('dfaga')
%     end
    co(i)=co(i)+diff;
    
end

end