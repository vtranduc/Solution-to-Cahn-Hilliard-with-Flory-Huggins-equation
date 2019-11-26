function co=ic_sphere(ci,include_fluc,ci_fluc,n,radius,neight,nx,ny,nz,nex,ney,nez)

co=zeros(1,neight);

if include_fluc==1
    index=-7;
    dx=1/nex;
    dy=1/ney;
    dz=1/nez;
    z=-dz;
    summation=0;
    for zth=1:1:nz
        z=z+dz;
        y=-dy;
        for yth=1:1:ny
            y=y+dy;
            x=-dx;
            for xth=1:1:nx
                x=x+dx;
                index=index+8;
                dist=sqrt((0.5-x)^2+(0.5-y)^2+(0.5-z)^2);
                if dist<=radius %&& dist>radius*0.8
%                     co(index)=ci_fluc*(2*rand(1,1)-1);
                    co(index)=-ci_fluc;
                    summation=summation+co(index);
                end
            end
        end
    end
    adjustment=ci-summation/n;
    for i=1:8:neight-7
        co(i)=co(i)+adjustment;
    end
else
    for i=1:8:neight-7
        co(i)=co(i)+ci;
    end
end


end