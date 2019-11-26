function [ac_max,ac_min,ac_ave]=test81(c,ne,nex,ney,ny,weights,diffT,Two_chi_n1,n1,n2)
nex3=nex*3;ney3=ney*3;
allen_cahn=zeros(ney3,nex3);
gp=[(-sqrt(3/5)+1)/2 0.5 (sqrt(3/5)+1)/2];
dx=1/nex;dy=1/ney;
x=zeros(1,nex3);y=zeros(1,ney3);
gpx=gp*dx;
index=0;
for iex=1:1:nex
    z1=(iex-1)*dx;
    for i=1:1:3
        index=index+1;
        x(index)=z1+gpx(i);
    end
end
gpy=gp*dy;
index=0;
for iey=1:1:ney
    z1=(iey-1)*dy;
    for i=1:1:3
        index=index+1;
        y(index)=z1+gpy(i);
    end
end
conc_=get_conc(ne,ney,ny,c,weights);
ino=-3;
inx=0;
for e=1:1:ne
    if mod(e,ney)==1
        imo=0;
        ino=ino+3;
        inx=inx+1;
    end
    for m=1:1:3
        for n=1:1:3
            allen_cahn(imo+m,ino+n)=conc_(m,n,4,e)+conc_(m,n,5,e)...
                -diffT(inx,n)*((log(conc_(m,n,1,e))+1)/n1...
                -Two_chi_n1(inx,n)*conc_(m,n,1,e)...
                +Two_chi_n1(inx,n)/2-(log(1-conc_(m,n,1,e))+1)/n2);
        end
    end
    imo=imo+3;
end
contourf(x,y,allen_cahn,'edgecolor','none')
axis([0 1 0 1])
xlabel('x');ylabel('y')
% colorbar
ac_max=max(max(allen_cahn));
ac_min=min(min(allen_cahn));
ac_ave=mean(mean(allen_cahn));
end