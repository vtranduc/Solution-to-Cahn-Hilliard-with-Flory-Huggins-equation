function solution=sf_oneElement(c_,co_,weights,ilocal,n1,n2,dx,dy,dt,chi,diffT)

w=[5/18 4/9 5/18];

con=conc(c_,weights,1);
cono=conc(co_,weights,1);
conx=conc(c_,weights,2);
cony=conc(c_,weights,3);
conxx=conc(c_,weights,4);
conyy=conc(c_,weights,5);

cont=(con-cono)/dt;

phi=weights(:,:,ilocal,1);
phixx=weights(:,:,ilocal,4);
phiyy=weights(:,:,ilocal,5);

solution=0;
for ix=1:1:3
    for iy=1:1:3
        solution=solution+w(ix)*w(iy)*(...
            phi(ix,iy)*cont(ix,iy)-diffT(ix)*phi(ix,iy)*...
            ((-1/(con(ix,iy)^2*n1)+1/((1-con(ix,iy))^2*n2))*...
            (conx(ix,iy)^2+cony(ix,iy)^2)+...
            (1/(con(ix,iy)*n1)+1/((1-con(ix,iy))*n2)-2*chi(ix))*...
            (conxx(ix,iy)+conyy(ix,iy)))+...
            (conxx(ix,iy)+conyy(ix,iy))*...
            (phixx(ix,iy)+phiyy(ix,iy))...
            );

    end
end
solution=dx*dy*solution;
end