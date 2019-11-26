function cont=get_cont_curved_3d(e,dt,conc_,conco_)
cont=zeros(3,3,3);
for ix=1:1:3
    for iy=1:1:3
        for iz=1:1:3
            cont(ix,iy,iz)=conc_(e,1,ix,iy,iz)-conco_(e,ix,iy,iz);
        end
    end
end
cont=cont/dt;
end