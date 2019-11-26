function interactions=interactions_rec_3D(node,ny,nz,nxny)

%THIS IS INVALID FOR NX/NY/NZ=1!!!!

floor_ID=mod(node,nxny);
if floor_ID==0
    floor_ID=nxny;
end

z_=(node-floor_ID)/nxny+1;

if z_==1
    z_=-1;
elseif z_==nz
    z_=1;
else
    z_=0;
end

if floor_ID<=ny
    x_=-1;
elseif floor_ID>nxny-ny
    x_=1;
else
    x_=0;
end

y_=mod(floor_ID,ny);
if y_==1
    y_=-1;
elseif y_==0
    y_=1;
else
    y_=0;
end

if z_==0
    if y_==0
        if x_==0
            interactions=zeros(1,27);
            interactions(1)=node-nxny-ny-1;
            interactions(2)=interactions(1)+1;
            interactions(3)=interactions(1)+2;
            interactions(4)=interactions(1)+ny;
            interactions(5)=interactions(4)+1;
            interactions(6)=interactions(4)+2;
            interactions(7)=interactions(4)+ny;
            interactions(8)=interactions(7)+1;
            interactions(9)=interactions(7)+2;
            i1=9;
            i2=18;
            for i=1:1:9
                i1=i1+1;
                i2=i2+1;
                interactions(i1)=interactions(i)+nxny;
                interactions(i2)=interactions(i)+2*nxny;
            end
        elseif x_==-1
            interactions=zeros(1,18);
            interactions(1)=node-nxny-1;
            interactions(2)=interactions(1)+1;
            interactions(3)=interactions(1)+2;
            interactions(4)=interactions(1)+ny;
            interactions(5)=interactions(4)+1;
            interactions(6)=interactions(4)+2;
            i1=6;
            i2=12;
            for i=1:1:6
                i1=i1+1;
                i2=i2+1;
                interactions(i1)=interactions(i)+nxny;
                interactions(i2)=interactions(i)+2*nxny;
            end
        elseif x_==1
            interactions=zeros(1,18);
            interactions(1)=node-nxny-ny-1;
            interactions(2)=interactions(1)+1;
            interactions(3)=interactions(1)+2;
            interactions(4)=interactions(1)+ny;
            interactions(5)=interactions(4)+1;
            interactions(6)=interactions(4)+2;
            i1=6;
            i2=12;
            for i=1:1:6
                i1=i1+1;
                i2=i2+1;
                interactions(i1)=interactions(i)+nxny;
                interactions(i2)=interactions(i)+2*nxny;
            end
        end
    elseif y_==-1
        if x_==0
            interactions=zeros(1,18);
            interactions(1)=node-nxny-ny;
            interactions(2)=interactions(1)+1;
            interactions(3)=interactions(1)+ny;
            interactions(4)=interactions(3)+1;
            interactions(5)=interactions(3)+ny;
            interactions(6)=interactions(5)+1;
            i1=6;
            i2=12;
            for i=1:1:6
                i1=i1+1;
                i2=i2+1;
                interactions(i1)=interactions(i)+nxny;
                interactions(i2)=interactions(i)+2*nxny;
            end
        elseif x_==-1
            interactions=zeros(1,12);
            interactions(1)=node-nxny;
            interactions(2)=interactions(1)+1;
            interactions(3)=interactions(1)+ny;
            interactions(4)=interactions(3)+1;
            i1=4;
            i2=8;
            for i=1:1:4
                i1=i1+1;
                i2=i2+1;
                interactions(i1)=interactions(i)+nxny;
                interactions(i2)=interactions(i)+2*nxny;
            end
        elseif x_==1
            interactions=zeros(1,12);
            interactions(1)=node-nxny-ny;
            interactions(2)=interactions(1)+1;
            interactions(3)=interactions(1)+ny;
            interactions(4)=interactions(3)+1;
            i1=4;
            i2=8;
            for i=1:1:4
                i1=i1+1;
                i2=i2+1;
                interactions(i1)=interactions(i)+nxny;
                interactions(i2)=interactions(i)+2*nxny;
            end
        end
    elseif y_==1
        if x_==0
            interactions=zeros(1,18);
            interactions(1)=node-nxny-ny-1;
            interactions(2)=interactions(1)+1;
            interactions(3)=interactions(1)+ny;
            interactions(4)=interactions(3)+1;
            interactions(5)=interactions(3)+ny;
            interactions(6)=interactions(5)+1;
            i1=6;
            i2=12;
            for i=1:1:6
                i1=i1+1;
                i2=i2+1;
                interactions(i1)=interactions(i)+nxny;
                interactions(i2)=interactions(i)+2*nxny;
            end
        elseif x_==-1
            interactions=zeros(1,12);
            interactions(1)=node-nxny-1;
            interactions(2)=interactions(1)+1;
            interactions(3)=interactions(1)+ny;
            interactions(4)=interactions(3)+1;
            i1=4;
            i2=8;
            for i=1:1:4
                i1=i1+1;
                i2=i2+1;
                interactions(i1)=interactions(i)+nxny;
                interactions(i2)=interactions(i)+2*nxny;
            end
        elseif x_==1
            interactions=zeros(1,12);
            interactions(1)=node-nxny-ny-1;
            interactions(2)=interactions(1)+1;
            interactions(3)=interactions(1)+ny;
            interactions(4)=interactions(3)+1;
            i1=4;
            i2=8;
            for i=1:1:4
                i1=i1+1;
                i2=i2+1;
                interactions(i1)=interactions(i)+nxny;
                interactions(i2)=interactions(i)+2*nxny;
            end
        end
    end
elseif z_==-1
    if y_==0
        if x_==0
            interactions=zeros(1,18);
            interactions(1)=node-ny-1;
            interactions(2)=interactions(1)+1;
            interactions(3)=interactions(1)+2;
            interactions(4)=node-1;
            interactions(5)=node;
            interactions(6)=node+1;
            interactions(7)=interactions(4)+ny;
            interactions(8)=interactions(7)+1;
            interactions(9)=interactions(8)+1;
            i1=9;
            for i=1:1:9
                i1=i1+1;
                interactions(i1)=interactions(i)+nxny;
            end
        elseif x_==-1
            interactions=zeros(1,12);
            interactions(1)=node-1;
            interactions(2)=node;
            interactions(3)=node+1;
            for i=1:1:3
                interactions(i+3)=interactions(i)+ny;
            end
            for i=1:1:6
                interactions(i+6)=interactions(i)+nxny;
            end
        elseif x_==1
            interactions=zeros(1,12);
            interactions(1)=node-ny-1;
            interactions(2)=interactions(1)+1;
            interactions(3)=interactions(1)+2;
            for i=1:1:3
                interactions(i+3)=interactions(i)+ny;
            end
            for i=1:1:6
                interactions(i+6)=interactions(i)+nxny;
            end
        end
    elseif y_==-1
        if x_==0
            interactions=zeros(1,12);
            interactions(1)=node-ny;
            interactions(2)=interactions(1)+1;
            interactions(3)=node;
            interactions(4)=node+1;
            interactions(5)=node+ny;
            interactions(6)=interactions(5)+1;
            for i=1:1:6
                interactions(i+6)=interactions(i)+nxny;
            end
        elseif x_==-1
            interactions=zeros(1,8);
            interactions(1)=1;
            interactions(2)=2;
            interactions(3)=1+ny;
            interactions(4)=interactions(3)+1;
            for i=1:1:4
                interactions(i+4)=interactions(i)+nxny;
            end
        elseif x_==1
            interactions=zeros(1,8);
            interactions(1)=node-ny;
            interactions(2)=interactions(1)+1;
            interactions(3)=node;
            interactions(4)=node+1;
            for i=1:1:4
                interactions(i+4)=interactions(i)+nxny;
            end
        end
    elseif y_==1
        if x_==0
            interactions=zeros(1,12);
            interactions(1)=node-ny-1;
            interactions(2)=node-ny;
            interactions(3)=node-1;
            interactions(4)=node;
            interactions(5)=interactions(3)+ny;
            interactions(6)=node+ny;
            for i=1:1:6
                interactions(i+6)=interactions(i)+nxny;
            end
        elseif x_==-1
            interactions=zeros(1,8);
            interactions(1)=node-1;
            interactions(2)=node;
            interactions(3)=interactions(1)+ny;
            interactions(4)=node+ny;
            for i=1:1:4
                interactions(i+4)=interactions(i)+nxny;
            end
        elseif x_==1
            interactions=zeros(1,8);
            interactions(1)=node-ny-1;
            interactions(2)=interactions(1)+1;
            interactions(3)=node-1;
            interactions(4)=node;
            for i=1:1:4
                interactions(i+4)=interactions(i)+nxny;
            end
        end
    end
elseif z_==1
    if y_==0
        if x_==0
            interactions=zeros(1,18);
            interactions(1)=node-nxny-ny-1;
            interactions(2)=interactions(1)+1;
            interactions(3)=interactions(1)+2;
            interactions(4)=interactions(1)+ny;
            interactions(5)=interactions(4)+1;
            interactions(6)=interactions(4)+2;
            interactions(7)=interactions(4)+ny;
            interactions(8)=interactions(7)+1;
            interactions(9)=interactions(7)+2;
            i1=9;
            for i=1:1:9
                i1=i1+1;
                interactions(i1)=interactions(i)+nxny;
            end
        elseif x_==-1
            interactions=zeros(1,12);
            interactions(1)=node-nxny-1;
            interactions(2)=interactions(1)+1;
            interactions(3)=interactions(1)+2;
            interactions(4)=interactions(1)+ny;
            interactions(5)=interactions(4)+1;
            interactions(6)=interactions(5)+1;
            for i=1:1:6
                interactions(i+6)=interactions(i)+nxny;
            end
        elseif x_==1
            interactions=zeros(1,12);
            interactions(1)=node-nxny-ny-1;
            interactions(2)=interactions(1)+1;
            interactions(3)=interactions(1)+2;
            interactions(4)=interactions(1)+ny;
            interactions(5)=interactions(4)+1;
            interactions(6)=interactions(5)+1;
            for i=1:1:6
                interactions(i+6)=interactions(i)+nxny;
            end
        end
    elseif y_==-1
        if x_==0
            interactions=zeros(1,12);
            interactions(1)=node-nxny-ny;
            interactions(2)=interactions(1)+1;
            interactions(3)=node-nxny;
            interactions(4)=interactions(3)+1;
            interactions(5)=interactions(3)+ny;
            interactions(6)=interactions(5)+1;
            for i=1:1:6
                interactions(i+6)=interactions(i)+nxny;
            end
        elseif x_==-1
            interactions=zeros(1,8);
            interactions(1)=node-nxny;
            interactions(2)=interactions(1)+1;
            interactions(3)=interactions(1)+ny;
            interactions(4)=interactions(3)+1;
            for i=1:1:4
                interactions(i+4)=interactions(i)+nxny;
            end
        elseif x_==1
            interactions=zeros(1,8);
            interactions(1)=node-nxny-ny;
            interactions(2)=interactions(1)+1;
            interactions(3)=interactions(1)+ny;
            interactions(4)=interactions(3)+1;
            for i=1:1:4
                interactions(i+4)=interactions(i)+nxny;
            end
        end
    elseif y_==1
        if x_==0
            interactions=zeros(1,12);
            interactions(1)=node-nxny-ny-1;
            interactions(2)=interactions(1)+1;
            interactions(3)=interactions(1)+ny;
            interactions(4)=interactions(3)+1;
            interactions(5)=interactions(3)+ny;
            interactions(6)=interactions(5)+1;
            for i=1:1:6
                interactions(i+6)=interactions(i)+nxny;
            end
        elseif x_==-1
            interactions=zeros(1,8);
            interactions(1)=node-nxny-1;
            interactions(2)=interactions(1)+1;
            interactions(3)=interactions(1)+ny;
            interactions(4)=interactions(3)+1;
            for i=1:1:4
                interactions(i+4)=interactions(i)+nxny;
            end
        elseif x_==1
            interactions=zeros(1,8);
            interactions(1)=node-nxny-ny-1;
            interactions(2)=interactions(1)+1;
            interactions(3)=interactions(1)+ny;
            interactions(4)=interactions(3)+1;
            for i=1:1:4
                interactions(i+4)=interactions(i)+nxny;
            end
        end
    end
end

end