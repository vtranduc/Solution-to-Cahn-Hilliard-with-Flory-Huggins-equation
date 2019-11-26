function [g_x_minus,g_y_minus,g_z_minus,g_x_plus,g_y_plus,g_z_plus,...
        h_x_minus,h_y_minus,h_z_minus,h_x_plus,h_y_plus,h_z_plus]=...
        generate_wetting_pattern(type,spec,hs,gs,nx,ny,nz,xlen,ylen,zlen)
    
% type==1 Stripe pattern
% type==2 Wetting on the top and bottom surface
% type==3 Square tile pattern
% type==4 Same wetting properties on all surfaces with respective normals
% type==5 Circular wetting on the top surface

if type==1
    if length(spec)~=1 || length(hs)~=2 || length(gs)~=2
        error('Specifications for type 1 is incorrect')
    end
    if spec>nx
        error('Too many patterns for such a small mesh!')
    end
    g_x_minus=zeros(ny,nz);
    g_y_minus=zeros(nx,nz);
    g_z_minus=zeros(nx,ny);
    g_x_plus=zeros(ny,nz);
    g_y_plus=zeros(nx,nz);
    g_z_plus=zeros(nx,ny);

    h_x_minus=zeros(ny,nz);
    h_y_minus=zeros(nx,nz);
    h_z_minus=zeros(nx,ny);
    h_x_plus=zeros(ny,nz);
    h_y_plus=zeros(nx,nz);
    h_z_plus=zeros(nx,ny);
    
    
    splitter=zeros(1,spec);
    splitter_assist=nx/spec;
    for i=1:1:spec-1
        splitter(i)=round(splitter_assist*i);
    end
    splitter(spec)=nx;
    switcher=1;
    xth=1;
    for i=1:1:spec
        switcher=mod(switcher+1,2);
        if switcher==0
            h=hs(1);
            g=gs(1);
        elseif switcher==1
            h=hs(2);
            g=gs(2);
        end
        for j=xth:1:splitter(i)
            
            for k=1:1:ny
                h_z_plus(j,k)=h;
                g_z_plus(j,k)=g;
            end
        end
        xth=splitter(i)+1;
    end
    
elseif type==2
    if length(hs)~=2 || length(gs)~=2
        error('Specifications for type 1 is incorrect')
    end
    g_x_minus=zeros(ny,nz);
    g_y_minus=zeros(nx,nz);
    g_z_minus=gs(1)*ones(nx,ny);
    g_x_plus=zeros(ny,nz);
    g_y_plus=zeros(nx,nz);
    g_z_plus=gs(2)*ones(nx,ny);

    h_x_minus=zeros(ny,nz);
    h_y_minus=zeros(nx,nz);
    h_z_minus=hs(1)*ones(nx,ny);
    h_x_plus=zeros(ny,nz);
    h_y_plus=zeros(nx,nz);
    h_z_plus=hs(2)*ones(nx,ny);
    
elseif type==3
    
    if length(spec)~=1 || length(hs)~=2 || length(gs)~=2
        error('Specifications for type 3 is incorrect')
    end
    if spec+1>nx || spec+1>ny
        error('Too many patterns for such a small mesh!')
    end
    g_x_minus=zeros(ny,nz);
    g_y_minus=zeros(nx,nz);
    g_z_minus=zeros(nx,ny);
    g_x_plus=zeros(ny,nz);
    g_y_plus=zeros(nx,nz);
    g_z_plus=zeros(nx,ny);

    h_x_minus=zeros(ny,nz);
    h_y_minus=zeros(nx,nz);
    h_z_minus=zeros(nx,ny);
    h_x_plus=zeros(ny,nz);
    h_y_plus=zeros(nx,nz);
    h_z_plus=zeros(nx,ny);
    
    xsplit=zeros(spec,2);
    ysplit=zeros(spec,2);
    
    dx=(nx-1)/spec;
    dy=(ny-1)/spec;
    
    xsplit(1,1)=1;
    ysplit(1,1)=1;
    
    xsplit(spec,2)=nx;
    ysplit(spec,2)=ny;
    
    for i=1:1:spec-1
        xsplit(i,2)=round(dx*i);
        xsplit(i+1,1)=xsplit(i,2)+1;
        ysplit(i,2)=round(dy*i);
        ysplit(i+1,1)=ysplit(i,2)+1;
    end
    
    for i=1:1:spec
        clause1=mod(i,2);
        for j=1:1:spec
            clause2=mod(j,2);
            
            if clause1+clause2==1
                for xth=xsplit(i,1):1:xsplit(i,2)
                    for yth=ysplit(j,1):1:ysplit(j,2)
                        g_z_plus(xth,yth)=hs(1);
                        h_z_plus(xth,yth)=gs(1);
                    end
                end
                
                
            else
                for xth=xsplit(i,1):1:xsplit(i,2)
                    for yth=ysplit(j,1):1:ysplit(j,2)
                        g_z_plus(xth,yth)=hs(2);
                        h_z_plus(xth,yth)=gs(2);
                    end
                end
            end
            
        end
    end
    
elseif type==4
    if length(hs)~=1 || length(gs)~=1
        error('Specifications for type 4 is incorrect')
    end
    
    g_x_minus=-gs*ones(ny,nz);
    g_y_minus=-gs*ones(nx,nz);
    g_z_minus=-gs*ones(nx,ny);
    g_x_plus=gs*ones(ny,nz);
    g_y_plus=gs*ones(nx,nz);
    g_z_plus=gs*ones(nx,ny);

    h_x_minus=-hs*ones(ny,nz);
    h_y_minus=-hs*ones(nx,nz);
    h_z_minus=-hs*ones(nx,ny);
    h_x_plus=hs*ones(ny,nz);
    h_y_plus=hs*ones(nx,nz);
    h_z_plus=hs*ones(nx,ny);
    
elseif type==5
    
    if length(spec)~=1 || length(hs)~=1 || length(gs)~=1
        error('Specifications for type 5 is incorrect')
    end
    
    x_half=xlen/2;
    y_half=ylen/2;
    
    if spec>x_half || spec>y_half
        error('spec is the radius of circle. It must not exceed the domain of surface')
    end
    
    
    g_x_minus=zeros(ny,nz);
    g_y_minus=zeros(nx,nz);
    g_z_minus=zeros(nx,ny);
    g_x_plus=zeros(ny,nz);
    g_y_plus=zeros(nx,nz);
    g_z_plus=zeros(nx,ny);

    h_x_minus=zeros(ny,nz);
    h_y_minus=zeros(nx,nz);
    h_z_minus=zeros(nx,ny);
    h_x_plus=zeros(ny,nz);
    h_y_plus=zeros(nx,nz);
    h_z_plus=zeros(nx,ny);
    
    dx=xlen/(nx-1);
    dy=ylen/(ny-1);
    
    
    for xth=1:1:nx
        for yth=1:1:ny
            xPos=dx*(xth-1);
            yPos=dy*(yth-1);
            
            if sqrt((xPos-x_half)^2+(yPos-y_half)^2)<=spec
                g_z_plus(xth,yth)=gs;
                h_z_plus(xth,yth)=hs;
            end
        end
    end
    
    
else
    error('This type of wetting pattern is not available')
end


end