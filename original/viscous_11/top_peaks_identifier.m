function peaks=top_peaks_identifier(d_dx,d_dy,nx,ny)

%This assumes boundary condition has been properly applied

peaks=zeros(ny,nx);
% grads=zeros(ny,nx);
% for j=1:1:ny
%     for i=1:1:nx
%         grads(j,i)=sqrt(d_dx(j,i)^2+d_dy(j,i)^2);
%     end
% end

for i=2:1:ny-1
    for j=2:1:nx-1
        if d_dx(i,j-1)>0 && d_dx(i,j+1)<0 &&...
                d_dy(i-1,j)>0 && d_dy(i+1,j)<0 &&...
                d_dx(i-1,j-1)>0 && d_dy(i-1,j-1)>0 &&...
                d_dx(i-1,j+1)<0 && d_dy(i-1,j+1)>0 &&...
                d_dx(i+1,j-1)>0 && d_dy(i+1,j-1)<0 &&...
                d_dx(i+1,j+1)<0 && d_dy(i+1,j+1)<0 
%                 grads(i,j)<=grads(i-1,j) && grads(i,j)<=grads(i+1,j) &&...
%                 grads(i,j)<=grads(i,j-1) && grads(i,j)<=grads(i,j+1) &&...
%                 grads(i,j)<=grads(i-1,j-1) && grads(i,j)<=grads(i-1,j+1) &&...
%                 grads(i,j)<=grads(i+1,j-1) && grads(i,j)<=grads(i+1,j+1)
            peaks(i,j)=1;
        end
    end
end

%=========================
% FIXED FOR COORDINATES IN DIAGONAL DIRECTIONS!!!!!!

%IMMEDIATELY! NOW!
                
for i=2:1:ny-1
    if d_dy(i-1,1)>0 && d_dy(i+1,1)<0 &&...
            d_dx(i,2)<0
        peaks(i,1)=1;
    end
    if d_dy(i-1,nx)>0 && d_dy(i+1,nx)<0 &&...
            d_dx(i,nx-1)>0
        peaks(i,nx)=1;
    end
end

for i=2:1:nx-1
    if d_dx(1,i-1)>0 && d_dx(1,i+1)<0 &&...
            d_dy(2,i)<0
        peaks(1,i)=1;
    end
    if d_dx(ny,i-1)>0 && d_dx(ny,i+1)<0 &&...
            d_dy(ny-1,i)>0
        peaks(ny,i)=1;
    end
end

%=======================================

if d_dx(1,2)<0 && d_dy(2,1)<0 &&...
        d_dx(2,2)<0 && d_dy(2,2)<0
    peaks(1,1)=1;
end

if d_dx(1,nx-1)>0 && d_dy(2,nx)<0 &&...
        d_dx(2,nx-1)>0 && d_dy(2,nx-1)<0
    peaks(1,nx)=1;
end

if d_dx(ny,2)<0 && d_dy(ny-1,1)>0 &&...
        d_dx(ny-1,2)<0 && d_dy(ny-1,2)>0
    peaks(ny,1)=1;
end

if d_dx(ny,nx-1)>0 && d_dy(ny-1,nx)>0 &&...
        d_dx(ny-1,nx-1)>0 && d_dy(ny-1,nx-1)>0
    peaks(ny,nx)=1;
end

end