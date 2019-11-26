function peaks=peaks_identifier(d_dx,d_dy,nx,ny)

% tol=10e-6;
format long

peaks=zeros(ny,nx);

grads=zeros(ny,nx);

size(d_dx)

ny
nx

for j=1:1:ny
    for i=1:1:nx
        grads(j,i)=sqrt(d_dx(j,i)^2+d_dy(j,i)^2);
        display('gradientssss')
        d_dx(j,i)
        d_dy(j,i)
%         if grads(j,i)==0
%             d_dx(j)
%             d_dy(i)
%             error('iust raise')
%         end
    end
end

d_dx

grads

for j=2:1:ny-1
    for i=2:1:nx-1
        if grads(j,i)<=grads(j-1,i) && grads(j,i)<=grads(j+1,i) &&...
                grads(j,i)<=grads(j,i-1) && grads(j,i)<=grads(j,i+1) &&...
                grads(j,i)<=grads(j-1,i-1) && grads(j,i)<=grads(j-1,i+1) &&...
                grads(j,i)<=grads(j+1,i-1) && grads(j,i)<=grads(j+1,i+1)
            grads(j,i)
            grads(j,i-1)
            grads(j,i+1)
%             error('iust raise')
            peaks(j,i)=1;
        end
    end
end

end