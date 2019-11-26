function co=co_linear_gradient_x(low_c,high_c,ci_fluc,nx,ny)

slope=(high_c-low_c)/(nx-1);

co_nodal=zeros(ny,nx);

for x_coord=1:1:nx-1
    c_specific=slope*(x_coord-1)+low_c;
    for y_coord=1:1:ny
        co_nodal(y_coord,x_coord)=c_specific+(ci_fluc*(2*rand(1,1)-1));
    end
    ci_diff=c_specific-mean(co_nodal(:,x_coord));
    for y_coord=1:1:ny
        co_nodal(y_coord,x_coord)=co_nodal(y_coord,x_coord)+ci_diff;
    end
end

for y_coord=1:1:ny
	co_nodal(y_coord,nx)=high_c+(ci_fluc*(2*rand(1,1)-1));
end
ci_diff=high_c-mean(co_nodal(:,nx));
for y_coord=1:1:ny
	co_nodal(y_coord,nx)=co_nodal(y_coord,nx)+ci_diff;
end

co_nodal=reshape(co_nodal,[nx*ny,1]);
co=4*nx*ny;

for i=1:1:nx*ny
    co((i-1)*4+1)=co_nodal(i);
end

end