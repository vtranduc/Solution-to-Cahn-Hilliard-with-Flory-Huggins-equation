function index=get_index_3D(node,type)

%1=absolute val
% 2=d/dx
% 3=d/dy
% 4=d/dz
% 5=d2/dxdy
% 6=d2/dxdz
% 7=d2/dydz
% 8=d3/dxdydz

index=8*node-8+type;

end