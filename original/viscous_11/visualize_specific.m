function visualize_specific

clear
clc

c=dlmread([pwd '/Results/Run_23/concentration/iteration_1199']);

c_nodal=extract_nodal_weights_2D(c,31,31);

x_coord=linspace(0,1,31);
y_coord=linspace(0,1,31);

size(c_nodal)

surf(x_coord,y_coord,c_nodal)
colorbar('southoutside')
axis([0 1 0 1 0 1])
xlabel('x'),ylabel('y'),zlabel('c')

end