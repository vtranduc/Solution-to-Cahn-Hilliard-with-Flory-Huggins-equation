function test9

clear
clc

nx=5;
ny=6;
ci_fluc = 0.01;

low_c=0.5;
high_c=0.8;

co=co_linear_gradient_x(low_c,high_c,ci_fluc,nx,ny)

% mean(co)

% A=rand(2,2)
% 
% B=reshape(A,[4,1])

abs_w=extract_nodal_weights_2D(co,nx,ny)

mean(abs_w)

end