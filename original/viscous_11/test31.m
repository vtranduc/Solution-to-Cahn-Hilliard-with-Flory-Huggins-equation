function test31

clear
clc

N = 100;
% sphere grid points
[X,Y,Z] = sphere(N);
% the amplitude error for every grid point
AmplError = rand(N+1)*0.2 + 1; % test data
% graphical rendering
surf(AmplError.*X,AmplError.*Y,AmplError.*Z);



end