function test100

clear
clc


alpha=0;

A=[0 0 0 1
    1 1 1 1
    0 0 1 0
    3 2 1 0];


b=[0;0;0;1];


A\b


% basis(vals,orientation,type,order)

vals=linspace(0,1,100);

plot(vals,basis(vals,0,0,0),vals,basis(vals,0,1,0),vals,basis(vals,1,0,0),vals,basis(vals,1,1,0))

grid on
grid minor




end