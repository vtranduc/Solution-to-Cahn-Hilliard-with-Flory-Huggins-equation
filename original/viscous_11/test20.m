function test20

clear
clc

a=nan;

a=[i 0];

isreal(a);

fminbnd(@(x)sin(x),0,pi/2)

linspace(1,0,10)

end