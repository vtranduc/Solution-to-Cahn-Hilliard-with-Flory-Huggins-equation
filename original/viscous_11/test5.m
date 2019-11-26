function test5

clear
clc

c=dlmread('fourier.txt');

nx=101;
ny=101;

n=nx*ny;

result=zeros(1,n);

length(c)

index=-3;
for i=1:1:n
    index=index+4;
    result(i)=c(index);
end

% shaper=reshape(result,nx,ny);

% surf(shaper)

myFourier=fft(c,n);

abs(myFourier(1))

% 
% plot(abs(myFourier))

k=zeros(1,n);
speed=3*10^8;
for i=1:1:n
    k(i)=speed/i;
end
for i=1:1:n
    k(i)=2*pi/k(i);
end

plot(k,abs(myFourier))





end