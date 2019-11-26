function test8

clear
clc

c=dlmread('fourier.txt');

[~,nfour]=size(c);

n=nfour/4;

sol=zeros(1,n);

for iii=1:1:n
    sol(iii)=c(iii*4-4+1);
end



% fft_=fft(sol,12)
% 
% fft_=abs(fft_)
% 
% fft_=fft_.^2

% plot(fft_)
% 
% Fs = 100;           % Sampling frequency
% t = -0.5:1/Fs:0.5;  % Time vector
% L = length(t);      % Signal length
% 
% X = 1/(4*sqrt(2*pi*0.01))*(exp(-t.^2/(2*0.01)));

% plot(t,X)

dev=zeros(1,n);
ave=mean(sol);

for iii=1:1:n
    dev(iii)=sol(iii)-ave;
end

fft_=fft(dev);

mag_=abs(fft_);
mag_sq=mag_.^2;

max(mag_sq);

length(mag_sq);

% plot(mag_sq(1:1000))

% dev=zeros(1,n)

q_len=100;

q_domain=linspace(1,12,q_len);

% k=0;
% 
% for q=q_domain
%     num=0.0;
%     for iii=1:1:n
%         num=num+exp(-1i*q*dev(iii));
%     end
%     num=abs(num)^2;
%     k=k+1;
%     S(k)=num/n;
% end
% 
% S=S./10^4

% plot(q_domain,S)
% axis([1 20 min(S) max(S)])
% 
% grid on

q_len=100;

q_domain=linspace(0,12,q_len);

S=zeros(1,q_len);

for iq=1:1:q_len
    
    cos_=0;
    sin_=0;
    
    for in=1:1:n
        cos_=cos_+(cos(q_domain(iq)*dev(in)));
        sin_=sin_+(sin(q_domain(iq)*dev(in)));
    end
    
    S(iq)=(abs(cos_)^2+abs(sin_)^2)/n;
    
end

% plot(q_domain,S)
% 
% S(1)

% x=linspace(0,2*pi,100);
% plot(x,sin(x))

grid on

fft_=fft(dev);

fft_=abs(fft_);

fft_=fft_.^2;

% max(fft_)

% plot(fft_)

sol_=reshape(sol,101,101);

% surf(sol_)

ave=mean(sol);

% dev'

sol2=fft2(sol_);

for iii=1:1:101
    for jjj=1:1:101
        sol2(iii,jjj)=abs(sol2(iii,jjj))^2;
    end
end
sol2(1,1)=0.7
surf( )

end