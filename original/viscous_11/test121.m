function test121

clear
clc

% Fs = 100;           % Sampling frequency
% t = -0.5:1/Fs:0.5;  % Time vector 
% L = length(t);      % Signal length
% 
% X = 1/(4*sqrt(2*pi*0.01))*(exp(-t.^2/(2*0.01)));
% 
% plot(t,X)
% title('Gaussian Pulse in Time Domain')
% xlabel('Time (t)')
% ylabel('X(t)')
% 
% n = 2^nextpow2(L);
% 
% n'
% 
% Y = fft(X,n);
% 
% f = Fs*(0:(n/2))/n;
% P = abs(Y/n);
% 
% plot(f,P(1:n/2+1)) 
% title('Gaussian Pulse in Frequency Domain')
% xlabel('Frequency (f)')
% ylabel('|P(f)|')

%=============================================================





Fs = 1000;            % Sampling frequency                    
T = 1/Fs;             % Sampling period       
L = 1500;             % Length of signal
t = (0:L-1)*T;        % Time vector

% size(t)
% 
% return

% f = Fs*(0:(L/2))/L;
% f'
% % size(f)
% % L/2
% % test=0:(L/2)*Fs/L;
% % test'
% return

S = 0.7*sin(2*pi*50*t) + sin(2*pi*120*t);

X = S + 2*randn(size(t));

plot(1000*t(1:50),X(1:50))
title('Signal Corrupted with Zero-Mean Random Noise')
xlabel('t (milliseconds)')
ylabel('X(t)')


Y = fft(X);

size(X)

P2 = abs(Y/L);
P1 = P2(1:L/2+1);


% size(P2)
% L/2+1
% 
% return

% take=1:L/2+1;
% take'
% L/2+1
% return


P1(2:end-1) = 2*P1(2:end-1);

f = Fs*(0:(L/2))/L;


f'

plot(f,P1) 
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')


end