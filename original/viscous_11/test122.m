function test122

clear
clc

 freqbin = 72/12;
freqbins = [freqbin 72-freqbin]+1;
tsfit = zeros(72,1);
tsfit(freqbins) = tsdft(freqbins);
tsfit = ifft(tsfit);
mu = mean(ts);
tsfit = mu+tsfit;


figure(1)
subplot(2,1,1)
plot(abs(tsdft))
hold on
plot([freqbins(1) freqbins(1)], ylim, '-r')
plot([freqbins(2) freqbins(2)], ylim, '-r')
hold off
xlabel('Frequency Bin')
ylabel('Amplitude')
title('Original Fourier Transform')
q = [-36:35];
subplot(2,1,2)
plot([-36:35], fftshift(abs(tsdft)))
hold on
plot(1-[freqbins(1) freqbins(1)], ylim, '-r')
plot(73-[freqbins(2) freqbins(2)], ylim, '-r')
hold off
xlabel('Frequency Bin')
ylabel('Amplitude')
title('Shifted Fourier Transform Showing ‘Negative Frequencies’')


end