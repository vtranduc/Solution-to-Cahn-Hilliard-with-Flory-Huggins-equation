function [data_fft,xfft,yfft,zfft]=fftn_new(data,lens,mesh)
data_size=size(data);
dimension=length(data_size);
if dimension==2 && (data_size(1)==1 || data_size(2)==1)
    dimension=1;
    nfft2=2^nextpow2(length(data));
    data_fft=fftshift(fft(data,nfft2));
else
    nfft2=(2*ones(1,dimension)).^nextpow2(data_size);
    data_fft=fftshift(fftn(data,nfft2));
end
fs=mesh./lens;
if dimension==2
    xfft=fs(1)*(-nfft2(1)/2:nfft2(1)/2-1)/nfft2(1);
    yfft=fs(2)*(-nfft2(2)/2:nfft2(2)/2-1)/nfft2(2);
    zfft=NaN;
elseif dimension==3
    
%     error('has not been coded')

    xfft=fs(1)*(-nfft2(1)/2:nfft2(1)/2-1)/nfft2(1);
    yfft=fs(2)*(-nfft2(2)/2:nfft2(2)/2-1)/nfft2(2);
    zfft=fs(3)*(-nfft2(3)/2:nfft2(3)/2-1)/nfft2(3);
    
    
elseif dimension==1
    xfft=fs*(-nfft2/2:nfft2/2-1)/nfft2;
    yfft=NaN;
    zfft=NaN;
end

end