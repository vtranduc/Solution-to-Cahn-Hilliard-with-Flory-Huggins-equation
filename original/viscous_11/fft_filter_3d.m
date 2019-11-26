function data_fft_new=fft_filter_3d(...
    data_fft,xfft,yfft,zfft,filter_type,filter_spec)


% if filter_type==1 Eliminate high frequency noise
% if filter_type==2 Eliminate low frequency noise
% if filter_type==3 Get the Frequency within a range

data_fft_new=data_fft;
% if filter_type==1
%     if length(filter_spec)~=1
%         error('This is the radius of the circle, and should only have one element')
%     end
%     for i=1:1:length(xfft)
%         for j=1:1:length(yfft)
%             if sqrt(xfft(i)^2+yfft(j)^2)>=filter_spec
%                 data_fft_new(i,j)=0;
%             end
%         end
%     end
% elseif filter_type==2
%     if length(filter_spec)~=1
%         error('This is the radius of the circle, and should only have one element')
%     end
%     for i=1:1:length(xfft)
%         for j=1:1:length(yfft)
%             if sqrt(xfft(i)^2+yfft(j)^2)<=filter_spec
%                 data_fft_new(i,j)=0;
%             end
%         end
%     end

if filter_type==3
    if length(filter_spec)~=2 || filter_spec(1)>=filter_spec(2)
        error('This is a hollow circle with filter_spec=[smaller_radius larger_radius]')
    end
    for i=1:1:length(xfft)
        for j=1:1:length(yfft)
            for z=1:1:length(zfft)
                r=sqrt(xfft(i)^2+yfft(j)^2+zfft(z)^2);
                if r<filter_spec(1) || r>filter_spec(2)
                    data_fft_new(i,j,z)=0;
                end
            end
        end
    end
    
end


end