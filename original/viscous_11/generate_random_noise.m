function c_new=generate_random_noise(c,ci_ave,spec,nx,ny,xlen,ylen)

c_new=c;

c_nodal=extract_nodal_weights_2D(c,nx,ny);
c_relative=c_nodal-ci_ave*ones(ny,nx);

[c_fft,xfft,yfft,~]=fftn_new(c_relative,[ylen xlen],[ny,nx]);

i_complex=sqrt(-1);

for iy=1:1:length(yfft)
    for ix=1:1:length(xfft)
        
        if abs(real(c_fft(ix,iy)))<spec(1)

%         if 1==1
%             error('fdafd')
            c_fft(ix,iy)=c_fft(ix,iy)-real(c_fft(ix,iy))+spec(3)*randn(1,1)+spec(2);
%             c_fft(ix,iy)=c_fft(ix,iy)+10000;
            
        end
        if abs(imag(c_fft(ix,iy)))<spec(4)
%         if 1==1
            c_fft(ix,iy)=c_fft(ix,iy)+i_complex*(-imag(c_fft(ix,iy))+spec(5)*randn(1,1)+spec(6));
        end
    end
end

% % c_fft=test;
% 
c_nodal_new=real(ifftn(ifftshift(c_fft),[ny,nx]));

c_nodal_new=c_nodal_new+ci_ave*ones(ny,nx);

% 
% size(c_nodal_new)
% 
% size(c_relative)
% 
% sum(sum(abs(c_nodal_new-c_relative)))
% 
% error('systematic love')
% % 
% % c_nodal_new=interp2(linspace(0,xlen,length(xfft)),linspace(0,ylen,length(yfft)),c_nodal_new,0.1,0.2);
% 
% % c_nodal_new=(ifftn(ifftshift(c_fft),[ny,nx]));
% 
% % sum(sum(c_nodal_new-c_nodal))
% warning('0fdaf')

index=-3;
for ix=1:1:nx
    for iy=1:1:ny
        index=index+4;
        c_new(index)=c_nodal_new(iy,ix);
    end
end

end