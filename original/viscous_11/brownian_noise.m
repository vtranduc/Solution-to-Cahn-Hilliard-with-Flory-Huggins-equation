function solution=brownian_noise(c,nx,ny,ci_fluc)

% clear
% clc
% 
% nx=100;
% ny=131;
% ci_fluc=0.02;
% c=0.7*ones(1,nx*ny*4);

solution=c;

noise=wgn(nx,ny,0);
max_noise=max(max(noise));
min_noise=min(min(noise));
noise=ci_fluc/(max_noise-min_noise)*noise;
ave_noise=mean(mean(noise));
noise=noise-ave_noise*ones(nx,ny);

index=-3;
for i=1:1:nx
    for j=1:1:ny
        index=index+4;
        solution(index)=solution(index)+noise(i,j);
    end
end

end