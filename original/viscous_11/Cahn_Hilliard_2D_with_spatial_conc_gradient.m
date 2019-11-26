function Cahn_Hilliard_2D_with_spatial_conc_gradient

clear
clc

% ======================= Start the counter ===============================

tic
counter_init=cputime;

% ===================== User defined variables ============================

nex=25;
ney=25;

dto=1.0e-5;
dt=1.0e-5;

obs_t=2.0e-1; %Observation time

diff=5000;

T=0.65;

ci = 0.7; %Overall initial concentrations
ci=[0.5,0.6];
% ci=0.7
include_fluc=1;
ci_fluc = 0.01; %Initial fluctuation

nr_tol=1.0e-6;
nr_max_iteration=20;

time_out=3600*(1/6);

entropy=1;
T_theta=1;

fps=10;

dt_down=0.5;
dt_up=1.2;
dt_min=1.0e-20;
dt_ideal=5.0e-8;
bypass_tol=5;

max_frame=500;

xlen=1;ylen=1;

alpha=200;
beta=8;

frame_size=[1920 1080];
frame_pos=[0 0];

limS=1; %Guessed value. This will be adjusted as necessary
limShigher=1;
limSlower=0.01;
limSstretcher=10;
limScompressor=0.1;

fc=15; %Put fc=0 if you wanna use Nyquist frequency

num_str_fac_1d=6;

limS1d=1;
limS1dhigher=1;
limS1dlower=0.01;
limSstretcher1d=10;
limScompressor1d=0.1;

nT_phase=1000000000;
T_min_phase=0.001;

n1=1;
n2=10;

ifig_phase=9;

proceed_request=1; % Yes=1

% === Validate user inputs ================================================

%COMPLETE IT LATER IN ANOTHER FUNCTION!

%CHECK TO MAKE SURE CI+/-CI_FLUC IS WITHIN 0 AND 1!!! 

% === Illustrate phase diagram ============================================

fprintf('Generating phase diagram...\n')

[xspinodal,yspinodal,xbinodal,ybinodal]=...
    phase_diagram(n1,n2,entropy,T_theta,nT_phase,T_min_phase);

figure(ifig_phase)
illustrate_simulation(xspinodal,yspinodal,xbinodal,ybinodal,T,ci)

if proceed_request==1
    proceed=input('f,%1.2f],fluc=%1.4f,diff=%d,T=%0.3f,elapsed time=%1.3fs.avi',nex,ney,ci(1),ci(2),ci_fluc*include_fluc,diff,T,elapsed_);
end
video=VideoWriter(vid_title, 'Uncompressed AVI');
video.FrameRate=fps;
open(video)
writeVideo(video,visual);
close(video)

% === Identify transition time ============================================

figure(2)
plot(logS(:,1),logS(:,2),'x')
grid on

% === Conclusion ==========================================================
% dlmwrite('fourier.txt',c)
toc

fig=figure(3);
set(fig,'Position',[frame_pos(1) frame_pos(2) frame_size(1) frame_size(2)])
movie(gcf,visual,1)

end