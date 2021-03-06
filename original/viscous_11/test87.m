function test87(ci,T,nex,ney,ci_ave,ci_fluc,include_fluc,diff,...
    counter_init,output_path,video_method,fps,duration,iframe,...
    ifig_raw,frame_pos,frame_size,minz,maxz,...
    figure_analysis,nx,ny,x_coord,y_coord,...
    frequency_domain,limS,...
    xspinodal,yspinodal,xbinodal,ybinodal,...
    ne,weights,adjusted_diffT,Two_chi_n1,n1,n2,...
    dxdy,diffTinx,chi_n1_inx,...
    alpha,beta,fc,limSstretcher,...
    nx_grad,ny_grad,x_grad,y_grad,contour_lgd_font_size)
fig=figure(ifig_raw);
set(fig,'Position',[frame_pos(1) frame_pos(2) frame_size(1) frame_size(2)])
output_conc=[output_path 'concentration/iteration_'];
output_time=[output_path 'time/iteration_'];
co=dlmread([output_conc '1']);
c=co;
time=0;
dt=inf;
visual=struct('cdata', cell(1,iframe),'colormap', cell(1,iframe));
[visual(1),E,ac_t,logS,limS]=visual_2d_2018(...
    figure_analysis,c,nx,ny,x_coord,y_coord,minz,maxz,time,...
    frequency_domain,ci_ave,limS,0,...
    xspinodal,yspinodal,xbinodal,ybinodal,T,ci,...
    ne,nex,ney,weights,adjusted_diffT,Two_chi_n1,n1,n2,...
    dxdy,iframe,0,diffTinx,chi_n1_inx,...
    co,dt,[0 0 0],alpha,beta,fc,limSstretcher,...
    nx_grad,ny_grad,x_grad,y_grad,contour_lgd_font_size);
for iVisual=2:1:iframe
    co=c;
    c=dlmread([output_conc num2str(iVisual)]);
    time=dlmread([output_time num2str(iVisual)]);
    [visual(iVisual),E,ac_t,logS,limS]=visual_2d_2018(...
        figure_analysis,c,nx,ny,x_coord,y_coord,minz,maxz,time,...
        frequency_domain,ci_ave,limS,logS,...
        xspinodal,yspinodal,xbinodal,ybinodal,T,ci,...
        ne,nex,ney,weights,adjusted_diffT,Two_chi_n1,n1,n2,...
        dxdy,iframe,E,diffTinx,chi_n1_inx,...
        co,dt,ac_t,alpha,beta,fc,limSstretcher,...
        nx_grad,ny_grad,x_grad,y_grad,contour_lgd_font_size);
end
close(fig)
if length(ci)==1 && length(T)==1
    vid_title=sprintf('2D Cahn Hilliard,%dx%d,ci=%1.3f,fluc=%1.4f,diff=%d,T=%0.3f,elapsed time=%1.3fs.avi',nex,ney,ci_ave,ci_fluc*include_fluc,diff,T,toc(counter_init));
elseif length(ci)==2 && length(T)==1
    vid_title=sprintf('2D Cahn Hilliard,%dx%d,ci=[%1.2f,%1.2f],fluc=%1.4f,diff=%d,T=%0.3f,elapsed time=%1.3fs.avi',nex,ney,ci(1),ci(2),ci_fluc*include_fluc,diff,T,toc(counter_init));
elseif length(ci)==1 && length(T)==2
    vid_title=sprintf('2D Cahn Hilliard,%dx%d,ci=%1.3f,fluc=%1.4f,diff=%d,T=[%1.2f,%1.2f],elapsed time=%1.3fs.avi',nex,ney,ci_ave,ci_fluc*include_fluc,diff,T(1),T(2),toc(counter_init));
elseif length(ci)==2 && length(T)==2
    vid_title=sprintf('2D Cahn Hilliard,%dx%d,ci=[%1.2f,%1.2f],fluc=%1.4f,diff=%d,T=[%1.2f,%1.2f],elapsed time=%1.3fs.avi',nex,ney,ci(1),ci(2),ci_fluc*include_fluc,diff,T(1),T(2),toc(counter_init));
end
vid_title=[output_path vid_title];
video=VideoWriter(vid_title, 'Uncompressed AVI');
fps_corrected=fps_correction(video_method,fps,duration,iframe);
video.FrameRate=fps_corrected;
open(video)
writeVideo(video,visual);
close(video)
fig=figure(ifig_raw);
set(fig,'Position',[frame_pos(1) frame_pos(2) frame_size(1) frame_size(2)])
movie(gcf,visual,1,fps_corrected)
close(fig)
end