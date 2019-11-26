function create_movie_3d(visual,vid_title,fps)
video=VideoWriter(vid_title, 'Uncompressed AVI');
video.FrameRate=fps;
open(video)
writeVideo(video,visual);
close(video)
end