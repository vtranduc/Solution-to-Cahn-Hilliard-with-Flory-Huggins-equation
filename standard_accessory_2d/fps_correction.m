function corrected_fps=fps_correction(video_method,fps,duration,iframe)
if video_method==1
    corrected_fps=fps;
elseif video_method==2
    corrected_fps=iframe/duration;
end
end