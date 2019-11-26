function main_multiple_2018

clear
clc

main_2d_2018(10000,0.3,[0.35 0.36],1,0,5,[4 0.1])
main_2d_2018(10000,0.3,[0.35 0.36],1,0,5,[1 0.2])
main_2d_2018(10000,0.3,[0.35 0.36],1,-100000,5,[4 0.1])
main_2d_2018(10000,0.3,[0.35 0.36],1,100000,5,[4 0.1])
main_2d_2018(10000,0.3,[0.35 0.36],1,100000,5,[1 0.2])

main_2d_2018(10000,0.3,[0.35 0.36],1,-100000,5,[1 0.2])

% %Order: diff,ci,T,n2
% 
% % User inputs
% 
% simuls=[10000 0.4 0.36 1
%     10000 0.4 0.3 1];
% 
% %-------------------------------
% 
% [nSimuls,~]=size(simuls)
% 
% % return
% 
% for simul=1:1:nSimuls
%     main_2d_2018(simuls(simul,1),simuls(simul,2),simuls(simul,3),simuls(simul,4));
% end

end