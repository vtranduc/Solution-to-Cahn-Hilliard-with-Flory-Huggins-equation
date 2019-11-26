function test10

clear
clc

% fig=figure(1);
% 
% set(fig,'Position',[0 0 1920 1080])
% 
% m=3;
% n=4;
% 
% for i=1:1:m*n
%     subplot(m,n,i)
% end

[X, Y] = meshgrid(-1:.1:1,-1:.1:1);
G1 = subs(g(1), [x y], {X,Y});
G2 = subs(g(2), [x y], {X,Y});
quiver(X, Y, G1, G2)



end