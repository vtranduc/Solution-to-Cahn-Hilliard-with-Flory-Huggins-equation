function test110

clear
clc

% f1=figure(1);clf;
% s1=subplot(1,2,1);
% surf(peaks(20));
% 
% s2=subplot(1,2,2);
% surf(peaks(20));
% s2Pos = get(s2,'position');
% 
% hb = colorbar('location','eastoutside');
% set(s2,'position',s2Pos);



% x = -5:0.1:5;
% [X,Y] = meshgrid(x);
% Z = X.^2 - Y.^2;
% figure(1)
% surf(X, Y, Z)
% grid on
% colorbar('EastOutside')
% figure_handle = gcf;
% % cbar_handle = findobj(figure_handle,'tag','Colorbar')
% % set(cbar_handle, 'YAxisLocation','right')

 t = 0:0.1:10 ; 
 hAxis(1) = subplot( 2, 1, 1 ) ; 
 plot( t, sin(t) ) ;
 hAxis(2) = subplot( 2, 1, 2 ) ;
 plot( t, cos(t) ) ;
 
 
pos(2) = 0.5 ;                         % Shift down.
pos(4) = 0.45 ;                        % Increase height.
set( hAxis(1), 'Position', pos ) ;



end