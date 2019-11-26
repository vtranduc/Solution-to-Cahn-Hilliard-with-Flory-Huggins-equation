function test99

%Example plotting script for generating high resolution Matlab plots
%suitable for a thesis or paper.
close all
clear all

%Parameters for saved images
ImageDPI=500;
ImageSizeX=6;
ImageSizeY=4;
ImageFontSize=9;
FileLabel='WithFormatting';
FontName='Garamond';
AxisFontName='CMU Serif';

a=[0:0.1:100];
b=0.005*(a);
c=sin(a).*exp(-a/5);
d=exp(-a/100);

figure2 = figure(2);

% axes1 = axes('FontSize',ImageFontSize,'FontName',AxisFontName);
% xlim(axes1,[min(a) max(a)]);
% ylim(axes1,[-0.4 1]);
% box(axes1,'on');
% hold(axes1,'all');

h=plot(a,b,'k',a,c,'k',a,d,'--k');

set(h(1), 'LineWidth',2)
set(gca,'FontName',AxisFontName,'FontSize',ImageFontSize)
%Legend entries
%IMPORTANT NOTE
%Matlab has problems drawing the box around the legend if the font has been
%changed. The only way that I have found to get around this is to force
%blank spaces at the end of the longest legend entry. In this case 
%the 'Another Signal{ }' has three blank spaces after it.
m2=legend('Signal','Total Signal','Another Signal{ }',...
'Location','NorthWest');

xlabel('Frequency (Hz)','fontsize',ImageFontSize)
ylabel('Current Noise (A/\surdHz)','fontsize',ImageFontSize)

%====================
%Set the image size and save to FileLabel.png where FileLabel is set at line 9. 
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 ImageSizeX ImageSizeY])
print('-dpng', strcat(FileLabel, '.png') , strcat('-r',num2str(ImageDPI)))

end