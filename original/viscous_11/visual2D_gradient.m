function frame_cap=visual2D_gradient(c,nx,ny,x_coord,y_coord,minz,maxz,t,...
    frequency_domain,structure_factor,ci_ave,limS,logS,...
    structure_factor_1d,freq_domain,logS1d,limS1d,...
    xspinodal,yspinodal,xbinodal,ybinodal,T,ci)

c_nodal=extract_nodal_weights_2D(c,nx,ny);

range_=[min(min(c_nodal)),max(max(c_nodal))];

plot_row=2;
plot_col=4;
nfig=0;

%Standard
nfig=nfig+1;
subplot(plot_row,plot_col,nfig)
% subplot(8,6,[1,20])
%===============================
% surf(x_coord,y_coord,c_nodal,'linestyle','none')
surf(x_coord,y_coord,c_nodal)
%===============================
% camlight(100,0)
% light('Position',[0 0 2],'Style','local')
% lighting gouraud
colorbar('southoutside')
axis([0 1 0 1 minz maxz])
xlabel('x'),ylabel('y'),zlabel('c')
try
    caxis(range_)
catch
    %Forget about adjusting axis
end


%Confour, filled
nfig=nfig+1;
subplot(plot_row,plot_col,nfig)
% subplot(8,6,[25,44])
contourf(x_coord,y_coord,c_nodal,'linestyle','none');
colorbar('southoutside')
xlabel('x');ylabel('y')
try
    caxis(range_)
catch
    %Forget about adjusting axis
end

%Droplets, 2 colors
nfig=nfig+1;
subplot(plot_row,plot_col,nfig)
% subplot(8,6,[3,22])
contourf(x_coord,y_coord,c_nodal,[0 ci_ave],'linestyle', 'none');
xlabel('x');ylabel('y')

%Gradient

%THERE HAS TO BE BETTER WAY TO DEFINE THE PEAKS!!!

%FIX THIS LATER!!!

% nfig=nfig+1;
% subplot(plot_row,plot_col,nfig)
% % subplot(8,6,[27,46])
% [d_dx,d_dy]=extract_nodal_gradient_2D(c,nx,ny);
% quiver(x_coord,y_coord,d_dx,d_dy)
% hold on
% contour(x_coord,y_coord,c_nodal);
% peaks=top_peaks_identifier(d_dx,d_dy,nx,ny);
% for i=1:1:ny
%     for j=1:1:nx
%         if peaks(i,j)==1
%             plot(x_coord(j),y_coord(i),'ro')
%         end
%     end
% end
% hold off
% grid on
% xlabel('x')
% ylabel('y')


%Structure factor, fixed
nfig=nfig+1;
subplot(plot_row,plot_col,nfig)
% subplot(8,6,[5,24])
plot(frequency_domain,structure_factor)
xlabel('Frequency domain, k')
ylabel({'Structure factor,';'S'})
ylim([0 limS])
grid on

%Structure factor, flexible
% nfig=nfig+1;
% subplot(plot_row,plot_col,nfig)
% plot(frequency_domain,structure_factor)
% xlabel('Frequency domain, k')
% ylabel('Structure factor, S')
% grid on

%logarithmic structure factor, all
nfig=nfig+1;
subplot(plot_row,plot_col,nfig)
% subplot(8,6,[29,48])
plot(logS(:,1),logS(:,2),'*')
% hello=logS(:,2);
% ylim([0 max(hello)])
% try
%     xlim([0 max(logS(:,1))])
% catch
%     xlim([0 10^100])
% end
grid on
xlabel('Time, t')
ylabel({'Natural log of';'structure factor,';'log(S)'})

%Structure factor 1D
nfig=nfig+1;
subplot(plot_row,plot_col,nfig)
hold on
[num_fig,~]=size(structure_factor_1d);
for i=1:1:num_fig
    plot(freq_domain,structure_factor_1d(i,:))
end
hold off
grid on
ylim([0 limS1d])

% structure_factor_1d,freq_domain,logS1d,limS1d

%Natural log of struture factor

nfig=nfig+1;
subplot(plot_row,plot_col,nfig)
hold on
for i=1:1:num_fig
    plot(logS(:,1),logS1d(i,:),'*')
end
hold off
grid on

nfig=nfig+1;
subplot(plot_row,plot_col,nfig)
illustrate_simulation(xspinodal,yspinodal,xbinodal,ybinodal,T,ci)

str=sprintf('Nonlinear Cahn Hilliard in 2D\n t=%d',t);
suptitle(str)
frame_cap=getframe(gcf);

end