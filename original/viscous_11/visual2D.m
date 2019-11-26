function frame_cap=visual2D(c_nodal,x_coord,y_coord,minz,maxz,t,...
    frequency_domain,structure_factor,ci_ave,limS,...
    logS)

range_=[min(min(c_nodal)),max(max(c_nodal))];

subplot(2,3,1)
surf(x_coord,y_coord,c_nodal);
% colorbar('eastoutside')
axis([0 1 0 1 minz maxz])
xlabel('x'),ylabel('y'),zlabel('c')
try
    caxis(range_)
catch
    %Forget about adjusting axis
end

subplot(2,3,2)
contourf(x_coord,y_coord,c_nodal);
colorbar('eastoutside')
xlabel('x');ylabel('y')
try
    caxis(range_)
catch
    %Forget about adjusting axis
end

subplot(2,3,3)
contourf(x_coord,y_coord,c_nodal,[0 ci_ave]);
xlabel('x');ylabel('y')

subplot(2,3,4)
plot(frequency_domain,structure_factor)
xlabel('Frequency domain, k')
ylabel('Structure factor, S')
ylim([0 limS])
grid on

subplot(2,3,5)
plot(frequency_domain,structure_factor)
xlabel('Frequency domain, k')
ylabel('Structure factor, S')
grid on

subplot(2,3,6)
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
ylabel('Natural log of structure factor, log(S)')

str=sprintf('Nonlinear Cahn Hilliard in 2D\n t=%d',t);
suptitle(str)
frame_cap=getframe(gcf);

end