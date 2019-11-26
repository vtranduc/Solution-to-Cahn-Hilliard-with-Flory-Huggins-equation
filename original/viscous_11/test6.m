function test6

clear
clc

n=100;

x=linspace(0,1,n);


% hold on
% 
% for type=0:1:3
%     for orientation=[0 1]
%         plot(x,septic_Hermite_basis(x,orientation,type))
%     end
% end
% 
% hold off
% 
% plot(x,septic_Hermite_basis(x,0,0),x,septic_Hermite_basis(x,0,2))
% 
% grid on

% septic_Hermite_basis(x,0,1,0)
for order=0:1:3
    figure(order+1)
    for type=0:1:3
        subplot(2,2,type+1)
        plot(x,septic_Hermite_basis(x,0,type,order),x,septic_Hermite_basis(x,1,type,order))
        grid on
    end
end


end