function plot_septic_Hermite_spline

clear
clc

n=100;

x=linspace(0,1,n);

for order=0:1:4
    figure(order+1)
    for type=0:1:3
        subplot(2,2,type+1)
        plot(x,septic_Hermite_basis(x,0,type,order),x,septic_Hermite_basis(x,1,type,order))
        grid on
    end
end

end