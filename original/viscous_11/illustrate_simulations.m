function illustrate_simulations(xspinodal,yspinodal,xbinodal,ybinodal,T_list,ci_list,nSimulations)

plot(xspinodal,yspinodal,xbinodal,ybinodal)
hold on

for iSimulation=1:1:nSimulations
    T=T_list(iSimulation,:);
    ci=ci_list(iSimulation,:);
    if T(2)~=0 || ci(2)~=0
        if T(2)==0
            T(2)=T(1);
        elseif ci(2)==0
            ci(2)=ci(1);
        end
        plot(ci,T,ci,T,'o')
    else
        plot(ci,T,'o')
    end
end
hold off

grid on
xlim([0 1])
xlabel('Concentration, c')
ylabel('Temperature, T')
lgnd=legend('Spinodal curve','Binodal curve','Simulation range','Location','southeast');
set(lgnd,'color','none')

end