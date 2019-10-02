function illustrate_simulation(xspinodal,yspinodal,xbinodal,ybinodal,T,ci)

if length(T)==2 || length(ci)==2
    if length(T)==1
        T_=[T T];
    else
        T_=T;
    end
    if length(ci)==1
        ci_=[ci ci];
    else
        ci_=ci;
    end
    plot(xspinodal,yspinodal,xbinodal,ybinodal,ci_,T_,ci_,T_,'o')
else
    plot(xspinodal,yspinodal,xbinodal,ybinodal,ci,T,'o')
end

grid on
grid minor
xlim([0 1])
xlabel('Concentration, c')
ylabel('Temperature, T')
lgnd=legend('Spinodal curve','Binodal curve','Simulation range','Location','southeast');
set(lgnd,'color','none')

end