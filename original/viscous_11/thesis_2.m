function thesis_2

clear
clc

n1=1;
% n2=1;

entropy=1;
T_theta=1;


nT=1000;
tol=1.0e-6;
nInitialGuessPts=100;
guess_accuracy=10^-3;
suppress_warning=1;

chi_a=-0.5;
chi_b=1;

axis_font='CMU Serif';
font_size=9;
ImageSizeX_half_page=6;
ImageSizeY_half_page=4;
ImageSizeX_profile_6_figs=6;
ImageSizeY_profile_6_figs=8;
dpi=150;

overlapping_colors={'k','b','r','m'};

exporting_folder=[pwd '/background/'];

if ~exist(exporting_folder,'dir')
    mkdir(exporting_folder)
end

exporting_folder

%---------------------------------------------
n2_list=[1 1.3 3 5];

multi_figure=2;
% multi_figure==0 prints only 1 phase diagram with n2=n2_list(1)
% multi_figure==1 prints 4 phase diagrams in n2_list as subplots
% multi_figure==2 prints 4 phase diagrams in n2_list overlapping
%--------------------------------------------

index=0;

my_fig=figure(1);

if multi_figure==1

    for n2=n2_list
        index=index+1;

        subplot(2,2,index)

        [xspinodal,yspinodal,xbinodal,ybinodal]=...
            spinodal_binodal_chi_type_1(...
            n1,n2,entropy,T_theta,nT,tol,1,nInitialGuessPts,...
            chi_a,chi_b,...
            guess_accuracy);

        plot(xspinodal,yspinodal,'k.',xbinodal,ybinodal,'kx')
        grid on
        grid minor

        legend('Spinodal line','Binodal line','Location','southeast');


        xlabel('c');ylabel('T')
        title(sprintf('n_2 = %i',n2),'FontWeight','Normal')

    end
    
    set(gca,'FontName',axis_font,'FontSize',font_size)
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 ImageSizeX_half_page ImageSizeY_half_page])
    print('-dpng', strcat(exporting_folder,'multiple_phase_diagrams.png') , strcat('-r',num2str(dpi)))
    close(my_fig)
    
elseif multi_figure==0
    
    n2=n2_list(1);
    
    [xspinodal,yspinodal,xbinodal,ybinodal]=...
        spinodal_binodal_chi_type_1(...
        n1,n2,entropy,T_theta,nT,tol,1,nInitialGuessPts,...
        chi_a,chi_b,...
        guess_accuracy);

    plot(xspinodal,yspinodal,'-k',xbinodal,ybinodal,'--k')
    grid on
    grid minor

    legend('Spinodal line','Binodal line','Location','southeast');


    xlabel('c');ylabel('T')
%     title(sprintf('n_2 = %i',n2),'FontWeight','Normal')

    set(gca,'FontName',axis_font,'FontSize',font_size)
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 ImageSizeX_half_page ImageSizeY_half_page])
    print('-dpng', strcat(exporting_folder,['phase_diagrams_n2=' num2str(n2) '.png']) , strcat('-r',num2str(dpi)))
    close(my_fig)
    
elseif multi_figure==2
    
    hold on
    
    for in2=1:1:4
        
        n2=n2_list(in2);
%         index=index+1;
% 
%         subplot(2,2,index)

        [xspinodal,yspinodal,xbinodal,ybinodal]=...
            spinodal_binodal_chi_type_1(...
            n1,n2,entropy,T_theta,nT,tol,1,nInitialGuessPts,...
            chi_a,chi_b,...
            guess_accuracy);
        
        style_spinodal=['-' char(overlapping_colors(in2))];
        style_binodal=['--' char(overlapping_colors(in2))];

        plot(xspinodal,yspinodal,style_spinodal,xbinodal,ybinodal,style_binodal)

    end
    
    hold off
    
    xlabel('c');ylabel('T')
    lgnd=legend(['n_2=' num2str(n2_list(1)) ',spinodal'],['n_2=' num2str(n2_list(1)) ',binodal'],...
        ['n_2=' num2str(n2_list(2)) ',spinodal'],['n_2=' num2str(n2_list(2)) ',binodal'],...
        ['n_2=' num2str(n2_list(3)) ',spinodal'],['n_2=' num2str(n2_list(3)) ',binodal'],...
        ['n_2=' num2str(n2_list(4)) ',spinodal'],['n_2=' num2str(n2_list(4)) ',binodal'],...
        'Location','northwest');
    
    
%     pos=get(my_fig,'Position')
%     
%     set(my_fig,'Position',pos.*[1 -0.2 1 1.3])
%     
% %     return
    
    
    set(lgnd,'color','none');
    
    grid on
    grid minor
    
%     set(gca,'MinorGridLineStyle','--')
    
    set(gca,'FontName',axis_font,'FontSize',font_size)
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 ImageSizeX_profile_6_figs ImageSizeY_profile_6_figs])
    print('-dpng', strcat(exporting_folder,'overlapping_phase_diagrams.png') , strcat('-r',num2str(dpi)))
    
%     close(my_fig)
    
end




% plot(trend(:,2),trend(:,1),'k')
% xlabel('Time, t');ylabel('Concentration, c')
% grid on;grid minor




end