function thesis_8_growth_factor

clear
clc
close all

diff=100000;
c_ave=0.3;
n1=1;
n2=1;
chi_a=-0.5;
chi_b=1;
fc=25;
n=100;
f_domain=linspace(0,fc,n);

diffs=[100000 20000 100000 100000 100000];
Ts=[0.3,0.23,0.35,0.36,0.4];
n2s=[1 10 1 1 1];
c_aves=[0.3 0.76 0.3 0.3 0.3];

dpi=150;
axis_font='CMU Serif';
font_size=9;

ImageSizeX_profile_6_figs=5;
ImageSizeY_profile_6_figs=8;

export_folder=[pwd '/figures/'];



if ~exist(export_folder,'dir')
    mkdir(export_folder);
end

hold on
for iT=1:1:5
    range=get_growth(diffs(iT),c_aves(iT),Ts(iT),n1,n2s(iT),chi_a,chi_b,f_domain);
    plot(f_domain,range)
end
hold off
axis([-inf inf -0.2*10^9 inf])
grid on
grid minor

xlabel('$\sqrt{{{k_1}}^2 + {k_2}^2}$','Interpreter','latex')
ylabel('R(k_1,k_2)')
title('Growth factor for D=100,000, c_{ave}=0.3, n_2=1,')
legend('T=0.3 (Unstable)','T=0.32 (Unstable)','T=0.35 (Metastable)','T=0.36 (Metastable)','T=0.4 (Stable)','Location','northwest')


return

set(gca,'FontName',axis_font,'FontSize',font_size)
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 ImageSizeX_profile_6_figs ImageSizeY_profile_6_figs])
print('-dpng', strcat(export_folder,'fc_demo.png') , strcat('-r',num2str(dpi)))

end

function sol=get_growth(diff,c_ave,T,n1,n2,chi_a,chi_b,f_domain)
sol=zeros(1,length(f_domain));
chi=chi_a+chi_b/T;
d2f_dc2=1/(c_ave*n1)+1/((1-c_ave)*n2)-2*chi/n1;
chunk=diff*T*d2f_dc2;
domain=(2*pi*f_domain).^2;
for i=1:1:length(f_domain)
    sol(i)=chunk*domain(i)+domain(i)^2;
end
sol=-sol;
end