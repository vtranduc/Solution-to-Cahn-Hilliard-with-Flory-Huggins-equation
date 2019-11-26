function [nx_grad,ny_grad,x_grad,y_grad]=setUpFigure_2d(...
    figure_analysis,nx,ny,...
    nx_grad,ny_grad)
% nex3=nex*3;ney3=ney*3;
% gp=[(-sqrt(3/5)+1)/2 0.5 (sqrt(3/5)+1)/2];
% dx=1/nex;dy=1/ney;
% x_fine=zeros(1,nex3);y_fine=zeros(1,ney3);
% gpx=gp*dx;
% index=0;
% for iex=1:1:nex
%     z1=(iex-1)*dx;
%     for i=1:1:3
%         index=index+1;
%         x_fine(index)=z1+gpx(i);
%     end
% end
% gpy=gp*dy;
% index=0;
% for iey=1:1:ney
%     z1=(iey-1)*dy;
%     for i=1:1:3
%         index=index+1;
%         y_fine(index)=z1+gpy(i);
%     end
% end

if figure_analysis==1
%     if length(T)==1
%         diffTinx=diff*ones(1,nx)*T;
%         chi_n1_inx=ones(1,nx)*get_chi(T,entropy,T_theta)/n1;
%     elseif length(T)==2
%         T_=linspace(T(1),T(2),nx);
%         diffTinx=diff*T_;
%         chi_n1_inx=get_chi(T_,entropy,T_theta)/n1;
%     end
    if nx<nx_grad
        nx_grad=nx;
    end
    if ny<ny_grad
        ny_grad=ny;
    end
    x_grad=linspace(0,1,nx_grad);
    y_grad=linspace(0,1,ny_grad);
    nx_grad=round(x_grad*(nx-1))+ones(1,length(x_grad));
    ny_grad=round(y_grad*(ny-1))+ones(1,length(y_grad));
else
    nx_grad=NaN;ny_grad=NaN;x_grad=NaN;y_grad=NaN;
end

end