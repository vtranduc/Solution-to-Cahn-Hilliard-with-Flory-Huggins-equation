function solution=sf_oneElement(chi,c_,co_,weights,ilocal,diffT,n1,n2,dx,dy,dt)

% weights=weight_adjuster(weights_,dx,dy);

% clear
% clc
% 
% %========inputs=============
% chi=0.9
% 
% c_=ones(1,16)*0.9;
% 
% weights=generateWeights();
% 
% %These will also be used in sj, thus should be computed in higher file
% % con=conc(c_,weights,1)
% % conx=conc(c_,weights,2)
% % cony=conc(c_,weights,3)
% % conxx=conc(c_,weights,4)
% % conyy=conc(c_,weights,5)
% 
% %These will be used throughout the main program, thus should be saved in
% %mainfile
w=[5/18 4/9 5/18];
% 
% % phi=conc(ones(1,16),weights,1)
% % phix=conc(ones(1,16),weights,2)
% % phiy=conc(ones(1,16),weights,3)
% 
% %This is a particular phi_i we are looking at. This depends on in which
% %direction the element lies relative to the node
% ilocal=12;
% 
% % order_type=1;
% 
% cont=0.5;
% diff=500;
% T=0.6;
% n1=1;
% n2=10;
% 
% dx=0.02;
% dy=0.02;
%===========================
% basis(0.112702,0,0,0)
% display('end the test')
% return
%===========================

con=conc(c_,weights,1);
cono=conc(co_,weights,1);
conx=conc(c_,weights,2);
cony=conc(c_,weights,3);
conxx=conc(c_,weights,4);
conyy=conc(c_,weights,5);

cont=(con-cono)/dt;

%No computation is required here, so it does not have to be input
phi=weights(:,:,ilocal,1);
% phix=weights(:,:,ilocal,2)
% phiy=weights(:,:,ilocal,3)
phixx=weights(:,:,ilocal,4);
phiyy=weights(:,:,ilocal,5);

solution=0;
for ix=1:1:3
    for iy=1:1:3
        solution=solution+w(ix)*w(iy)*(...
            phi(ix,iy)*cont(ix,iy)-diffT(ix)*phi(ix,iy)*...
            ((-1/(con(ix,iy)^2*n1)+1/((1-con(ix,iy))^2*n2))*...
            (conx(ix,iy)^2+cony(ix,iy)^2)+...
            (1/(con(ix,iy)*n1)+1/((1-con(ix,iy))*n2)-2*chi(ix))*...
            (conxx(ix,iy)+conyy(ix,iy)))+...
            (conxx(ix,iy)+conyy(ix,iy))*...
            (phixx(ix,iy)+phiyy(ix,iy))...
            );

    end
end
solution=dx*dy*solution;
end