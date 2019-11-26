function total_energy=test84(conc_abs,diffT,n1,n2,ne,ney,Two_chi_n1,x_fine,y_fine,nex3,ney3,w_set,dxdy)
energy=zeros(ney3,nex3);
inxo=-3;
inxoo=0;
total_energy=0;
for e=1:1:ne
    if mod(e,ney)==1
        inxo=inxo+3;
        inyo=-3;
        inxoo=inxoo+1;
    end
    elemental_energy=0;
    inyo=inyo+3;
    for iny=1:1:3
        for inx=1:1:3
            c2=1-conc_abs(inx,iny,e);
            energy(inyo+iny,inxo+inx)=diffT(inxoo,inx)...
                *(conc_abs(inx,iny,e)*log(conc_abs(inx,iny,e))/n1...
                +c2*log(c2)/n2...
                +conc_abs(inx,iny,e)*c2*Two_chi_n1(inxoo,inx)/2);
            elemental_energy=elemental_energy+w_set(inx,iny)*energy(inyo+iny,inxo+inx);
        end
    end
    total_energy=total_energy+elemental_energy;
end
total_energy=dxdy*total_energy;
contourf(x_fine,y_fine,energy,'edgecolor','none');
axis([0 1 0 1])
xlabel('x');ylabel('y')
colorbar

end

% function test84(conc_abs,diffT,n1,n2,ne,nex,ney)
% 
% nex3=nex*3;ney3=ney*3;
% energy=zeros(ney3,nex3);
% inxo=-3;
% for e=1:1:ne
%     if mod(e,ney)==1
%         inxo=inxo+3;
%         inyo=-3;
%     end
%     inyo=inyo+3;
%     for iny=1:1:3
%         for inx=1:1:3
%             energy(inyo+iny,inxo+inx)=conc_abs(inx,iny,e);
%         end
%     end
% end
% 
% surf(energy);
% % axis([0 1 0 1 0 1])
% 
% end