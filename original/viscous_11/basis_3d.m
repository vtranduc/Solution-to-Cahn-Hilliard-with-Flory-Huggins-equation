function basis_3d

clear
clc

n=10;

order=0;

x=linspace(0,1,n);
y=linspace(0,1,n);
z=linspace(0,1,n);

V=zeros(n,n,n);

isovals=[-0.5 -0.01 -0.005 -0.001 0.001 0.005 0.01 0.5];
face_colors={'yellow','magenta','cyan','red','green','blue','white','black'};
edge_color='none';

ifig=0;

for itype_z=[0 1]
    for iorientation_z=[0 1]
        ifig=ifig+1;
        figure(ifig)
        ipos=0;
        for itype_x=[0 1]
            for iorientation_x=[0 1]
                for itype_y=[0 1]
                    for iorientation_y=[0 1]
                        for i=1:n
                            for j=1:n
                                V(j,i,:)=...
                                    basis(x(i),iorientation_x,itype_x,order)*...
                                    basis(y(j),iorientation_y,itype_y,order)*...
                                    basis(z,iorientation_z,itype_z,order);
                            end
                        end
                        ipos=ipos+1;
                        subplot(4,4,ipos)
                        isosurf(x,y,z,V,3,0.7,isovals,face_colors,edge_color)
                        xlabel('x'),ylabel('y'),zlabel('z')
                        plotTitle=strcat(...
                            '\phi_',num2str(iorientation_x),'_',num2str(itype_x),'(\alpha)',...
                            '\phi_',num2str(iorientation_y),'_',num2str(itype_y),'(\beta)',...
                            '\phi_',num2str(iorientation_z),'_',num2str(itype_z),'(\gamma)');
                        title(plotTitle)
                        camlight 
                        lighting gouraud
                    end
                end
            end
        end
    end
end

end