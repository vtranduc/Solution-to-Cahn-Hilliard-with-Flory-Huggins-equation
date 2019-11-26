function test2

clear
clc

n=100;

x=linspace(0,1,n);
y=linspace(0,1,n);
z=linspace(0,1,n);

% V=zeros(4,3,5);
% 
% % for i=1:n
% %     for j=1:n
% %         for k=1:n
% %             V(k,j,i)=basis(x(i),0,0)*basis(y(j),0,0)*basis(z(k),0,0);
% %         end
% %     end
% % end
% kk=linspace(0,1,5);
% for i=1:4
%     for j=1:3
%         for k=1:5
%             V(i,j,k)=kk(k);
%         end
%     end
% end
% 
% size(V)
% 
% V
% 
% isosurf(x,y,z,V,3,0.3,[0.5],{'cyan'},'black')
% xlabel('x'),ylabel('y'),zlabel('z')
% 
% hold on
% 
% scatter3(0.2,0.2,0.2,0.5)
figure(1)
plot(x,basis(x,0,1,1),x,basis(x,1,1,1))
grid on

figure(2)
plot(x,basis(x,0,0,1),x,basis(x,1,0,1))
grid on

figure(3)
plot(x,basis(x,0,0,2),x,basis(x,1,0,2))
grid on

figure(4)
plot(x,basis(x,0,1,2),x,basis(x,1,1,2))
grid on

end