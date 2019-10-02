function solution=weight_adjuster(weights,dx,dy)

% clear
% clc

%=-----------
% weights=generateWeights();
% dx=0.25;
% dy=0.2;
% size(weights)
%=----------------
solution=weights;
for i=[1 3 5]
    for j=[2 4 6 8 10 12 14 16]
        solution(:,:,j,i)=solution(:,:,j,i)*dx;
    end
end
for i=[1 2 4]
    for j=[3 4 7 8 11 12 15 16]
        solution(:,:,j,i)=solution(:,:,j,i)*dy;
    end
end
for i=[1 3 5 7 9 11 13 15]
    solution(:,:,i,2)=solution(:,:,i,2)/dx;
end
for i=[1 2 5 6 9 10 13 14]
    solution(:,:,i,3)=solution(:,:,i,3)/dy;
end
for i=[1 3 5 7 9 11 13 15]
    solution(:,:,i,4)=solution(:,:,i,4)/(dx^2);
end
for i=[2 4 6 8 10 12 14 16]
    solution(:,:,i,4)=solution(:,:,i,4)/dx;
end

for i=[1 2 5 6 9 10 13 14]
    solution(:,:,i,5)=solution(:,:,i,5)/(dy^2);
end
for i=[3 4 7 8 11 12 15 16]
    solution(:,:,i,5)=solution(:,:,i,5)/dy;
end

end