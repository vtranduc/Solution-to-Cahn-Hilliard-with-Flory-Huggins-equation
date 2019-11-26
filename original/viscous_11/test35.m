function test35

clear
clc

w=[5/18 4/9 5/18];
gp=[(-sqrt(3/5)+1)/2 0.5 (sqrt(3/5)+1)/2];

% solution=basis(vals,orientation,type,order)

sol=0;
for i=1:1:3
    sol=sol+w(i)*basis(gp(i),0,0,0);
end

sol

end