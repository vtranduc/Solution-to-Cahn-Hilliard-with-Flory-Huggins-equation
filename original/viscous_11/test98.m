function test98

clear
clc

x=0.42;
tester(x,1,10)


end

function sol=tester(x,n1,n2)
sol=(log(x)+1)/n1-(log(1-x)+1)/n2-0.5*(1-2*x);
end