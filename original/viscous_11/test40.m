function test40

clear
clc

gps=[(-sqrt(3/5)+1)/2 0.5 (sqrt(3/5)+1)/2];
w=[5/18 4/9 5/18];
gps_plus=zeros(3,3);
for i=1:1:3
    gps_plus(i,:)=gps*gps(i);
end

summation=0;
for i=1:1:3
	summation=summation+w(i)*basis(gps(i),0,1,0);
end
summation

summation=0;
for i=1:1:3
	summation=summation+w(i)*basis(gps(i),0,0,1);
end
summation

summation=0;
for i=1:1:3
	summation=summation+w(i)*basis(gps(i),0,1,1);
end
summation

summation=0;
for i=1:1:3
	summation=summation+w(i)*basis(gps(i),0,1,2);
end
summation

3^6*91

summation=0;
point=3;
for i=1:1:3
	summation=summation+w(i)*basis(gps_plus(point,i),0,1,1);
end
summation*gps(point)
basis(gps(point),0,1,0)
end