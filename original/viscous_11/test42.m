function [X,Y,Z]=test42(X,Y,Z,x_rotation)
n=length(X);
xyz=zeros(3,1);
Rx=get_Rx(x_rotation);
for i=1:1:n
    xyz(1)=X(i);
    xyz(2)=Y(i);
    xyz(3)=Z(i);
    xyz_=Rx*xyz;
    X(i)=xyz_(1);
    Y(i)=xyz_(2);
    Z(i)=xyz_(3);
end

end

function sol=get_Rx(x_rotation)
sol=zeros(3,3);
sol(1,1)=1;
sol(2,2)=cos(x_rotation);
sol(3,3)=sol(2,2);
sol(3,2)=sin(x_rotation);
sol(2,3)=-sol(3,2);
end