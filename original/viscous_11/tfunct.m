function [PHI,PHIX,PHIY,PHIXX,PHIYY,PHIXY]=tfunct(GP1,GP2,DX,DY)
PHI=fPHI(GP1,GP2,DX,DY);
PHIX=fPHIX(GP1,GP2,DX,DY);
PHIY=fPHIY(GP1,GP2,DX,DY);
PHIXX=fPHIXX(GP1,GP2,DX,DY);
PHIYY=fPHIYY(GP1,GP2,DX,DY);
PHIXY=fPHIXY(GP1,GP2,DX,DY);
end

function solution=F1(Z)
solution=1-3*Z^2+2*Z^3;
end
function solution=F2(Z,DZ)
solution=(Z-2*Z^2+Z^3)*DZ;
end
function solution=F3(Z)
solution=3*Z^2-2*Z^3;
end
function solution=F4(Z,DZ)
solution=(Z^3-Z^2)*DZ;
end


function solution=F1Z(Z,DZ)
solution=(-6*Z+6*Z^2)/DZ;
end
function solution=F2Z(Z)
solution=1-4*Z+3*Z^2;
end
function solution=F3Z(Z,DZ)
solution=(6*Z-6*Z^2)/DZ;
end
function solution=F4Z(Z)
solution=3*Z^2-2*Z;
end


function solution=F1ZZ(Z,DZ)
solution=(-6+12*Z)/DZ^2;
end
function solution=F2ZZ(Z,DZ)
solution=(-4+6*Z)/DZ;
end
function solution=F3ZZ(Z,DZ)
solution=(6-12*Z)/DZ^2;
end
function solution=F4ZZ(Z,DZ)
solution=(6*Z-2)/DZ;
end

function solution=fPHI(GP1,GP2,DX,DY)
solution=zeros(1,16);
solution(1)=F1(GP1)*F1(GP2);
solution(2)=F2(GP1,DX)*F1(GP2);
solution(3)=F1(GP1)*F2(GP2,DY);
solution(4)=F2(GP1,DX)*F2(GP2,DY);
solution(5)=F1(GP1)*F3(GP2);
solution(6)=F2(GP1,DX)*F3(GP2);
solution(7)=F1(GP1)*F4(GP2,DY);
solution(8)=F2(GP1,DX)*F4(GP2,DY);
solution(9)=F3(GP1)*F1(GP2);
solution(10)=F4(GP1,DX)*F1(GP2);
solution(11)=F3(GP1)*F2(GP2,DY);
solution(12)=F4(GP1,DX)*F2(GP2,DY);
solution(13)=F3(GP1)*F3(GP2);
solution(14)=F4(GP1,DX)*F3(GP2);
solution(15)=F3(GP1)*F4(GP2,DY);
solution(16)=F4(GP1,DX)*F4(GP2,DY);
end

function solution=fPHIX(GP1,GP2,DX,DY)
solution=zeros(1,16);
solution(1)=F1Z(GP1,DX)*F1(GP2);
solution(2)=F2Z(GP1)*F1(GP2);
solution(3)=F1Z(GP1,DX)*F2(GP2,DY);
solution(4)=F2Z(GP1)*F2(GP2,DY);
solution(5)=F1Z(GP1,DX)*F3(GP2);
solution(6)=F2Z(GP1)*F3(GP2);
solution(7)=F1Z(GP1,DX)*F4(GP2,DY);
solution(8)=F2Z(GP1)*F4(GP2,DY);
solution(9)=F3Z(GP1,DX)*F1(GP2);
solution(10)=F4Z(GP1)*F1(GP2);
solution(11)=F3Z(GP1,DX)*F2(GP2,DY);
solution(12)=F4Z(GP1)*F2(GP2,DY);
solution(13)=F3Z(GP1,DX)*F3(GP2);
solution(14)=F4Z(GP1)*F3(GP2);
solution(15)=F3Z(GP1,DX)*F4(GP2,DY);
solution(16)=F4Z(GP1)*F4(GP2,DY);
end

function solution=fPHIY(GP1,GP2,DX,DY)
solution=zeros(1,16);
solution(1)=F1(GP1)*F1Z(GP2,DY);
solution(2)=F2(GP1,DX)*F1Z(GP2,DY);
solution(3)=F1(GP1)*F2Z(GP2);
solution(4)=F2(GP1,DX)*F2Z(GP2);
solution(5) =F1(GP1)*F3Z(GP2,DY);
solution(6)=F2(GP1,DX)*F3Z(GP2,DY);
solution(7)=F1(GP1)*F4Z(GP2);
solution(8)=F2(GP1,DX)*F4Z(GP2);
solution(9)=F3(GP1)*F1Z(GP2,DY);
solution(10)=F4(GP1,DX)*F1Z(GP2,DY);
solution(11)=F3(GP1)*F2Z(GP2);
solution(12)=F4(GP1,DX)*F2Z(GP2);
solution(13)=F3(GP1)*F3Z(GP2,DY);
solution(14)=F4(GP1,DX)*F3Z(GP2,DY);
solution(15)=F3(GP1)*F4Z(GP2);
solution(16)=F4(GP1,DX)*F4Z(GP2);
end

function solution=fPHIXY(GP1,GP2,DX,DY)
solution=zeros(1,16);
solution(1)=F1Z(GP1,DX)*F1Z(GP2,DY);
solution(2)=F2Z(GP1)*F1Z(GP2,DY);
solution(3)=F1Z(GP1,DX)*F2Z(GP2);
solution(4)=F2Z(GP1)*F2Z(GP2);
solution(5)=F1Z(GP1,DX)*F3Z(GP2,DY);
solution(6)=F2Z(GP1)*F3Z(GP2,DY);
solution(7)=F1Z(GP1,DX)*F4Z(GP2);
solution(8)=F2Z(GP1)*F4Z(GP2);
solution(9)=F3Z(GP1,DX)*F1Z(GP2,DY);
solution(10)=F4Z(GP1)*F1Z(GP2,DY);
solution(11) = F3Z(GP1,DX)*F2Z(GP2);
solution(12)=F4Z(GP1)*F2Z(GP2);
solution(13)=F3Z(GP1,DX)*F3Z(GP2,DY);
solution(14)=F4Z(GP1)*F3Z(GP2,DY);
solution(15)=F3Z(GP1,DX)*F4Z(GP2);
solution(16)=F4Z(GP1)*F4Z(GP2);
end

function solution=fPHIXX(GP1,GP2,DX,DY)
solution=zeros(1,16);
solution(1)=F1ZZ(GP1,DX)*F1(GP2);
solution(2)=F2ZZ(GP1,DX)*F1(GP2);
solution(3)=F1ZZ(GP1,DX)*F2(GP2,DY);
solution(4)=F2ZZ(GP1,DX)*F2(GP2,DY);
solution(5)=F1ZZ(GP1,DX)*F3(GP2);
solution(6)=F2ZZ(GP1,DX)*F3(GP2);
solution(7)=F1ZZ(GP1,DX)*F4(GP2,DY);
solution(8)=F2ZZ(GP1,DX)*F4(GP2,DY);
solution(9)=F3ZZ(GP1,DX)*F1(GP2);
solution(10)=F4ZZ(GP1,DX)*F1(GP2);
solution(11)=F3ZZ(GP1,DX)*F2(GP2,DY);
solution(12)=F4ZZ(GP1,DX)*F2(GP2,DY);
solution(13)=F3ZZ(GP1,DX)*F3(GP2);
solution(14)=F4ZZ(GP1,DX)*F3(GP2);
solution(15)=F3ZZ(GP1,DX)*F4(GP2,DY);
solution(16)=F4ZZ(GP1,DX)*F4(GP2,DY);
end

function solution=fPHIYY(GP1,GP2,DX,DY)
solution=zeros(1,16);
solution(1)=F1(GP1)*F1ZZ(GP2,DY);
solution(2)=F2(GP1,DX)*F1ZZ(GP2,DY);
solution(3)=F1(GP1)*F2ZZ(GP2,DY);
solution(4)=F2(GP1,DX)*F2ZZ(GP2,DY);
solution(5)=F1(GP1)*F3ZZ(GP2,DY);
solution(6)=F2(GP1,DX)*F3ZZ(GP2,DY);
solution(7)=F1(GP1)*F4ZZ(GP2,DY);
solution(8)=F2(GP1,DX)*F4ZZ(GP2,DY);
solution(9)=F3(GP1)*F1ZZ(GP2,DY);
solution(10)=F4(GP1,DX)*F1ZZ(GP2,DY);
solution(11)=F3(GP1)*F2ZZ(GP2,DY);
solution(12)=F4(GP1,DX)*F2ZZ(GP2,DY);
solution(13)=F3(GP1)*F3ZZ(GP2,DY);
solution(14)=F4(GP1,DX)*F3ZZ(GP2,DY);
solution(15)=F3(GP1)*F4ZZ(GP2,DY);
solution(16)=F4(GP1,DX)*F4ZZ(GP2,DY);
end