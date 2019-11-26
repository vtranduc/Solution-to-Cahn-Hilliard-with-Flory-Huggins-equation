function test83

clear
clc

% A=[0 1 0 0 0
%     1 0 1 1 1
%     0 0 1 0 0
%     0 0 0 2 3
%     1 0 0 0 0];
% 
% b=zeros(5,1);
% b(1)=1;
% b(5)=1;
% 
% x=A\b

A=[1 0 0 0
    1 1 1 1
    0 1 0 0
    0 1 2 3];

b=zeros(4,1);

b(2)=1;

x=A\b


AA=zeros(10,8);

AA(1:4,1:4)=A;
AA(5:8,5:8)=A;

bb=zeros(10,1);

bb(1)=1;
bb(7)=1;

xx=AA\bb

AA*xx

%-------------

AA(9,3)=2;
AA(9,7)=2;
AA(10,3)=2;
AA(10,4)=6;
AA(10,7)=2;
AA(10,8)=6;

bb(9)=1;

xx=AA\bb

AA*xx


end