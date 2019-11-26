function test44

clear
clc

format long

a=rand(1,1)

dlmwrite('you_can_delete_this',a,'precision',16)

b=dlmread('you_can_delete_this');

if a-b==0
    
    error('good error')
end
a-b

end