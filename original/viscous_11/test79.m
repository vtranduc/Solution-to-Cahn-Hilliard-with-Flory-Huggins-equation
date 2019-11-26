function test79

clear
clc

coxwell=["Ran","Kaori","Koo","Marika","Chie","Hiro","Kaori friend"];

ages=[38,18,29,24,21,40,21];

kojiko=7;

for jk=1:kojiko
    
    fprintf('\n%s is %i years old.\n',coxwell(jk),ages(jk))
    
    if ages(jk)>25
        fprintf('He or She is old.\n')
    else
        fprintf('He or She is young.\n')
    end
    
end


a=[1 2; 3 4; 5 6]

reshape(a,[1,6])

end