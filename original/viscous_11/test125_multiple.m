function test125_multiple

clear
clc



% test125(diff,ci,T,n2,thermophoresis,distribution_type,T_distribution_spec)


test125(300,0.3,0.3,1,NaN,NaN,NaN)
test125(7300,0.3,0.3,1,NaN,NaN,NaN)
test125(65400,0.3,0.3,1,NaN,NaN,NaN)
test125(diff,0.3,0.3,1,NaN,NaN,NaN)



test125(200000,0.3,0.3,1,1,NaN,NaN)
test125(200000,0.3,0.3,1,2,NaN,NaN)
test125(200000,0.3,0.3,1,2,NaN,NaN)


end