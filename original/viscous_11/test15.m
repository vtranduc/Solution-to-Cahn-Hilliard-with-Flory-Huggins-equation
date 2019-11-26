function test15(T)

chi=get_chi(T,1,1); %===============TO BE FIXED===========================
sol=vpasolve(


end

function solution=freeE(c1,n1,n2,chi)
c2=1-c1;
solution=c1*log(c1)/n1+c2*log(c2)/n2+c1*c2*chi;
end

function solution=dfreeE_dc1(c1,n1,n2,chi)
solution=(1/n1-1/n2)+log(c1)/n1-log(1-c1)/n2+(1-2*c1)*chi;
end