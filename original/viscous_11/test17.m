function test17

clear
clc


a=8+4i
b=646i
c=425

d=NaN

A=[a b c d]

real(A)

any(imag(A))
any(isnan(A))

if -5646
    disp('yes')
else
    disp('no')
end

e=inf

imag(e)
isnan(e)
isinf(e)

A=[a b c d e]

any(isinf(A))

end