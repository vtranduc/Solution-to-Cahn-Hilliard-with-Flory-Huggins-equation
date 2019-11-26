function differentiate_septic_Hermite_spline

clear
clc

syms t

H00=1-35*t^4+84*t^5-70*t^6+20*t^7
H10=35*t^4-84*t^5+70*t^6-20*t^7

H01=t-20*t^4+45*t^5-36*t^6+10*t^7
H11=-15*t^4+39*t^5-34*t^6+10*t^7

H02=0.5*t^2-5*t^4+10*t^5-7.5*t^6+2*t^7
H12=2.5*t^4-7*t^5+6.5*t^6-2*t^7

H03=(1/6)*t^3-(2/3)*t^4+t^5-(2/3)*t^6+(1/6)*t^7
H13=(-1/6)*t^4+0.5*t^5-0.5*t^6+(1/6)*t^7

%First derivative

dH00_dt=diff(H00,t)
dH10_dt=diff(H10,t)
dH01_dt=diff(H01,t)
dH11_dt=diff(H11,t)
dH02_dt=diff(H02,t)
dH12_dt=diff(H12,t)
dH03_dt=diff(H03,t)
dH13_dt=diff(H13,t)

%Second derivative

d2H00_dt2=diff(diff(H00,t))
d2H10_dt2=diff(diff(H10,t))
d2H01_dt2=diff(diff(H01,t))
d2H11_dt2=diff(diff(H11,t))
d2H02_dt2=diff(diff(H02,t))
d2H12_dt2=diff(diff(H12,t))
d2H03_dt2=diff(diff(H03,t))
d2H13_dt2=diff(diff(H13,t))

%Third derivative

d3H00_dt3=diff(diff(diff(H00,t)))
d3H10_dt3=diff(diff(diff(H10,t)))
d3H01_dt3=diff(diff(diff(H01,t)))
d3H11_dt3=diff(diff(diff(H11,t)))
d3H02_dt3=diff(diff(diff(H02,t)))
d3H12_dt3=diff(diff(diff(H12,t)))
d3H03_dt3=diff(diff(diff(H03,t)))
d3H13_dt3=diff(diff(diff(H13,t)))

%Fourth derivative

d4H00_dt4=diff(diff(diff(diff(H00,t))))
d4H10_dt4=diff(diff(diff(diff(H10,t))))
d4H01_dt4=diff(diff(diff(diff(H01,t))))
d4H11_dt4=diff(diff(diff(diff(H11,t))))
d4H02_dt4=diff(diff(diff(diff(H02,t))))
d4H12_dt4=diff(diff(diff(diff(H12,t))))
d4H03_dt4=diff(diff(diff(diff(H03,t))))
d4H13_dt4=diff(diff(diff(diff(H13,t))))


end