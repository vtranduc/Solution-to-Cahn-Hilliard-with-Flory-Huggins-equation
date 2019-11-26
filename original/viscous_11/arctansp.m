function theta=arctansp(x_val,y_val)
%y_val is y-val, x_val is x-val
if x_val<0
    theta=pi-atan(y_val/-x_val);
elseif x_val>0
    if y_val>=0
        theta=atan(y_val/x_val);
    elseif y_val<0
        theta=atan(y_val/x_val)+2*pi;
    end
elseif x_val==0
    if y_val>0
        theta=pi/2;
    elseif y_val<0
        theta=1.5*pi;
    elseif y_val==0
        theta=0;
    end
end
end