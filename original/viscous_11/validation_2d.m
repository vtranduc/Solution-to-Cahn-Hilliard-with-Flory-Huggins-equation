function validation_2d(diff,ci,T,n1,n2,time_dependent_temperature,...
    temperature_type,temperature_spec)


if time_dependent_temperature
    if temperature_type==1
        if length(temperature_spec)~=1 || temperature_spec<=0
            error('temperature_spec is the duration in which the temperature of T(2) varies')
        end
        if length(T)~=2
            error('You must specify two temperature ends for this type')
        end
    end
end


end