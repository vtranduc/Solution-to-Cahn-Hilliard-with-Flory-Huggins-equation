function sol=convert_to_unit_vector(v)
sol=v/sqrt(sum(v.^2));
end