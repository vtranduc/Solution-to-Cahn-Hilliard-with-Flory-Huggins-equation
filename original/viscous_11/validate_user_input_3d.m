function validation=validate_user_input_3d(input)
inputSize=size(input);
if inputSize(1)~=1 || inputSize(2)~=1 || ~isreal(input) || isnan(input) || isinf(nan)
    validation=0;
else
    validation=1;
end
end