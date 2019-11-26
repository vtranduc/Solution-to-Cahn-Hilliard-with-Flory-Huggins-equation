function solution=peak_adjust(peak,currentMax,lower_lim,higher_lim,...
    stretcher,compressor)

%Lower and higher limits must be betweeen 0 and 1

%Stretcher must be larger than 1

%Compressor must be between but not equal to 0 and 1

solution=currentMax;

current_pos=peak/solution;

if current_pos>higher_lim
    while current_pos>higher_lim
        solution=solution*stretcher;
        current_pos=peak/solution;
    end
elseif current_pos<lower_lim
    while current_pos<lower_lim
        solution=solution*compressor;
        current_pos=peak/solution;
    end
end

end