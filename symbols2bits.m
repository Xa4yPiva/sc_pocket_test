function [bits] = symbols2bits(symbols)

bits = zeros(1, 2*length(symbols));

for i = 0 : length(symbols)-1
    switch symbols(i+1)
        case -3
            bits(2*i+1) = 0;
            bits(2*i+2) = 1;
        case -1
            bits(2*i+1) = 0;
            bits(2*i+2) = 0;
        case 1
            bits(2*i+1) = 1;
            bits(2*i+2) = 0;
        case 3
            bits(2*i+1) = 1;
            bits(2*i+2) = 1;
    end
end

end

