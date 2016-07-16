% Data Generation
N = 256;   % IQ-Complex Data
uncoded_data = zeros(1,N);

for i=1:N
    real_sign = rand;
    if (real_sign <= 0.5) 
        real_sign = 1;
    else
        real_sign = -1;
    end
    img_sign  = rand;
    if (img_sign <= 0.5) 
        img_sign = 1;
    else
        img_sign = -1;
    end
    uncoded_data(i) =  real_sign * (sqrt(1/2)) + img_sign * ( (sqrt(-1) * sqrt(1/2)) );
end

% Encoding Process:
[ coded_data, chaotic_sequence ] = chaos_coding(uncoded_data, 1);

% Decoding Process:
[ decoded_data ] = chaos_decoding(coded_data, chaotic_sequence);

% Calculate Errors:
n_errors = 0;
uncoded_data = round(uncoded_data,8);
decoded_data = round(decoded_data,8);
for i=1:N
    if( (real(uncoded_data(i)) - real(decoded_data(i))) ~= 0 && ...
        (imag(uncoded_data(i)) - imag(decoded_data(i))) ~= 0    )
        n_errors = n_errors + 1;
    end
end

if (n_errors > 0 )
    display('Errors were found...')
end