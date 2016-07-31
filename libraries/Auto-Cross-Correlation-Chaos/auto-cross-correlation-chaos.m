% -- * -- * -- * -- * -- * -- * -- * -- * -- * -- * -- * -- * -- * -- *
% Rho Corporation @ 2016 - PhD Communications
% -- * -- * -- * -- * -- * -- * -- * -- * -- * -- * -- * -- * -- * -- *
% This is a test for the chaotic coding-decoding process.
% Parameters:
%       N   = Number of modulated IQ Complex data
%  
%       C_N = Defines the lenght of the chaotic sequence, for Real (I) and
%             Imaginary (Q) portions.
%
%       Mod_Map = Symbol Modulation Map Used for transmission

clear
clc

% Data Generation
N       = 256;   % IQ-Complex Data
C_N     = 256;
Mod_Map = 'QPSK';
uncoded_data = zeros(1,N);

% QPSK-Signal Generation 
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
[ coded_data, chaotic_seq_set ] = chaos_coding_1( uncoded_data, ...
                                                  C_N,          ...
                                                  Mod_Map,      ...
                                                  1);

% Decoding Process:
[ decoded_data ] = chaos_decoding_1( coded_data,      ...
                                     N,               ...
                                     chaotic_seq_set, ...
                                     Mod_Map                );

% Calculate Errors:
n_errors = 0;
uncoded_data = round(uncoded_data,8);
decoded_data = round(decoded_data,8);
for i=1:N
    if( uncoded_data(i) ~= decoded_data(i) )
        n_errors = n_errors + 1;
    end
end

if (n_errors > 0 )
    display('Errors were found...')
else
    display('Successful Coding! :-) No Errors Were Found...')
end