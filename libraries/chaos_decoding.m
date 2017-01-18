% -- * -- * -- * -- * -- * -- * -- * -- * -- * -- * -- * -- * -- * -- *
% Rho Corporation @ 2016 - PhD Communications
% -- * -- * -- * -- * -- * -- * -- * -- * -- * -- * -- * -- * -- * -- * 

function  [ decoded_data ] = chaos_decoding(    coded_sequence,     ...
                                                Tx_Symbol_N, ...
                                                chaotic_seq_set,   ...
                                                Mod_Symbol_Map             )
%
% Description:
% 
%  chaos_decoding decodes complex data based on the received sequence, and
%  its correlation with the sent symbol, based on real and Img IQ-data.
%
%  Inputs:
%      coded_sequence  = complex IQ-data with lenght Tx_Symbol_N
%      Tx_Symbol_N     = Number of IQ-Transmitted Symbols
%      chaotic_seq_set = set of chaotic sequences used for encoding
%      Mod_Symbol_Map  = the symbol mapping it's essential for decoding
%    
%  Outputs:
%      decoded_data     = estimated complex IQ-data
%      

% IQ-Data Sequence Length:
N = Tx_Symbol_N;

switch Mod_Symbol_Map
    case 'QPSK'
        Symbol_Map = [ sqrt(0.5), (-1) * sqrt(0.5) ];
    otherwise
        % QPSK Modulation
        Symbol_Map = [ sqrt(0.5), (-1) * sqrt(0.5) ];
end

% Separate IQ Signals:
Chaos_Seq_Real = real(coded_sequence);
Chaos_Seq_Img  = imag(coded_sequence);

% Circular Matix:
circ_real = zeros(N,N);
circ_img  = zeros(N,N);

%Initialization:
circ_real(1,:) = chaotic_seq_set{1};
circ_img (1,:) = chaotic_seq_set{2};

for i=2:N
    circ_real(i,:) = [circ_real(i-1,2:end), circ_real(i-1,1)];
    circ_img (i,:) = [circ_img(i-1,2:end),  circ_img(i-1,1) ];
end

Chaos_Seq_Real = circ_real * Chaos_Seq_Real;
Chaos_Seq_Real = Chaos_Seq_Real'/N;
Chaos_Seq_Img  = circ_img  * Chaos_Seq_Img;
Chaos_Seq_Img  = Chaos_Seq_Img'/N;

% Data to be decoded
decoded_data   = zeros(1,N);

% Decoding Process:
Partial_Dencoding_Real = zeros(1,N);
Partial_Dencoding_Img  = zeros(1,N);

for i=1:Tx_Symbol_N
    switch Mod_Symbol_Map
        case 'QPSK'
            % Real Part (I)
            Chaos_Seq_Real(i) = Chaos_Seq_Real(i) - Symbol_Map(1);  % 1st Symbol
            [first_symbol_corr, lags_1, bounds_1]  = crosscorr(Chaos_Seq_Real,chaotic_seq_set{1}, length(chaotic_seq_set{1})-1);
            Chaos_Seq_Real(i) = Chaos_Seq_Real(i) - Symbol_Map(2);  % 2nd Symbol
            [second_symbol_corr, lags_2, bounds_2] = crosscorr(Chaos_Seq_Real,chaotic_seq_set{1}, length(chaotic_seq_set{1})-1);
            if (first_symbol_corr(length(chaotic_seq_set{1})) > second_symbol_corr(length(chaotic_seq_set{1})) )    
                Partial_Dencoding_Real(i) = Symbol_Map(1);
            else
                Partial_Dencoding_Real(i) = Symbol_Map(2);
            end
            % Imaginary Part (Q)
            Chaos_Seq_Img(i) =  Chaos_Seq_Img(i) - Symbol_Map(1);  % 1st Symbol
            [first_symbol_corr, lags_1, bounds_1]  = crosscorr(Chaos_Seq_Img,chaotic_seq_set{2}, length(chaotic_seq_set{1})-1);
            Chaos_Seq_Img(i) =  Chaos_Seq_Img(i) - Symbol_Map(2);  % 2nd Symbol
            [second_symbol_corr, lags_2, bounds_2] = crosscorr(Chaos_Seq_Img,chaotic_seq_set{2}, length(chaotic_seq_set{1})-1);
            if (first_symbol_corr(length(chaotic_seq_set{1})) > second_symbol_corr(length(chaotic_seq_set{1})) )            
                Partial_Dencoding_Img(i) = Symbol_Map(1);
            else
                Partial_Dencoding_Img(i) = Symbol_Map(2);
            end
    end
    % Complex Decoded Data
    decoded_data(i) = Partial_Dencoding_Real(i) + (sqrt(-1) * Partial_Dencoding_Img(i));
end

end