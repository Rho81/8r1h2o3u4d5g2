% -- * -- * -- * -- * -- * -- * -- * -- * -- * -- * -- * -- * -- * -- *
% Rho Corporation @ 2016 - PhD Communications
% -- * -- * -- * -- * -- * -- * -- * -- * -- * -- * -- * -- * -- * -- * 

function  [ decoded_data ] = chaos_decoding_1( coded_sequence,   ...
                                               Tx_Symbol_N,      ...
                                               chaotic_seq_set,  ...
                                               Mod_Symbol_Map       )
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
        N_Chaotic_Sequences = 2;
        Symbol_Map          = [ sqrt(0.5), (-1) * sqrt(0.5) ];
    otherwise
        % QPSK Modulation
        N_Chaotic_Sequences = 2;
        Symbol_Map          = [ sqrt(0.5), (-1) * sqrt(0.5) ];
end

% Separate IQ Signals:
Chaos_Seq_Real    = real(coded_sequence);
Chaos_Seq_Img     = imag(coded_sequence);

% Data to be encoded
decoded_data = zeros(1,N);

% Decoding Process:
Partial_Dencoding_Real = zeros(1,N);
Partial_Dencoding_Img  = zeros(1,N);

for i=1:Tx_Symbol_N
    switch Mod_Symbol_Map
        case 'QPSK'
            % Real Part (I)
            [first_symbol_corr, lags_1, bounds_1]  = crosscorr(Chaos_Seq_Real(i,:),chaotic_seq_set{1}, length(chaotic_seq_set{1})-1);
            [second_symbol_corr, lags_2, bounds_2] = crosscorr(Chaos_Seq_Real(i,:),chaotic_seq_set{2}, length(chaotic_seq_set{1})-1);
            if (first_symbol_corr(length(chaotic_seq_set{1})) > second_symbol_corr(length(chaotic_seq_set{1})) )
                Partial_Dencoding_Real(i) = Symbol_Map(1);
            else
                Partial_Dencoding_Real(i) = Symbol_Map(2);
            end
            % Imaginary Part (Q)
            [first_symbol_corr, lags_1, bounds_1]  = crosscorr(Chaos_Seq_Img(i,:),chaotic_seq_set{1}, length(chaotic_seq_set{1})-1);
            [second_symbol_corr, lags_2, bounds_2] = crosscorr(Chaos_Seq_Img(i,:),chaotic_seq_set{2}, length(chaotic_seq_set{1})-1);
            if (first_symbol_corr(length(chaotic_seq_set{1})) > second_symbol_corr(length(chaotic_seq_set{1})) )
                Partial_Dencoding_Img(i) = Symbol_Map(1);
            else
                Partial_Dencoding_Img(i) = Symbol_Map(2);
            end
    end
    % Complex Decoded Data
    decoded_data(i) = Partial_Dencoding_Real(i) + sqrt(-1) * Partial_Dencoding_Img(i);
end

end