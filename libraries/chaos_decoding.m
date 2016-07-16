function  [ decoded_data ] = chaos_decoding(coded_sequence, chaotic_sequence)

% Chaotic Sequence Length:
N = length(coded_sequence);

% Chaotic Sequences:
Chaos_Seq_Real    = zeros(1,N);
Chaos_Seq_Img     = zeros(1,N);

for i=1:N
   Chaos_Seq_Real(i) = real(chaotic_sequence(i));
   Chaos_Seq_Img(i)  = imag(chaotic_sequence(i)); 
end

% Data to be encoded
decoded_data = zeros(1,N);

% Decoding Process:
for i=1:N
    decoded_data(i) = (Chaos_Seq_Real(i) + real(coded_sequence(i))) + (sqrt(-1) * (Chaos_Seq_Img(i) + imag(coded_sequence(i))));
end

end