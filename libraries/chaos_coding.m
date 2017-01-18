% -- * -- * -- * -- * -- * -- * -- * -- * -- * -- * -- * -- * -- * -- *
% Rho Corporation @ 2016 - PhD Communications
% -- * -- * -- * -- * -- * -- * -- * -- * -- * -- * -- * -- * -- * -- * 
function  [ encoded_data,  ...
            Chaos_Sequence     ] = chaos_coding(  IQ_sequence,     ...
                                                  sequence_length, ...
                                                  Mod_Symbol_Map,  ...
                                                  display_on)
% Description:
% 
%  chaos_coding encodes complex IQ data based on the following sequences:
%  Chaotic Equation Setup:  
%      x_(n+1) = f(x_n) = 1 - mu * (x_n)^2 , mu ? [1.40015,2], x ? (-1,1)
%
%      with the initial conditions ex.:
%  
%            x_o = 0.655007;   % Real (I)
%            x_1 = 0.655008;   % Imaginary (Q)
%
%  Inputs:
%      IQ_sequence     = complex IQ-data
%      sequence_length = Chaotic sequence length to be generated = N
%      Mod_Symbol_Map  = Modulation Symbol Map Type
%      display_on      = enables the chaotic sequence preview
%
%  Outputs:
%      encoded_data     = complex IQ-data with the chaotic offset included
%      chaotic_sequence = This complex IQ-data it is meant to be used for
%                         decoding.

unit_constant = 1;

% Bifurcation Parameter:
mu = 2;

% Chaotic Sequence Length:
N = sequence_length;

% IQ-Data Length:
IQ_Lenght = length(IQ_sequence);

switch Mod_Symbol_Map
    case 'QPSK'
        N_Chaotic_Sequences = 2;
    otherwise
        % QPSK Modulation
        N_Chaotic_Sequences = 2;
end

% Initial Value:
x_o = rand(1, N_Chaotic_Sequences);
%x_o = [ 0.655007, 0.655008 ];

% Chaotic Sequences Initial Conditions:
parfor i = 1:N_Chaotic_Sequences
    Chaos_Sequence{i} = zeros(1,N);
end

% Iterative Process to generate Chaotic Sequences:
parfor i = 1:N_Chaotic_Sequences
    Chaos_Seq    = Chaos_Sequence{i};
    Chaos_Seq(1) = x_o(i);
    for j=1:N-1
       Chaos_Seq(j+1) = (unit_constant - (mu *  Chaos_Seq(j)^2));
    end
    Chaos_Sequence{i} = Chaos_Seq;
end

% Encoding Process:
Chaos_Seq_Real   = zeros(1,N);
Chaos_Seq_Img    = zeros(1,N);
Partial_Encoding = zeros(1,N);
for i=1:IQ_Lenght
    switch Mod_Symbol_Map
        case 'QPSK'
            % Real Part (I)
            Chaos_Seq_Real = ((real(IQ_sequence(i))) + (Chaos_Sequence{1}(i)));
            % Imaginary Part (Q)
            Chaos_Seq_Img  = sqrt(-1) * (((imag(IQ_sequence(i))) + (Chaos_Sequence{2}(i)))); % Complex Value
    end
    Partial_Encoding(i) = Chaos_Seq_Real + Chaos_Seq_Img;
end

% Assignments:
encoded_data = Partial_Encoding;

if(display_on == 1)
    parfor i=1:N_Chaotic_Sequences
        % Auto-Correlation Calculations:
        acf_funct   = autocorr(Chaos_Sequence{i},N-1);

        figure(i)
        subplot(1,2,1);
        plot(acf_funct);
        grid on
        legend('Auto-Correlation Function');
        title('Auto-Correlation Function of Chaotic Sequence');
        ylabel('Correlation Percentage');
        xlabel('Number of Samples');

        subplot(1,2,2);
        for k=1:N
            plot(encoded_data(k),'o','MarkerSize',6);
            hold on
            grid on
        end
        legend('Complex Mapping (Modulation)');
        title('Coded Data with Chaotic Sequences');
        ylabel('Img');
        xlabel('Real');
    end
    figure(3000);
    crosscorr(Chaos_Sequence{1},Chaos_Sequence{2},N-1);
    figure(4000);
    crosscorr(real(Partial_Encoding), Chaos_Sequence{1}, N-1);
end

end