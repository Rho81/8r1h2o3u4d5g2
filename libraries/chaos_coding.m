function  [encoded_data, chaotic_sequence] = chaos_coding(uncoded_sequence, display_on)

% Chaotic Equation Setup:  x_(n+1) = f(x_n) = 1 - mu * (x_n)^2 , mu ? [1.40015,2], x ? (-1,1) 
unit_constant = 1;

% Bifurcation Parameter:
mu = 2;

% Chaotic Sequence Length:
N = length(uncoded_sequence);

% Initial Value:
x_o = 0.655007;   % Real
x_1 = 0.655008;   % Imaginary

% Chaotic Sequences:
Chaos_Seq_Real    = zeros(1,N);
Chaos_Seq_Real(1) = x_o;

Chaos_Seq_Img     = zeros(1,N);
Chaos_Seq_Img(1)  = x_1;

% Iterative Process to generate Chaotic Sequence:
for i=1:N
   Chaos_Seq_Real(i+1) = unit_constant - (mu *  Chaos_Seq_Real(i)^2);
   Chaos_Seq_Img(i+1)  = unit_constant - (mu *  Chaos_Seq_Img(i)^2); 
end

coded_data       = zeros(1,N);
chaotic_sequence = zeros(1,N);
% Encoding Process:
for i=1:N
    coded_data(i) = (Chaos_Seq_Real(i) + real(uncoded_sequence(i))) + (sqrt(-1) * (Chaos_Seq_Img(i) + imag(uncoded_sequence(i))));
    chaotic_sequence(i) = (-1) * (Chaos_Seq_Real(i) + ((sqrt(-1) * Chaos_Seq_Img(i))));
end

% Assignments:
encoded_data = coded_data;

if(display_on == 1)
    % Calculations:
    seq_autocorr_real = acf(Chaos_Seq_Real(:),N);
    seq_autocorr_img  = acf(Chaos_Seq_Img(:),N);

    % Re-mapping:
    Sign_val = zeros(1,N);
    for i=1:N
        Sign_val(i) = sign(Chaos_Seq_Real(i)); 
    end
    
    figure(1)
    subplot(3,2,1);
    plot(Chaos_Seq_Real);

    subplot(3,2,2);
    plot(seq_autocorr_real);

    subplot(3,2,3);
    plot(seq_autocorr_img);

    subplot(3,2,4);
    for i=1:N
        plot(Sign_val(i),'o','MarkerSize',6);
        hold on
        grid on
    end

    subplot(3,2,5);
    for i=1:N
        plot(uncoded_sequence(i),'o','MarkerSize',6);
        hold on
        grid on
    end

    subplot(3,2,6);
    for i=1:N
        plot(coded_data(i),'o','MarkerSize',6);
        hold on
        grid on
    end
end

end