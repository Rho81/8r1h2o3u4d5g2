% -------------------------------------------------------------------------
% OFDM System Over AWGN channel or TAP Exponential Channel
% -------------------------------------------------------------------------
% Author: Rodrigo Calderon Rico
%         Phd. Candidate,
%         University of Guadalajara,
%         Guadalajara, Mexico.
% -------------------------------------------------------------------------
%
% Description:
%   Rx & Tx Communication blocks for an
%
% Parameters:
%
% Number of Carriers          : CARRIER_NUM
% Coding                      : Convolutional Coding
% Coding Flag                 : CODING_EN
% Single Frame Size(Bits)     : FRAME_SIZE
% Total Number of Frames      : FRAME_NUM 100
% Modulation                  : M_QAM
% Number of Transmit Antennas : N_T
% Number of Receive Antennas  : N_R
% Channel Estimation Enable F : CHANNEL_EST_EN
% Channel Estimation Algrtm T : CH_EST_TYPE
% Channel Type Selection      : CHANNEL_TYPE
% Enable OFDM TX-RX           : OFDM_EN
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
%  Clean Workspace:
% -------------------------------------------------------------------------
close all;
clear all;
clc

% -------------------------------------------------------------------------
%  Simulation Parameters:
% -------------------------------------------------------------------------
% Initialize Random Generators:
% rng('shuffle', 'twister')

% =========================================================================
% Parameters Definitions:
% =========================================================================
SNRdB    = 10;
SNR      = 10.^(SNRdB/10);
BER      = zeros(1,length(SNRdB));
bits     = zeros(1,length(SNRdB));
N_FRAMES = 50;

% -------------------------------------------------------------------------
% OFDM Parameters:
% -------------------------------------------------------------------------
CARRIER_NUM    = 256;    % Number of SubCarriers
CYCLIC_EXT_PER = 25;     % Cyclic Extension Percentaje i.e. 25% (16)

% -------------------------------------------------------------------------
% Modulation Parameters:
% -------------------------------------------------------------------------
ldM       = 2;       % Bits Per Symbol
mod_type  = 2;       % Modulation Type
MAP       = 0;       % Gray == 0

% -------------------------------------------------------------------------
%Channel Parameters:
% -------------------------------------------------------------------------
CHANNEL_TYPE = 'FIR_EXP';  % Channel Options:   
			               %	a) AWGN
                           %    b) FIR_EXP : X Taps Exponential Profile
                           %                 Power.
FIR_TAPS     = 16;

% -------------------------------------------------------------------------
% Channel Estimation Parameters:
% -------------------------------------------------------------------------
CH_EST_TYPE    = 'LS';    % Channel Estimation Methods:
                           %      a) PERFECT_CSI
                           %      b) MSE_LS
                           %      c) LS
                           
if ( strcmp(CH_EST_TYPE, 'PERFECT_CSI') || strcmp(CHANNEL_TYPE, 'AWGN'))                        
    PILOT_NUM   = 0;       % No need of Pilots
else
    PILOT_NUM   = 16;      % Number of Pilots 
end

PILOT_REF   = (1+1i)*sqrt(0.5); % For Channel Est. QPSK Symbols

% Methods Availables:
        %   a)  'SYMETRIC_PILOT_PLACEMENT' : DEACTIVATED_SUB_CARRIERS []
        %   b)  'PILOT_DESIGN_METHOD'      :
        %   c)  'PILOT_FIXED'      :

PILOT_DESIGN_MODE = 'SYMETRIC_PILOT_PLACEMENT';

% -------------------------------------------------------------------------
% MIMO Channel Parameters:  (FUTURE USE)
% -------------------------------------------------------------------------
N_R = 1;
N_T = 1;

% -------------------------------------------------------------------------
%  Cognitive Radio Parameters:
% -------------------------------------------------------------------------
if ( strcmp(CH_EST_TYPE, 'PERFECT_CSI') || strcmp(CHANNEL_TYPE, 'AWGN'))  
   NSubCarriersDeactivated_p = 0;  % XX% of Subcarriers Deactivated
else
   NSubCarriersDeactivated_p = 0.25;  % XX% of Subcarriers Deactivated 
end

N_DeactivatedSubCarriers  = floor(CARRIER_NUM * NSubCarriersDeactivated_p);


tic
fprintf('Simulation Started at %s.\n',datestr(now));

for snr_index = 1:length(SNR);
    
    for tt = 1:N_FRAMES
        
        % -----------------------------------------------------------------
        % Spectrum Sensing Agent:
        % -----------------------------------------------------------------       
        % Methods Availables:
        %   a)  'SYMETRIC_PILOT_PLACEMENT' : DEACTIVATED_SUB_CARRIERS []
        %   b)  'PILOT_DESIGN_METHOD'      : 
        [DSubCarriers] = spectrum_sensing_agnt( CARRIER_NUM, ...
                                                NSubCarriersDeactivated_p, ...
                                                PILOT_DESIGN_MODE);
        
        % -----------------------------------------------------------------
        % Pilot Design Agent:
        % -----------------------------------------------------------------
        % Methods Availables:
        %   a)  'SYMETRIC_PILOT_PLACEMENT' : DEACTIVATED_SUB_CARRIERS []
        %   b)  'PILOT_DESIGN_METHOD'      : 
        if ( strcmp(CH_EST_TYPE, 'PERFECT_CSI') || strcmp(CHANNEL_TYPE, 'AWGN'))
            PILOT_POSITION = [];
        else
            [ PILOT_POSITION ] = PilotDesign(   PILOT_NUM, ...
                                                CARRIER_NUM, ...
                                                DSubCarriers, ...
                                                FIR_TAPS, ...
                                                PILOT_DESIGN_MODE);
        end
        % -----------------------------------------------------------------
        % Generating Data:
        % -----------------------------------------------------------------
        x = floor( randn(1,(CARRIER_NUM - PILOT_NUM - N_DeactivatedSubCarriers) * ldM)) < 1;
        
        % -----------------------------------------------------------------
        % Serial-Parallel Conversion:
        % -----------------------------------------------------------------
        x_parallel = reshape(x, (CARRIER_NUM - PILOT_NUM - N_DeactivatedSubCarriers), ldM);
        
        % Zero_padding for pilot insertion:
        if ( strcmp(CH_EST_TYPE, 'PERFECT_CSI') || strcmp(CHANNEL_TYPE, 'AWGN'))
            x_p = x_parallel;
        else 
            CARRIER_i_data                 = 1:CARRIER_NUM;
            CARRIER_i_data([PILOT_POSITION, DSubCarriers]) = [];   
            x_p                            = zeros(CARRIER_NUM, ldM);
            x_p(CARRIER_i_data,:)          = x_parallel;           
        end
        
        % -----------------------------------------------------------------
        % Modulation
        % -----------------------------------------------------------------
        [x_mod,ldM, maxout] = sigmap_map(x_p,mod_type,MAP);
        
        % -----------------------------------------------------------------
        % Pilot Insertion:
        % -----------------------------------------------------------------
        %  Normally N_PILOT == length(PILOT_POSITION), but if not; 
        %  length(PILOT_POSITION) will have priority, giving the possibility
        %  to change default pilot scheme. Always the first transmission 
        %  will add the PILOT_NUM number of Pilots.
        %  [data_with_pilots, d_pilot]  = addpilots( x_mod,...
        % 					                         PILOT_NUM,...
        % 					                         PILOT_REF,...
        % 				                             PILOT_POSITION,...
        % 					                         N_T,...
        % 					                         CARRIER_NUM);
        %
        x_data_pilot = zeros(CARRIER_NUM, 1);
        if ( strcmp(CH_EST_TYPE, 'PERFECT_CSI') || strcmp(CHANNEL_TYPE, 'AWGN') )
            x_data_pilot = x_mod;           
        else
            % Remove Placeholders to Insert Pilots:
            data_pos = 1:CARRIER_NUM;
            data_pos([PILOT_POSITION, DSubCarriers]) = [];
            x_data_pilot(data_pos)       = x_mod(data_pos);
            x_data_pilot(PILOT_POSITION) = PILOT_REF;
        end
        
        % -----------------------------------------------------------------
        % IFFT
        % -----------------------------------------------------------------
        ifft_sig = ifft(x_mod);
        
        % -----------------------------------------------------------------
        % Adding Cyclic Extension
        % -----------------------------------------------------------------
        % Cyclic Extension:
        CYCLIC_EXT  = ceil((size(ifft_sig,1) * CYCLIC_EXT_PER)/100);  
        cext_data   = zeros(CARRIER_NUM + CYCLIC_EXT, (size(ifft_sig,2)));
        
        cext_data(1:CYCLIC_EXT)     = ifft_sig((end-CYCLIC_EXT+1):end);
        cext_data(CYCLIC_EXT+1:end) = ifft_sig(1:end);

    
        % -----------------------------------------------------------------
        % Parallel-Serial Conversion:
        % -----------------------------------------------------------------
        x_hat   = reshape(cext_data, 1, size(cext_data,1));
        tx_data = x_hat;
    
        % -----------------------------------------------------------------
        % Channel Type Selection:
        % -----------------------------------------------------------------
        H     = zeros(1,FIR_TAPS);     % Channel Initialization
        
        switch CHANNEL_TYPE
            case 'FIR_EXP'
                % Generate channel matrix (Exponential Fading):
                for i = 1:FIR_TAPS
                    H(i) = (randn(N_R,N_T) + sqrt(-1)*randn(N_R,N_T)) * exp((-0.2)*(i));                   
                end
                
                % Power Profile:
                H = abs(H).^2;

                % Signal Filtering:
                % The filter is a direct form II transposed implementation
                % of the standard difference equation.
                tx_data_fir = filter(H,1,tx_data);
                
                noise       = (randn(1,length(tx_data_fir)) + sqrt(-1)*randn(1,length(tx_data_fir))) * (1/sqrt(2));
                
                %y           = tx_data_fir + noise * sqrt(1/(SNR(snr_index) * ldM));
                y           = tx_data_fir + noise * (1/(SNR(snr_index)));
                
                
            case 'AWGN'
                % ---------------------------------------------------------
                % AWGN Channel:
                % ---------------------------------------------------------
                % Uncorrelated Zero-Mean Gaussian Random noise, affected by
                % sqrt(Power Profile) to have equal energy in real and 
                % imaginary part:
                noise = (randn(1,length(tx_data)) + sqrt(-1)*randn(1,length(tx_data))) * 1/sqrt(2);
                y = tx_data + noise * (1/(SNR(snr_index) * ldM ));
        end


        % -----------------------------------------------------------------
        %                            RECEIVER
        % -----------------------------------------------------------------
        
        % -----------------------------------------------------------------
        % Serial-Parallel Conversion:
        % -----------------------------------------------------------------
        y_parallel = reshape(y, size(y,2), 1);
        
        % -----------------------------------------------------------------
        % Removing Cyclic Extension
        % -----------------------------------------------------------------
        y_cycle = y_parallel(CYCLIC_EXT+1:end);
        
        % -----------------------------------------------------------------
        % FFT
        % -----------------------------------------------------------------
        ff_sig = fft(y_cycle);
                
        % -----------------------------------------------------------------
        % Channel Estimation:
        % -----------------------------------------------------------------
        if ( strcmp(CH_EST_TYPE, 'PERFECT_CSI') && ~strcmp(CHANNEL_TYPE, 'AWGN'))
            % Perfect channel knowledge + Zero_Padding, FFT Interpolation used:          
            [ H_est ] = channel_est_rx(H, SNR, ff_sig, 'PERFECT_CSI', PILOT_POSITION, PILOT_REF, CARRIER_NUM);
            
            % -------------------------------------------------------------
            % Equalization:
            % -------------------------------------------------------------
            for kk = 1:CARRIER_NUM
                % Equalization Per-SubCarrier:
                y_hat(kk,:) = ff_sig(kk,:)./H_est(kk);
            end
                
        elseif (strcmp(CHANNEL_TYPE, 'AWGN'))
            y_hat = ff_sig;
            
        else
            % -------------------------------------------------------------
            % Get Pilots:
            % -------------------------------------------------------------
            H_Pilot = ff_sig(PILOT_POSITION);
            
            if ( strcmp(CH_EST_TYPE, 'LS') )
                [ H_est ] = channel_est_rx(H_Pilot, SNR(snr_index), ff_sig, 'LS', PILOT_POSITION, PILOT_REF, CARRIER_NUM);
            elseif ( strcmp(CH_EST_TYPE, 'MSE_LS') )
                [ H_est ] = channel_est_rx(H_Pilot, SNR(snr_index), ff_sig, 'MSE_LS', PILOT_POSITION, PILOT_REF, CARRIER_NUM);
            end
            
            % -------------------------------------------------------------
            % Equalization:
            % -------------------------------------------------------------
            for kk = 1:CARRIER_NUM
                y_hat(kk,:) = ff_sig(kk,:)./H_est(kk);
            end
            
            
        end
    
        % -----------------------------------------------------------------
        % Demodulation:
        % -----------------------------------------------------------------
        y_hat([DSubCarriers PILOT_POSITION]) = [];
        [soft_dem_data, dem_data]            = sigdemap(y_hat,mod_type,MAP);
                
        % -----------------------------------------------------------------
        % Parallel-Serial Conversion:
        % -----------------------------------------------------------------
        if ( strcmp(CH_EST_TYPE, 'PERFECT_CSI') || strcmp(CHANNEL_TYPE, 'AWGN'))
            dem_data_s = reshape(dem_data, 1, size(dem_data,2)*(CARRIER_NUM));
        else
            dem_data_s = reshape(dem_data, 1, size(dem_data,2)*(CARRIER_NUM - PILOT_NUM - N_DeactivatedSubCarriers));
        end
        
        % -----------------------------------------------------------------
        % Calculating BER
        % -----------------------------------------------------------------
        % Remove from x the pilot positions:

        bit_diff = zeros(1,length(x));
        for bit_index = 1:length(x)
            if ( x(bit_index) ~= dem_data_s(bit_index) )
                bit_diff(bit_index) = 1;
            else
                bit_diff(bit_index) = 0;
            end
        end
        errors    = nnz(bit_diff);
        BER(snr_index)    =  BER(snr_index) + errors;
        bits(snr_index)   =  bits(snr_index) + length(dem_data_s);
       
    end
    
end

% -----------------------------------------------------------------------
% Graphic BER
% -----------------------------------------------------------------------

ber = BER./bits;
figure(100)
semilogy(SNRdB,ber,'ko-', 'linewidth' , 1.0);
% This only to check if the proyection is correct based on its distribution
% hold on
% semilogy(SNRdB, (0.5 * erfc(sqrt(SNR)/sqrt(2)) ),'m-','linewidth',2.0);
title('BER vs SNR');
ylabel('BER');
xlabel('SNR (dB)');
grid on

% -------------------
% Equalization loos in power
% -------------------
% figure(905);
% plot(SNRdB,channel_ideal_power, 'k-', 'linewidth' , 1.0);
% hold on
% plot(channel_measured_power, 'b-', 'linewidth' , 1.0 );
% title('Channel Ideal Power vs Channel Estimated Power');
% grid on

if (length(SNRdB) == 1)
    H_fft = abs(fft(H, CARRIER_NUM));
    H_pp  = zeros(1, CARRIER_NUM);
    H_pp(PILOT_POSITION) = H_est(PILOT_POSITION);
    
    %H_Deactivated = zeros(1:CARRIER_NUM);
    H_Deactivated(DSubCarriers) = abs(H_est(DSubCarriers));

    figure(1);
    semilogy(H_fft, 'o', 'MarkerSize', 6);
    grid on
    hold on
    semilogy(abs(H_pp), 'kp', 'MarkerSize', 8, 'MarkerFaceColor','k');
    hold all
    semilogy(abs(H_est),'m*');
    hold on
    semilogy(H_Deactivated,'g*');
    hleg = legend('Perfect Channel Knowledge','Selected Pilot Tone Positions','Interpolated Estimated Channel', 'Unused Sub-Carrriers');
    xlabel('Number of OFDM Subcarriers');
    ylabel('ABS(H)');
    
    
    figure(2);
    semilogy(real(fft(H,CARRIER_NUM)), 'o', 'MarkerSize', 6);
    grid on
    hold on
    semilogy(real(H_pp), 'kp', 'MarkerSize', 8, 'MarkerFaceColor','k');
    hold all
    semilogy(real(H_est),'m*');
    hold on
    semilogy(real(H_est(DSubCarriers)),'g*');
    hleg = legend('Perfect Channel Knowledge','Selected Pilot Tone Positions','Interpolated Estimated Channel', 'Unused Sub-Carrriers');
    xlabel('Number of OFDM Subcarriers');
    ylabel('Real(H)');
    
    figure(3);
    semilogy(imag(fft(H,CARRIER_NUM)), 'o', 'MarkerSize', 6);
    grid on
    hold on
    semilogy(imag(H_pp), 'kp', 'MarkerSize', 8, 'MarkerFaceColor','k');
    hold all
    semilogy(imag(H_est),'m*');
    hold on
    semilogy(imag(H_est(DSubCarriers)),'g*');
    hleg = legend('Perfect Channel Knowledge','Selected Pilot Tone Positions','Interpolated Estimated Channel', 'Unused Sub-Carrriers');
    xlabel('Number of OFDM Subcarriers');
    ylabel('Imaginary(H)');
    
    
else
%     figure(1);
%     semilogy(channel_ideal_power(SNRdB,:), 'o', 'MarkerSize', 6);
%     grid on
%     hold on
%     semilogy(channel_measured_power(SNRdB,:), 'kp', 'MarkerSize', 8, 'MarkerFaceColor','k');
%     hleg = legend('Perfect Channel Knowledge','Selected Pilot Tone Positions');
%     xlabel('Number of OFDM Subcarriers');
%     ylabel('ABS(H)');
end

toc
fprintf('Simulation ended at %s.\n',datestr(now));

%% EOF