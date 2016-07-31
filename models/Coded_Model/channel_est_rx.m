% -- * -- * -- * -- * -- * -- * -- * -- * -- * -- * -- * -- * -- * -- * -- *
% Rho Corporation @ 2016 - PhD Communications
% -- * -- * -- * -- * -- * -- * -- * -- * -- * -- * -- * -- * -- * -- * -- * 
%
% -------------------------------------------------------------------------
% Author: Rodrigo Calderon Rico
%         Phd. Candidate,
%         University of Guadalajara,
%         Guadalajara, Mexico.
% -------------------------------------------------------------------------
%
% Description:  This function is aimed to provide the estimated channel,
% based on the received data in Freq-Domain and the a-priori knowledge of
% the pilot tones.
%
%     
% Inputs:
%
%          d_pilot        = channel information in the form of pilot tones
%          sigma_n        = Signal-to-Noise-Ratio (SNR)
%          data           = received data in Frequency-Domain
%          ce_type        = Channel Estimation Methods:
%                            a) PERFECT_CSI
%                            b) MSE_LS
%                            c) LS
%          PILOT_POSITION = Pilot Tone´s Position
%          PILOT_REF      = Reference Pilot Tone
%          CARRIER_NUM    = Number of Sub-Carriers
%
% Outputs:
%
%          H_est  = Interpolated Estimated Channel
% -------------------------------------------------------------------------

function [ H_est ] = channel_est_rx(d_pilot, sigma_n, data, ce_type, PILOT_POSITION, PILOT_REF, CARRIER_NUM)

    switch ce_type
        case 'PERFECT_CSI'
            % Perfect channel knowledge + Zero_Padding
            H_hat = [d_pilot zeros(1,(length(data) - length(d_pilot)))];     
            H_est = fft(H_hat);
            %figure(5000);
            %plot(H_est, 'x', 'MarkerSize', 8);
        case 'LS'
           
            % Channel knowledge + Zero_Padding
            H_hat_pos  = 1:CARRIER_NUM+1;
            
            H_est      = d_pilot./PILOT_REF; 
            
            % Interpolate by FFT method:
%             H_hat_interp = ifft(H_est, (length(data)- length(d_pilot)));
%             H_hat_interp = fft(H_hat_interp);
            
            H_hat_interp = interp1(PILOT_POSITION, H_est, H_hat_pos, 'spline', 'extrap');
            
            H_est = H_hat_interp(1:256);          
        case 'MSE_LS'
            % Channel knowledge + Zero_Padding
            H_hat_pos  = 1:CARRIER_NUM+1;

            H_est   = (sigma_n^2 / PILOT_REF) * d_pilot;
            
            % Interpolate Method:        
            H_hat_interp = interp1(PILOT_POSITION,H_est, H_hat_pos,'spline','extrap');
            
            H_est        = H_hat_interp(1:256);           
        case 'MMSE'
            H_est   = data*d_pilot'*inv(sigma_n^2*eye(size(d_pilot,1))+d_pilot*d_pilot');
        case 'ML'
            H_est   = data*d_pilot'*inv(d_pilot*d_pilot');
        otherwise
            error('Not a correct channel estimation type specified! (see ''ce_type'')');
    end
end
