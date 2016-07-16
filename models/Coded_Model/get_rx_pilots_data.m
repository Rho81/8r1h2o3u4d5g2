% -------------------------------------------------------------------------
% Author: Rodrigo Calderon Rico
%         Phd. Candidate,
%         University of Guadalajara,
%         Guadalajara, Mexico.
% -------------------------------------------------------------------------
function [rx_pilots, rx_d] = get_rx_pilots_data(y, N_Pilot, Pilot_Position, N_T, N_R, CARRIER_NUM)
% get_rx_pilots_data: Rx Pilots and Rx Data is received into one single
% rx vector based on the Pilot_Position vector indication.
%   y              == Received Data (Pilots + Data)
%   N_Pilot        == Number of Pilots Being Inserted
%   Pilot_Position == It Contains the Pilot Location Reference 
%   N_T            == Number of Tx Antennas
%   N_R            == Number of Rx Antennas
%   CARRIER_NUM    == 

% Getting Pilot Length:
[row_i, rx_pilot_length]  = size(Pilot_Position);

% Getting the Data Length:
[row_i, data_length] = size(y');
rx_data_length       = data_length - rx_pilot_length;

% ------------------------------------------------------------------------
% Checking for Errors:
% ------------------------------------------------------------------------
if (N_Pilot > CARRIER_NUM) || ( rx_pilot_length > CARRIER_NUM)
   error('addpilots.m: Number of Pilots Exceed the Number of Carriers!');
end

if ( (rx_pilot_length + rx_data_length) > CARRIER_NUM) || ...
     ((N_Pilot + rx_data_length) > CARRIER_NUM)
   error('addpilots.m: Number of Pilots + Data Symbols Exceed the Number of Carriers!');
end

% Placeholder Vector to Store Rx Pilots: 
pilots    = zeros(1,N_Pilot);
pilot_cnt = 0;
data_cnt  = 0;

% ------------------------------------------------------------------------
% Get Pilots:
% ------------------------------------------------------------------------
for k = 1:data_length
    % By looking every time at all pilot positions
    % we ensure that if the Pilot_Position vector
    % is out of order we locate right pilots:
    for kk = 1:rx_pilot_length
        if ( Pilot_Position(kk) == k)
            rm_pilot_flag = 1;
            pilot_cnt     = pilot_cnt + 1;
            break;
        else
            rm_pilot_flag = 0;
        end
    end

    if (rm_pilot_flag)
        pilots(pilot_cnt) = y(k);
    else
        data_cnt  = data_cnt + 1;
    end
         
end

% ------------------------------------------------------------------------
% Checking for Errors:
% ------------------------------------------------------------------------
if (rx_data_length  ~=  data_cnt)
   error('get_rx_pilots_data.m: Final data length is different from input data length!');
end
if (rx_pilot_length ~=  pilot_cnt)
   error('get_rx_pilots_data.m: Final data is different from input data!');
end

% Remove Pilots from Rx Data:
% Separate from get pilots loop in order to be simple in removing the 
% extra data information for rx_d:
for jj =  Pilot_Position
   y(jj) = []; 
end
             
rx_d      = y;
rx_pilots = pilots;
end

