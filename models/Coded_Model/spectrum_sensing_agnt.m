function [DSubC] = spectrum_sensing_agnt(   proposed_inactive_pilot_spaces, ...
                                            CARRIER_NUM, ...
                                            DeactivatedSubPercent, ...
                                            METHOD)
% -------------------------------------------------------------------------
% Author: Rodrigo Calderon Rico
%         Phd. Candidate,
%         University of Guadalajara,
%         Guadalajara, Mexico.
% -------------------------------------------------------------------------
% Description Goes Here....
% Methods Availables:
%   a)  'SYMETRIC_PILOT_PLACEMENT' : DEACTIVATED_SUB_CARRIERS []
%   b)  'PILOT_DESIGN_METHOD'      : 
% -------------------------------------------------------------------------
% Initializate:
N_DeactivatedSubCarriers   = floor (CARRIER_NUM * DeactivatedSubPercent);
Selected_Items             = zeros(1,N_DeactivatedSubCarriers);

switch METHOD
    case 'PILOT_DESIGN_METHOD'
        
        if (N_DeactivatedSubCarriers > 0)
            
            PilotPositionsSet = 1:CARRIER_NUM;
            
            for i = 1:N_DeactivatedSubCarriers
               random_items         = randperm(CARRIER_NUM - (i - 1));
               Selected_Items(i)    = PilotPositionsSet(random_items(1));
               PilotPositionsSet(random_items(1)) = [] ;   % Remove Selected Item from Set
            end

        else
            Selected_Items  = [];
        end

        DSubC = sort(Selected_Items);
    
    case 'CHAOTIC_PILOT_DESIGN_METHOD'
        
        if (N_DeactivatedSubCarriers > 0)
            
            PilotPositionsSet = 1:CARRIER_NUM;
            
            for i = 1:N_DeactivatedSubCarriers
               random_items         = randperm(CARRIER_NUM - (i - 1));
               Selected_Items(i)    = PilotPositionsSet(random_items(1));
               PilotPositionsSet(random_items(1)) = [] ;   % Remove Selected Item from Set
            end

        else
            Selected_Items  = [];
        end

        DSubC = sort(Selected_Items);

    case 'FULL_WAVELET_PILOT_DESIGN_METHOD'
        
        if (N_DeactivatedSubCarriers > 0)
            
            PilotPositionsSet = 1:CARRIER_NUM;
            
            for i = 1:N_DeactivatedSubCarriers
               random_items         = randperm(CARRIER_NUM - (i - 1));
               Selected_Items(i)    = PilotPositionsSet(random_items(1));
               PilotPositionsSet(random_items(1)) = [] ;   % Remove Selected Item from Set
            end

        else
            Selected_Items  = [];
        end

        DSubC = sort(Selected_Items);        
        
     case 'WAVELET_ENERGY_BASED'
        
        if (N_DeactivatedSubCarriers > 0)
            
            PilotPositionsSet = 1:CARRIER_NUM;
            
            for i = 1:N_DeactivatedSubCarriers
               random_items         = randperm(CARRIER_NUM - (i - 1));
               Selected_Items(i)    = PilotPositionsSet(random_items(1));
               PilotPositionsSet(random_items(1)) = [] ;   % Remove Selected Item from Set
            end

        else
            Selected_Items  = [];
        end

        DSubC = sort(Selected_Items);
        
     case 'USE_FEEDBACK_CHANNEL_ESTIMATION'
        
        if (N_DeactivatedSubCarriers > 0)
            
            PilotPositionsSet = 1:CARRIER_NUM;
            PilotPositionsSet = PilotPositionsSet(proposed_inactive_pilot_spaces);
            
            for i = 1:N_DeactivatedSubCarriers
               random_items         = randperm(length(PilotPositionsSet));
               Selected_Items(i)    = PilotPositionsSet(random_items(1));
               PilotPositionsSet(random_items(1)) = [] ;   % Remove Selected Item from Set
            end

        else
            Selected_Items  = [];
        end

        DSubC = sort(Selected_Items);

    case 'SYMETRIC_PILOT_PLACEMENT'
        if (N_DeactivatedSubCarriers > 0)
            
            PilotPositionsSet = 1:CARRIER_NUM;
            
            for i = 1:N_DeactivatedSubCarriers
               random_items         = randperm(CARRIER_NUM - (i - 1));
               Selected_Items(i)    = PilotPositionsSet(random_items(1));
               PilotPositionsSet(random_items(1)) = [] ;   % Remove Selected Item from Set
            end

        else
            Selected_Items  = [];
        end

        DSubC = sort(Selected_Items);
        
    case 'ASSOCIATIVE_PILOT_ASSIGMENT_WITH_FEEDBACK_CHANNEL_INFORMATION'
        
        if (N_DeactivatedSubCarriers > 0)
            
            PilotPositionsSet = 1:CARRIER_NUM;
            
            for i = 1:N_DeactivatedSubCarriers
               random_items         = randperm(CARRIER_NUM - (i - 1));
               Selected_Items(i)    = PilotPositionsSet(random_items(1));
               PilotPositionsSet(random_items(1)) = [] ;   % Remove Selected Item from Set
            end

        else
            Selected_Items  = [];
        end

        DSubC = sort(Selected_Items);
        
    case 'SYMETRIC_PILOTS_SET_AND_INFLUENCE_OF_HIGHEST_CH_VARIATIONS'
        
        if (N_DeactivatedSubCarriers > 0)
            
            PilotPositionsSet = 1:CARRIER_NUM;
            
            for i = 1:N_DeactivatedSubCarriers
               random_items         = randperm(CARRIER_NUM - (i - 1));
               Selected_Items(i)    = PilotPositionsSet(random_items(1));
               PilotPositionsSet(random_items(1)) = [] ;   % Remove Selected Item from Set
            end

        else
            Selected_Items  = [];
        end

        DSubC = sort(Selected_Items);
        
    case 'PILOTS_ONLY_ON_THE_HIGHEST_CH_VARIATIONS'
        if (N_DeactivatedSubCarriers > 0)
            
            PilotPositionsSet = 1:CARRIER_NUM;
            
            for i = 1:N_DeactivatedSubCarriers
               random_items         = randperm(CARRIER_NUM - (i - 1));
               Selected_Items(i)    = PilotPositionsSet(random_items(1));
               PilotPositionsSet(random_items(1)) = [] ;   % Remove Selected Item from Set
            end

        else
            Selected_Items  = [];
        end

        DSubC = sort(Selected_Items);
        
    case 'PILOT_FIXED'

        if (N_DeactivatedSubCarriers <= 16)
            DeactivatedSubCarriers = [14 15 10 5 6 16 11 3 12 2 ...
                                      4 9 1 13 7 8];
            DeactivatedSubCarriers = DeactivatedSubCarriers(1:N_DeactivatedSubCarriers)+1-1;
        elseif ( N_DeactivatedSubCarriers <= 32 )
            DeactivatedSubCarriers = [24 18 32 20 8 12 7 5 21 19 22 ...
                                      1 6 29 15 31 17 23 3 28 27 16 ...
                                      11 10 14 2 13 9 30 25 4 26];
            DeactivatedSubCarriers = DeactivatedSubCarriers(1:N_DeactivatedSubCarriers)+1-1;

        elseif ( N_DeactivatedSubCarriers <= 64 )
            DeactivatedSubCarriers = [59 35 28 27 55 48 63 57 32 4 5 ...
                                      51 37 41 13 49 10 14 8  6 43 22 ...
                                      39 36 12 17 25 64 56 47 34 16 40 ...
                                      29 53 3 20 26 33 19 42 15 44 45 ...
                                      46 24 23 60 30 38 9 61 52 18 7 ...
                                      62 1 50 21 11 31 2 58 54 ];
            DeactivatedSubCarriers = DeactivatedSubCarriers(1:N_DeactivatedSubCarriers)+1-1;

        else
            DeactivatedSubCarriers  = [];
            error('N_DeactivatedSubCarriers Specified, is not: 16 32 or 64...');
        end

        DSubC = sort(DeactivatedSubCarriers);
            
    otherwise
        error('No METHOD Specified, see: spectrum_sensing_agnt.METHOD');
end
%% EOF
