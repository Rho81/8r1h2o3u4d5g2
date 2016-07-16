function [ ObjFunct_P_MN] = PilotDesignObjFunct( PILOT_M, ...
                                                 PILOT_N, ...
                                                 CARRIER_NUM, ...
                                                 CHANNEL_LENGTH, ...
                                                 METHOD)
% -------------------------------------------------------------------------
% Author: Rodrigo Calderon Rico
%         Phd. Candidate,
%         University of Guadalajara,
%         Guadalajara, Mexico.
% -------------------------------------------------------------------------
% Description goes here...
%
% -------------------------------------------------------------------------
 switch METHOD
     case 'PILOT_DESIGN_METHOD'
        ObjFunct_P_MN = (sin((pi * CHANNEL_LENGTH * (PILOT_M - PILOT_N)) / CARRIER_NUM ))^2 ...
                      / (sin( (pi * (PILOT_M - PILOT_N)) / CARRIER_NUM))^2;
     otherwise
        error('Not a correct pilot design method supported! (see ''PilotDesignObjFunct.METHOD'')');    
 end

end