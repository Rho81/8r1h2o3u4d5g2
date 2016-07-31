function [ ObjFunct_P_MN] = PilotDesignObjFunct( PILOT_M, ....
                                                 PILOT_N, ....
                                                 CARRIER_NUM, ...
                                                 CHANNEL_LENGTH)
% -------------------------------------------------------------------------
% Author: Rodrigo Calderon Rico
%         Phd. Candidate,
%         University of Guadalajara,
%         Guadalajara, Mexico.
% -------------------------------------------------------------------------
% Description goes here...
%
% -------------------------------------------------------------------------

ObjFunct_P_MN = (sin((pi * CHANNEL_LENGTH * (PILOT_M - PILOT_N)) / CARRIER_NUM ))^2 ...
                / (sin( (pi * (PILOT_M - PILOT_N)) / CARRIER_NUM))^2;
     

end