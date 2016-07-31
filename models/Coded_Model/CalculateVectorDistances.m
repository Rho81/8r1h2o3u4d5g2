% -------------------------------------------------------------------------
% Author: Rodrigo Calderon Rico
%         		Phd. Candidate,
%         		University of Guadalajara,
%         		Guadalajara, Mexico.
% ---------------------------------------------------------------------------
% Description: Calculate vector distances with respect to
% the first element based on PilotDesignObjFunct and 
% METHOD specified.
% ---------------------------------------------------------------------------
function [  vector_dis ] = CalculateVectorDistances( vector, selected_indx, ...
                                                     CARRIER_NUM, ...
                                                     CHANNEL_LENGTH)
		vector_dis = zeros( 1, length(vector) );
         for i = 1:length(vector)
         		vector_dis( i ) = PilotDesignObjFunct(  vector( selected_indx ), ...
         												vector( i ), ...
                                                        CARRIER_NUM, ...
                                                        CHANNEL_LENGTH);
         end
end