function [ PILOT_VECTOR ] = PilotDesign( PILOT_NUM, ...
                                         CARRIER_NUM, ...
                                         DEACTIVATED_SUB_CARRIERS, ...
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
PROB_PILOT_VEC = 1:CARRIER_NUM;

% Discard the randomly deactivated Sub-Carriers to hold a Pilot tone:
PROB_PILOT_VEC(DEACTIVATED_SUB_CARRIERS) = [];

PILOT_VECTOR  = zeros(1, PILOT_NUM);
selected_indx = zeros(1, PILOT_NUM);

 switch METHOD
     case 'PILOT_DESIGN_METHOD'
         % Step 1: Obtain p0 and p1 by finding the min {p0, p1} of G1
         initial_indx = 1;
         [ vector_dis ] = CalculateVectorDistances( PROB_PILOT_VEC, ...
                                                    initial_indx, ...
                                                    CARRIER_NUM, ...
                                                    CHANNEL_LENGTH, ...
                                                    METHOD);
         
         vector_dis(initial_indx) = 0;  % Avoids the operations with NaN results.

         [ObjFunct_P_zero_min, ObjFunct_P_zero_indx] = min(vector_dis((initial_indx + 1):end));
         
         second_indx = ObjFunct_P_zero_indx;
         
         % Step 2: Initialize set P = { p_zero, p_one}. Set n = 2:
         selected_indx(1) = initial_indx;
         selected_indx(2) = second_indx;
         
         % Step 3: given{p0, ..., pn}, obtain pn by solving the optimization
         % objective function with respect to that posiition:
         for indx = 2:length(PROB_PILOT_VEC)    
            [ vector_dis_all(indx,:)] = CalculateVectorDistances( PROB_PILOT_VEC, ...
                                                                indx, ...
                                                                CARRIER_NUM, ...
                                                                CHANNEL_LENGTH, ...
                                                                METHOD);
            vector_dis_all(indx,indx) = 0;% Avoids the operations with NaN results.
         end
         
         vector_dis_all(1,:) = vector_dis;
         
         index = 1;
         % Step 4: Continue until kk = PILOT_NUM
         for i = 3:PILOT_NUM
           %Preparing to take into account the selected Pilot positions:
           only_sum_these_positions                   = 1:length(PROB_PILOT_VEC);
           only_sum_these_positions([selected_indx(1:(i-1))]) = [];
           
           base_sum = sum( vector_dis_all( selected_indx(1:(i-1)), :));
           copy_base_sum = base_sum;
           copy_base_sum(selected_indx(1:(i-1))) = [];
           [evaluate_next_min_val index_garbage]= min(copy_base_sum);
           dis_vector_min = evaluate_next_min_val;
           
           for k = 1:length(PROB_PILOT_VEC)
               if (nnz(k == selected_indx) >= 1)
                   dis_vector_min = dis_vector_min;
               else
                    if(dis_vector_min >= base_sum(k))
                        dis_vector_min = base_sum(k);
                        index = k;
                    end
               end
           end
      
           selected_indx(i) = index;
         end
         
         PILOT_VECTOR = sort(PROB_PILOT_VEC(selected_indx));
     case 'SYMETRIC_PILOT_PLACEMENT'
        % Contorl on Pilots: FIXME, nof very accurate for more Pilots > 32
        SEL_PILOT_VECTOR = [1:ceil((CARRIER_NUM - length(DEACTIVATED_SUB_CARRIERS))/PILOT_NUM):(CARRIER_NUM-length(DEACTIVATED_SUB_CARRIERS))];
        
        PILOT_VECTOR = sort(PROB_PILOT_VEC(SEL_PILOT_VECTOR)); 
     otherwise
        error('Not a correct pilot design method supported! (see ''PilotDesign.METHOD'')');    
 end

end

