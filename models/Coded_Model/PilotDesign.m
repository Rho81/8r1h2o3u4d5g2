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
% Description goes here...
%
% -------------------------------------------------------------------------

function [ PILOT_VECTOR ] = PilotDesign( PILOT_NUM, ...
                                         CARRIER_NUM, ...
                                         DEACTIVATED_SUB_CARRIERS, ...
                                         CHANNEL_LENGTH, ...
                                         METHOD, ...
                                         proposed_active_pilot_spaces, ...
                                         FORCE_PILOT_NUM)

PROB_PILOT_VEC = 1:CARRIER_NUM;

PILOT_VECTOR  = zeros(1, PILOT_NUM);
selected_indx = zeros(1, PILOT_NUM);

 switch METHOD
     case 'PERFECT_CSI_KNOWLEDGE'
         PILOT_VECTOR = [];

     case 'PILOT_DESIGN_METHOD'
         % Discard the randomly deactivated Sub-Carriers to hold a Pilot tone:
         PROB_PILOT_VEC(DEACTIVATED_SUB_CARRIERS) = [];
         
         % Step 1: Obtain p0 and p1 by finding the min {p0, p1} of G1
         initial_indx = 1;
         [ vector_dis ] = CalculateVectorDistances( PROB_PILOT_VEC, ...
                                                    initial_indx, ...
                                                    CARRIER_NUM, ...
                                                    CHANNEL_LENGTH);
         
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
                                                                  CHANNEL_LENGTH);
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
     
     case 'USE_FEEDBACK_CHANNEL_ESTIMATION'
         
         % Discard the randomly deactivated Sub-Carriers to hold a Pilot tone:
         PROB_PILOT_VEC(DEACTIVATED_SUB_CARRIERS) = [];
         
         % Step 1: Obtain p0 and p1 by finding the min {p0, p1} of G1
         initial_indx = 1;
         [ vector_dis ] = CalculateVectorDistances( PROB_PILOT_VEC, ...
                                                    initial_indx, ...
                                                    CARRIER_NUM, ...
                                                    CHANNEL_LENGTH);
         
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
                                                                  CHANNEL_LENGTH);
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
        % Discard the randomly deactivated Sub-Carriers to hold a Pilot tone:
        PROB_PILOT_VEC(DEACTIVATED_SUB_CARRIERS) = []; 
         
        % Contorl on Pilots: FIXME, nof very accurate for more Pilots > 32
        SEL_PILOT_VECTOR = [1:ceil((CARRIER_NUM - length(DEACTIVATED_SUB_CARRIERS))/PILOT_NUM):(CARRIER_NUM-length(DEACTIVATED_SUB_CARRIERS))];
        
        PILOT_VECTOR = sort(PROB_PILOT_VEC(SEL_PILOT_VECTOR));
        
     case 'CHAOTIC_PILOT_DESIGN_METHOD'
        % Discard the randomly deactivated Sub-Carriers to hold a Pilot tone:
        PROB_PILOT_VEC(DEACTIVATED_SUB_CARRIERS) = [];
        
        % In this method all remainning subcarriers will have a duality,
        % namely they will be pilots and data combined!
        PILOT_VECTOR = PROB_PILOT_VEC;
     
     case 'FORD_CIRCLES'
         
        % Discard the randomly deactivated Sub-Carriers to hold a Pilot tone:
        PROB_PILOT_VEC(DEACTIVATED_SUB_CARRIERS) = [];
        
        % Step 2: Initialize set P = { p_zero, p_one}. Set n = 2:
        a(1)   = PROB_PILOT_VEC(1);
        ford_circle_factor = 1;
        
        %Forward Circles
        for u=2:(PILOT_NUM/2)
            
            if( selected_indx(u-1) <= length(PROB_PILOT_VEC))
                % Next position
                next_index         = 2*(ford_circle_factor)^2 + selected_indx(u-1);
                selected_indx(u)   = next_index;
                ford_circle_factor = ford_circle_factor + 1;
            else
                selected_indx(u)   = length(PROB_PILOT_VEC);
                break;
            end
            
        end
        
        %Backward Circles
        ford_circle_factor = 1;
        for u=((PILOT_NUM/2))+1:PILOT_NUM
            
            if( selected_indx(u-1) <= length(PROB_PILOT_VEC))
                % Next position
                next_index         = selected_indx(u-1) - 2*(ford_circle_factor)^2;
                selected_indx(u)   = next_index;
                ford_circle_factor = ford_circle_factor + 1;
            end
            
            if(selected_indx(u) <= 1)
                selected_indx(u) = 1;
                break;
            end
            
        end
        
        % In this method all remainning subcarriers will have a duality,
        % namely they will be pilots and data combined!
        PILOT_VECTOR = sort(PROB_PILOT_VEC(selected_indx));
        
     case 'ASSOCIATIVE_PILOT_ASSIGMENT_WITH_FEEDBACK_CHANNEL_INFORMATION'
        
        
        
        % 1) Check if any proposed item is deactivated already:
         save_repeated_index=[];
         for jjj=1:length(proposed_active_pilot_spaces)
             for deactivated_index=1:length(DEACTIVATED_SUB_CARRIERS)
                 if( proposed_active_pilot_spaces(jjj) == DEACTIVATED_SUB_CARRIERS(deactivated_index))
                     save_repeated_index = [save_repeated_index, jjj];
                 end
             end
         end
         proposed_active_pilot_spaces(save_repeated_index) =[];
         
         % 2) Discard the randomly deactivated Sub-Carriers to hold a Pilot tone:
        PROB_PILOT_VEC([DEACTIVATED_SUB_CARRIERS proposed_active_pilot_spaces]) = [];
        
        % 3) Build sub-carrier sections with variable widths and agrupate them:
        SubCarrier_Groups = {};
        carrier_group_indx      = 1;
        tmp_grp_elmnts_indx     = 1;
        non_priority_items_indx = 1;
        low_priority_subcarrier_indexes = PROB_PILOT_VEC;
        
        
        for u=1:length(PROB_PILOT_VEC)
            
            if (proposed_active_pilot_spaces)
            
                if ( tmp_grp_elmnts_indx == 1)
                    % Form the groups with higher priority:
                    temporal_group_elements(tmp_grp_elmnts_indx) = proposed_active_pilot_spaces(1);
                    priority  = 'high';
                    tmp_grp_elmnts_indx = tmp_grp_elmnts_indx + 1;
                else
                    % Consecutive elements???
                    if ((temporal_group_elements(tmp_grp_elmnts_indx-1) + 1) == proposed_active_pilot_spaces(1))
                        temporal_group_elements(tmp_grp_elmnts_indx) = proposed_active_pilot_spaces(1);
                        priority  = 'high';
                        tmp_grp_elmnts_indx = tmp_grp_elmnts_indx + 1;
                    else
                        %When they do not become consecutive elements load the
                        %group and start the new one:
                        % a) Load Previous Data:
                        SubCarrier_Groups{carrier_group_indx}.data     = temporal_group_elements;
                        SubCarrier_Groups{carrier_group_indx}.priority = priority;
                        carrier_group_indx = carrier_group_indx + 1;
                        % b) Start New Group:
                        tmp_grp_elmnts_indx     = 1;
                        temporal_group_elements = [];
                        temporal_group_elements(tmp_grp_elmnts_indx) = proposed_active_pilot_spaces(1);                   
                    end
                end

                proposed_active_pilot_spaces(1) = [];  
            else
                if (low_priority_subcarrier_indexes )
                    if (non_priority_items_indx == 1)
                        % Start New Group:
                        temporal_group_elements = [];                
                        temporal_group_elements(non_priority_items_indx) = low_priority_subcarrier_indexes(1);
                        priority  = 'low';
                        non_priority_items_indx = non_priority_items_indx + 1;
                    else
                        % Consecutive elements???
                        if ((temporal_group_elements(non_priority_items_indx-1) + 1) == low_priority_subcarrier_indexes(1))
                            temporal_group_elements(non_priority_items_indx) = low_priority_subcarrier_indexes(1);
                            priority  = 'low';
                            non_priority_items_indx = non_priority_items_indx + 1;
                        else
                            %When they do not become consecutive elements load the
                            %group and start the new one:
                            % a) Load Previous Data:
                            SubCarrier_Groups{carrier_group_indx}.data     = temporal_group_elements;
                            SubCarrier_Groups{carrier_group_indx}.priority = priority;
                            carrier_group_indx = carrier_group_indx + 1;
                            % b) Start New Group:
                            non_priority_items_indx     = 1;
                            temporal_group_elements     = [];
                            temporal_group_elements(non_priority_items_indx) = low_priority_subcarrier_indexes(1);                   
                        end
                    end
                    low_priority_subcarrier_indexes(1) = [];
                end
            end
            
        end
        
        % 4) Pilot Assigment:
        %    a) Know the relationship between group's number and number of
        %       Pilots configured:
        group_number_vs_pilot_number = PILOT_NUM / length(SubCarrier_Groups);
        
        if (group_number_vs_pilot_number == 1)
            % Distribute the pilots accordingly to the (tan) function, so
            % that the sine is the location of the current pilot, while the
            % cosine is the pilot from the previous signal.
            
            % Set the first element of the first group:
            pre_pilot_vector    = zeros(1,PILOT_NUM);
            
            for n=1:PILOT_NUM
                group_length        = length(SubCarrier_Groups{n}.data);
                
                % Try this idea later... :
                %selected_index = sin(SubCarrier_Groups{n}.data(1)) / cos(pre_pilot_vector(n-1));
                
                selected_index      = ceil(group_length/2); %Look for the middle point
                pre_pilot_vector(n) = SubCarrier_Groups{n}.data(selected_index);
                
            end
                      
        elseif(group_number_vs_pilot_number < 1)
            % We have less pilots than the number of groups: how to assign
            % the pilots?
            % First take the groups with only one member
            % Then look for the groups that has bigger number of members
            % and order them taking into account its priority
            %    Priority --> number of members
            pilot_index  = 0;
            pre_pilot_vector = zeros(1,PILOT_NUM);
            for group_index = 1:length(SubCarrier_Groups)
                group_length = length(SubCarrier_Groups{group_index}.data);
                
                if( group_length > 1 && ... 
                    group_length <= CHANNEL_LENGTH)
                    
                    pilot_index = pilot_index + 1;
                    pre_pilot_vector(pilot_index) = SubCarrier_Groups{group_index}.data(1);
                    
                end
                
                if (pilot_index == PILOT_NUM)
                    break;
                end
            end
            
            % Attend groups with more than 1 group memeber(if required), high priority first:
            if(pilot_index < PILOT_NUM)                
                for group_indx_w_more_members = 1:length(SubCarrier_Groups)
                    group_length = length(SubCarrier_Groups{group_indx_w_more_members}.data);
                    
                    if ( strcmp(SubCarrier_Groups{group_indx_w_more_members}.priority, 'high') && ...
                         group_length > CHANNEL_LENGTH )
                        
                        selected_index = ceil(group_length/2); %Look for the middle point
                        pilot_index    = pilot_index + 1;
                        pre_pilot_vector(pilot_index) = SubCarrier_Groups{group_indx_w_more_members}.data(selected_index);                       
                        
                    end
                    
                    if (pilot_index == PILOT_NUM)
                        break;
                    end
                end                
            end
            
            
            % Attend groups with more than 1 group memeber(if required), low priority last:
            if(pilot_index < PILOT_NUM)                
                for group_indx_w_more_members = 1:length(SubCarrier_Groups)
                    group_length = length(SubCarrier_Groups{group_indx_w_more_members}.data);
                    
                    if ( strcmp(SubCarrier_Groups{group_indx_w_more_members}.priority, 'low') && ...
                         group_length > CHANNEL_LENGTH  )
                        
                        selected_index = ceil(group_length/2); %Look for the middle point
                        pilot_index    = pilot_index + 1;
                        pre_pilot_vector(pilot_index) = SubCarrier_Groups{group_indx_w_more_members}.data(selected_index);
                                                
                    end
                    
                    if (pilot_index == PILOT_NUM)
                        break;
                    end
                end                
            end

        else  % (group_number_vs_pilot_number > 1)
            % We have more pilots than the number of groups: 

            % First asssing priority to the groups, meaning the number of
            % pilots are proportional to the number of members (on group's
            % length).
            % Then take the groups with only one member
            % Then look for the groups that has bigger number of members
            % and order them taking into account its priority
            %    Priority --> number of members            
            for group_index = 1:length(SubCarrier_Groups)
                SubCarrier_Groups{group_index}.length  = length(SubCarrier_Groups{group_index}.data);
                SubCarrier_Groups{group_index}.isTaken = 0;
            end
            
            % Make a list of number of pilots assigned per group:
            for group_index = 1:length(SubCarrier_Groups)-1
                if ( SubCarrier_Groups{group_index}.length > SubCarrier_Groups{group_index+1}.length )
                    max_index = group_index;
                    max_value = SubCarrier_Groups{group_index}.length;
                else
                    max_index = group_index + 1;
                    max_value = SubCarrier_Groups{group_index + 1}.length;
                end
            end
            
            % Calculate how many pilots to each group will be given:
            number_of_pilots_taken = 0;
            
            % Initialize the number of pilots assigned:
            for group_index = 1:length(SubCarrier_Groups)                
                SubCarrier_Groups{group_index}.number_of_pilots = 0;              
            end
            
            if (number_of_pilots_taken < PILOT_NUM)   
                %symetric_spacing = ceil( (CARRIER_NUM - length(DEACTIVATED_SUB_CARRIERS))/(PILOT_NUM - number_of_pilots_taken)); 
                symetric_spacing = CHANNEL_LENGTH;
                for group_index = 1:length(SubCarrier_Groups)
                    if( SubCarrier_Groups{group_index}.length < symetric_spacing && ...
                        SubCarrier_Groups{group_index}.length > 1                      )  % Ignore groups with only one member!
                        SubCarrier_Groups{group_index}.number_of_pilots = 1;
                        number_of_pilots_taken = number_of_pilots_taken + 1;
                    end
                    
                    if (number_of_pilots_taken == PILOT_NUM)
                        break;
                    end
                end
                
            end
            
            
            if (number_of_pilots_taken < PILOT_NUM)
            
                for group_index = 1:length(SubCarrier_Groups)
                    if((SubCarrier_Groups{group_index}.length >= symetric_spacing) && ...
                        strcmp(SubCarrier_Groups{group_index}.priority, 'high'))
                        % Calculate the number of pilots:
                        numberOfPilots_in_group = ceil(SubCarrier_Groups{group_index}.length / (symetric_spacing));
                        next_number_of_Pilots_taken = numberOfPilots_in_group + number_of_pilots_taken;
                        if ( next_number_of_Pilots_taken == PILOT_NUM)
                            SubCarrier_Groups{group_index}.number_of_pilots = numberOfPilots_in_group;
                            number_of_pilots_taken = number_of_pilots_taken + numberOfPilots_in_group;
                        elseif (next_number_of_Pilots_taken > PILOT_NUM)
                            % Just add the remaining number of pilots:
                            SubCarrier_Groups{group_index}.number_of_pilots = numberOfPilots_in_group - (next_number_of_Pilots_taken - PILOT_NUM);
                            number_of_pilots_taken = number_of_pilots_taken + SubCarrier_Groups{group_index}.number_of_pilots;
                        else
                            SubCarrier_Groups{group_index}.number_of_pilots = numberOfPilots_in_group;
                            number_of_pilots_taken = number_of_pilots_taken + numberOfPilots_in_group;
                        end                        
                    end
                    
                    if (number_of_pilots_taken == PILOT_NUM)
                        break;
                    end
                end
                
            end
            
            if (number_of_pilots_taken < PILOT_NUM)
            
                for group_index = 1:length(SubCarrier_Groups)
                    if((SubCarrier_Groups{group_index}.length >= symetric_spacing) && ...
                        strcmp(SubCarrier_Groups{group_index}.priority, 'low'))
                        % Calculate the number of pilots:
                        numberOfPilots_in_group = ceil(SubCarrier_Groups{group_index}.length / symetric_spacing);
                        next_number_of_Pilots_taken = numberOfPilots_in_group + number_of_pilots_taken;
                        if ( next_number_of_Pilots_taken == PILOT_NUM)
                            SubCarrier_Groups{group_index}.number_of_pilots = numberOfPilots_in_group;
                            number_of_pilots_taken = number_of_pilots_taken + numberOfPilots_in_group;
                        elseif (next_number_of_Pilots_taken > PILOT_NUM)
                            % Just add the remaining number of pilots:
                            SubCarrier_Groups{group_index}.number_of_pilots = numberOfPilots_in_group - (next_number_of_Pilots_taken - PILOT_NUM);
                            number_of_pilots_taken = number_of_pilots_taken + SubCarrier_Groups{group_index}.number_of_pilots;
                        else
                            SubCarrier_Groups{group_index}.number_of_pilots = numberOfPilots_in_group;
                            number_of_pilots_taken = number_of_pilots_taken + numberOfPilots_in_group;
                        end                        
                    end
                    
                    if (number_of_pilots_taken == PILOT_NUM)
                        break;
                    end
                end
                
            end
                   
            % Final Step: Accomodate the pilots in the rigth place,
            % according to what they have in their specification.
            
            pilot_index  = 0;
            pre_pilot_vector = zeros(1,PILOT_NUM);
            for group_index = 1:length(SubCarrier_Groups)
                
                NumberOfPilotsInGroup = SubCarrier_Groups{group_index}.number_of_pilots;
                
                if( NumberOfPilotsInGroup == 1)
                    
                    pilot_index = pilot_index + 1;
                    pre_pilot_vector(pilot_index)          = SubCarrier_Groups{group_index}.data(1);
                    SubCarrier_Groups{group_index}.isTaken = 1;
                    
                elseif ( NumberOfPilotsInGroup > 1)
                    GroupSpacing  = floor(SubCarrier_Groups{group_index}.length / NumberOfPilotsInGroup);
                    
                    for group_assign_indx = 1:NumberOfPilotsInGroup
                        
                        pilot_index = pilot_index + 1;
                        pre_pilot_vector(pilot_index) = SubCarrier_Groups{group_index}.data(group_assign_indx + (GroupSpacing * (group_assign_indx -1)));
                        SubCarrier_Groups{group_index}.isTaken = 1;                        
                        
                    end
                    
                else
                    
                end
                
            end          
            
        end
        
        %Remove the zero elements if pre_pilot_vector size < PILOT_NUM:
        %Optimization Factor!!!
        save_zero_index = [];
        for rm_index = 1:length(pre_pilot_vector)
            if (pre_pilot_vector(rm_index) == 0)
                save_zero_index = [save_zero_index rm_index];
            end
        end
        pre_pilot_vector(save_zero_index)=[];
        
        PILOT_VECTOR = sort(pre_pilot_vector);
        
        
     case 'SYMETRIC_PILOTS_SET_AND_INFLUENCE_OF_HIGHEST_CH_VARIATIONS'
        
        % -----------------------------------------------------------------
        % 1) First symetric pilot placement:
        % -----------------------------------------------------------------
        % Discard the randomly deactivated Sub-Carriers to hold a Pilot tone:
        PROB_PILOT_VEC(DEACTIVATED_SUB_CARRIERS) = []; 
         
        % Contorl on Pilots: FIXME, nof very accurate for more Pilots > 32
        SEL_PILOT_VECTOR = [1:ceil((CARRIER_NUM - length(DEACTIVATED_SUB_CARRIERS))/PILOT_NUM):(CARRIER_NUM-length(DEACTIVATED_SUB_CARRIERS))];
        
        PILOT_VECTOR      = sort(PROB_PILOT_VEC(SEL_PILOT_VECTOR));        
        PILOT_VECTOR_COPY = PILOT_VECTOR;
        % -----------------------------------------------------------------
        % 2) Second Influence of the highest channel variations:
        % -----------------------------------------------------------------
        
        % Check if any proposed item is deactivated already:
        save_repeated_index=[];
         for jjj=1:length(proposed_active_pilot_spaces)
             for deactivated_index=1:length(DEACTIVATED_SUB_CARRIERS)
                 if( proposed_active_pilot_spaces(jjj) == DEACTIVATED_SUB_CARRIERS(deactivated_index))
                     save_repeated_index = [save_repeated_index, jjj];
                 end
             end
         end
        proposed_active_pilot_spaces(save_repeated_index) =[];
        
        % a) Select a window sizeto check if more proposed pilot positions
        % are present
        window_size = CARRIER_NUM/length(PILOT_VECTOR_COPY);
        
        for pilot_indx = 1:length(PILOT_VECTOR_COPY)
            
            [value, index] = min(abs(proposed_active_pilot_spaces - PILOT_VECTOR_COPY(pilot_indx)));
            
            PILOT_VECTOR_COPY(pilot_indx) = proposed_active_pilot_spaces(index);
            proposed_active_pilot_spaces(index) = [];
        end
        
        PILOT_VECTOR = PILOT_VECTOR_COPY;
        
     case 'PILOTS_ONLY_ON_THE_HIGHEST_CH_VARIATIONS'
         
         % Discard the randomly deactivated Sub-Carriers to hold a Pilot tone:
         PROB_PILOT_VEC(DEACTIVATED_SUB_CARRIERS) = [];
         
         %Proposed Algorithm steps:
         
         % Check if any proposed item is deactivated already:
         save_repeated_index=[];
         for jjj=1:length(proposed_active_pilot_spaces)
             for deactivated_index=1:length(DEACTIVATED_SUB_CARRIERS)
                 if( proposed_active_pilot_spaces(jjj) == DEACTIVATED_SUB_CARRIERS(deactivated_index))
                     save_repeated_index = [save_repeated_index, jjj];
                 end
             end
         end
         proposed_active_pilot_spaces(save_repeated_index) =[];
         
         % 1) Select first and last element from the remaining set
         selected_indx(1)   = PROB_PILOT_VEC(1);
         selected_indx(end) = PROB_PILOT_VEC(end);
         
         % 2) Remove Selected Items from the set:
         %PROB_PILOT_VEC([1,end]) = [];
               
         %Checking order is increasing...
         proposed_active_pilot_spaces = sort(proposed_active_pilot_spaces);
         %checking if the first and last are still present to avoid
         %repeated items:
         save_repeated_index=[];
         for j=1:length(proposed_active_pilot_spaces)
             if( (proposed_active_pilot_spaces(j) == selected_indx(1)) || ...
                 (proposed_active_pilot_spaces(j) == selected_indx(end)) )
                 save_repeated_index = [save_repeated_index, j];
             end       
         end         
         proposed_active_pilot_spaces(save_repeated_index) = [];
         
         % 3) From the remaining proposed, set equally distribute the pilots:
         space_between_pilots = floor( length(proposed_active_pilot_spaces) / (PILOT_NUM-2));
         for i=1:PILOT_NUM-2
             
             if (i == 1)
                 selected_indx(i+1) = proposed_active_pilot_spaces(1);
             else
                 selected_indx(i+1) = proposed_active_pilot_spaces( (space_between_pilots*(i-1)) - 1 );
             end
             
         end         
         
         PILOT_VECTOR = selected_indx;
         
     case 'FULL_WAVELET_PILOT_DESIGN_METHOD'
         % 1) Discard the randomly deactivated Sub-Carriers to hold a Pilot tone:
        %PROB_PILOT_VEC(DEACTIVATED_SUB_CARRIERS) = [];        
        
        beta_coef_lim = 0.1600; %min(min(proposed_active_pilot_spaces)); %0.0;  % 0.1600
        
        % 2) Calculate the number of Wavelet Scales based on the highest
        % frequency present Sca, row. Only 25% above of the highest wavelet scale
        % to be accounted.
        
        number_of_wavelet_scales = size(proposed_active_pilot_spaces,1);
               
        for wavelet_indx_scale = 1:number_of_wavelet_scales
            
            wavelet_scale_coeficients = proposed_active_pilot_spaces(wavelet_indx_scale, :);
            
            % Apply Min & MAX Approach on search of values:
            max_min_coefficient_indexes = [];
            
            for coeficient_indx = 1:length(wavelet_scale_coeficients)
    
                % Is the only one maximum?
                %if ((wavelet_scale_coeficients(coeficient_indx) / (beta_coef_lim) >= wavelet_scale_coeficients(coeficient_indx)) && ...
                %    (wavelet_scale_coeficients(coeficient_indx) / ((beta_coef_lim) - ((10*beta_coef_lim)/100)) >= wavelet_scale_coeficients(coeficient_indx)) )
                %    max_min_coefficient_indexes = [max_min_coefficient_indexes coeficient_indx];
                %end
                    
                % Is the only one maximum?
                if ((wavelet_scale_coeficients(coeficient_indx) / (beta_coef_lim) <= wavelet_scale_coeficients(coeficient_indx)))
                    max_min_coefficient_indexes = [max_min_coefficient_indexes coeficient_indx];
                end            
%                 % Is the only one minimum?
%                 if (wavelet_scale_coeficients(coeficient_indx) / min_coefficient_val < 0.01)
%                     max_min_coefficient_indexes = [max_min_coefficient_indexes coeficient_indx];
%                 end
                
            end
            
            selected_coefficients_from_all_scales{wavelet_indx_scale} = max_min_coefficient_indexes;
                                 
        end
        
        
%         Now desition has to be made of all those components which ones take,
%         so the rule could be:
%         a) The indexes that are repeated the more, are choosen
%         b) should I weigth some scales more??? which ones? how the weigth
%         should be choosen??
        SELECTED_PILOT_VECTOR = [];
        for wavelet_indx_scale = 1:number_of_wavelet_scales
            
            SELECTED_PILOT_VECTOR = [SELECTED_PILOT_VECTOR selected_coefficients_from_all_scales{wavelet_indx_scale}];
            
        end
        
        % removes redundancy from items that were repeadedly selected across several scales.
        % but no weithing is employed to aply this reduction, just
        % repeatition.
        merged_selected_positions_indexes = unique(SELECTED_PILOT_VECTOR);
        
        
        % Remove the randomly deactivated subcarriers:
        items_to_remove = [];
        for deactivated_carrier_num_indx = 1:length(DEACTIVATED_SUB_CARRIERS)
            
            for carrier_num_indx = 1:length(merged_selected_positions_indexes)
            
                if(DEACTIVATED_SUB_CARRIERS(deactivated_carrier_num_indx) == merged_selected_positions_indexes(carrier_num_indx))
                    items_to_remove = [items_to_remove  carrier_num_indx];
                end
            end
            
        end
        
        merged_selected_positions_indexes(items_to_remove) = [];
        
        
        if(length(merged_selected_positions_indexes) > PILOT_NUM)
            
            %Calculate the saparation space:
            spacing = ceil(length(merged_selected_positions_indexes) / PILOT_NUM); 
            
            selected_positions_indexes = merged_selected_positions_indexes(1:spacing:end);
            
        elseif(length(merged_selected_positions_indexes) == PILOT_NUM)
            selected_positions_indexes = merged_selected_positions_indexes;
        else
            selected_positions_indexes = merged_selected_positions_indexes;
        end
        
        %Debug:
%         figure(3500)
%         stem(selected_positions_indexes, PROB_PILOT_VEC(selected_positions_indexes))
        
        % Last) Organize in order the proposed pilot positiions:
        PILOT_VECTOR = selected_positions_indexes;
              

%Alternative code:

% 
%         % 1) Discard the randomly deactivated Sub-Carriers to hold a Pilot tone:
%         PROB_PILOT_VEC(DEACTIVATED_SUB_CARRIERS) = [];        
%         
%         % 2) Calculate the number of Wavelet Scales based on the highest
%         % frequency present Sca, row. Only 25% above of the highest wavelet scale
%         % to be accounted.
%         
%         number_of_wavelet_scales = size(proposed_active_pilot_spaces,1);
%         max_coefficient_val      = roundn(max(max(proposed_active_pilot_spaces)), -4);
%         min_coefficient_val      = roundn((max_coefficient_val * 0.15), -4);
%         
%         % Apply Min & MAX Approach on search of values:
%         max_coefficient_scales = [];
%             
%         for wavelet_indx_scale = 1:number_of_wavelet_scales
%             
%             wavelet_scale_coeficients = proposed_active_pilot_spaces(wavelet_indx_scale, :);
%             
%             % Apply Min & MAX Approach on search of values:
%             max_min_coefficient_indexes = [];
%             
%             for coeficient_indx = 1:length(wavelet_scale_coeficients)
%     
%                 % Is the only one maximum?
%                 if ( wavelet_scale_coeficients(coeficient_indx) >= (max_coefficient_val - (max_coefficient_val * 0.10)) )
%                     max_min_coefficient_indexes = [max_min_coefficient_indexes coeficient_indx];
%                 end
%                     
%                             
%                 % Is the only one minimum?
%                 if (min_coefficient_val >= wavelet_scale_coeficients(coeficient_indx) <= (min_coefficient_val + (min_coefficient_val * 0.20)))
%                     max_min_coefficient_indexes = [max_min_coefficient_indexes coeficient_indx];
%                 end
%                 
%             end
%             
%             selected_coefficients_from_all_scales{wavelet_indx_scale} = max_min_coefficient_indexes;
% 
%             
%             
%             for coeficient_indx = 1:length(wavelet_scale_coeficients)
%     
%                 % Is this scale the one that has one of the maximum value?
%                 if ( roundn(wavelet_scale_coeficients(coeficient_indx),-4) == max_coefficient_val )
%                     max_coefficient_scales = [max_coefficient_scales wavelet_indx_scale];
%                     break;
%                 end
%                 
%             end
%                                  
%         end
%         
%         
%         for wavelet_indx_scale = 1:length(max_coefficient_scales)
%             
%             wavelet_scale_coeficients = proposed_active_pilot_spaces(max_coefficient_scales(wavelet_indx_scale), :);
%             
%             % Apply Min & MAX Approach on search of values:
%             max_min_coefficient_indexes = [];
%             
%             for coeficient_indx = 1:length(wavelet_scale_coeficients)
%     
%                 % Is the only one maximum?
%                 if ( roundn(wavelet_scale_coeficients(coeficient_indx),-4) >= roundn((max_coefficient_val - (max_coefficient_val * 0.20)),-4) )
%                     max_min_coefficient_indexes = [max_min_coefficient_indexes coeficient_indx];
%                 end
%                     
%                             
%                 % Is the only one minimum?
%                 if ((min_coefficient_val                                    <= roundn(wavelet_scale_coeficients(coeficient_indx),-4)) && ...
%                     ( roundn(wavelet_scale_coeficients(coeficient_indx),-4) <= (roundn((min_coefficient_val + (min_coefficient_val * 0.25)), -4))))
%                     max_min_coefficient_indexes = [max_min_coefficient_indexes coeficient_indx];
%                 end
%                 
%             end
%             
%             selected_coefficients_from_all_scales{wavelet_indx_scale} = max_min_coefficient_indexes;
% 
%                             
%         end
%         
%         
%         
%         
%         
%         
%         Now desition has to be made of all those components which ones take,
%         so the rule could be:
%         a) The indexes that are repeated the more, are choosen
%         b) should I weigth some scales more??? which ones? how the weigth
%         should be choosen??
%         SELECTED_PILOT_VECTOR = [];
%         %for wavelet_indx_scale = 1:number_of_wavelet_scales
%         for wavelet_indx_scale = 1:length(max_coefficient_scales)   
%             SELECTED_PILOT_VECTOR = [SELECTED_PILOT_VECTOR selected_coefficients_from_all_scales{wavelet_indx_scale} ];
%             
%         end
%         
%         % removes redundancy from items that were repeadedly selected across several scales.
%         % but no weithing is employed to aply this reduction, just
%         % repeatition.
%         merged_selected_positions_indexes = unique(SELECTED_PILOT_VECTOR);
%         
%         % Remove the randomly deactivated subcarriers:
%         items_to_remove = [];
%         for deactivated_carrier_num_indx = 1:length(DEACTIVATED_SUB_CARRIERS)
%             
%             for carrier_num_indx = 1:length(merged_selected_positions_indexes)
%             
%                 if(DEACTIVATED_SUB_CARRIERS(deactivated_carrier_num_indx) == merged_selected_positions_indexes(carrier_num_indx))
%                     items_to_remove = [items_to_remove  carrier_num_indx];
%                 end
%             end
%             
%         end
%         
%         merged_selected_positions_indexes(items_to_remove) = [];
%         
%         if(length(merged_selected_positions_indexes) > PILOT_NUM)
%             
%             %Calculate the saparation space:
%             %First Approach:
%             %spacing = ceil(length(merged_selected_positions_indexes) / PILOT_NUM); 
%             %selected_positions_indexes = merged_selected_positions_indexes(1:spacing:end);
% 
%             %Second Approach:
%             spacing = ceil(length(merged_selected_positions_indexes) / PILOT_NUM); 
%             selected_positions_indexes = merged_selected_positions_indexes(1:spacing:end);
%             
%         elseif(length(merged_selected_positions_indexes) == PILOT_NUM)
%             selected_positions_indexes = merged_selected_positions_indexes;
%         else
%             selected_positions_indexes = merged_selected_positions_indexes;
%         end
%         
%         % Last) Organize in order the proposed pilot positiions:
%         PILOT_VECTOR = selected_positions_indexes;

     case 'WAVELET_ENERGY_BASED'
         
         % 1) Discard the randomly deactivated Sub-Carriers to hold a Pilot tone:
        PROB_PILOT_VEC(DEACTIVATED_SUB_CARRIERS) = [];
        
        wavelet_coeficients     = proposed_active_pilot_spaces;
        wavelet_step_energy_sum = sum(wavelet_coeficients);

        % 2) Algorithm:
        % Find Max and Min values as per derivative == 0
        for ii = 1:length(wavelet_step_energy_sum)
            if (wavelet_step_energy_sum(ii) < 0 )
              xi_min_sig(ii) = -1 * (wavelet_step_energy_sum(ii));
            else
              xi_min_sig(ii) = wavelet_step_energy_sum(ii);
            end
        end
        [~,selected_set]     = findpeaks(xi_min_sig);

        
%         %Debug:
%         figure(3500)
%         stem(wavelet_step_energy_sum)
        
        % 3) Remove the Randomly Deactivated Subcarriers if any choosen from the algorithm:
        items_to_remove = [];
        for deactivated_carrier_num_indx = 1:length(DEACTIVATED_SUB_CARRIERS)
            
            for carrier_num_indx = 1:length(selected_set)      
                if(DEACTIVATED_SUB_CARRIERS(deactivated_carrier_num_indx) == selected_set(carrier_num_indx))
                    items_to_remove = [items_to_remove  carrier_num_indx];
                end
            end
            
        end
        
        selected_set(items_to_remove) = [];
       
        selected_set_length = length(selected_set);
        
        % 4) Warranty at least two items for interpolation:
        if(selected_set_length > PILOT_NUM)
            
            %Calculate the saparation space:
            spacing                    = ceil(length(selected_set) / PILOT_NUM);      
            selected_positions_indexes = selected_set(1:spacing:end);
            
        elseif(selected_set_length == PILOT_NUM)
            
            selected_positions_indexes = selected_set;
            
        elseif(selected_set_length <= 1)
            
            % Add first and last subcarrier positions:
            selected_set = [ PROB_PILOT_VEC(1) selected_set  PROB_PILOT_VEC(end) ];
            selected_positions_indexes = selected_set;
            
        elseif( (selected_set_length >  1)         && ...
                (selected_set_length <  PILOT_NUM) && ...
                (FORCE_PILOT_NUM     == 1))
            
            xi_min_si_temp               = xi_min_sig;
            xi_min_si_temp(selected_set) = 0; % No min nor max of already selected items.
            
            while(selected_set_length ~= PILOT_NUM)
                
                [~,selected_set_temp] = findpeaks(xi_min_si_temp);
                
                %Remove the Deactivated Subcarriers if any choosen from the algorithm:
                items_to_remove = [];
                for deactivated_carrier_num_indx = 1:length(DEACTIVATED_SUB_CARRIERS)
            
                    for carrier_num_indx = 1:length(selected_set_temp)      
                        if(DEACTIVATED_SUB_CARRIERS(deactivated_carrier_num_indx) == selected_set_temp(carrier_num_indx))
                            items_to_remove = [items_to_remove  carrier_num_indx];
                        end
                    end
            
                end
                
                
                if(length(selected_set_temp) == length(items_to_remove))
                    % No min nor max of already selected items:
                    xi_min_si_temp(selected_set_temp) = 0; 
                else
                
                    selected_set_temp(items_to_remove) = [];
                    
                    if( (length(selected_set_temp) + selected_set_length) > PILOT_NUM )
                        % Calculate how many additional Pilots are required:
                        N_required_additional_p = PILOT_NUM - selected_set_length;
                        selected_set_temp       = selected_set_temp(1:N_required_additional_p);
                    end

                    selected_set = [selected_set selected_set_temp];

                    selected_set_length = length(selected_set);
                    % No min nor max of already selected items:
                    xi_min_si_temp(selected_set_temp) = 0;
                end
            end
            
            selected_positions_indexes = selected_set;
            
        else
            selected_positions_indexes = selected_set;
        end
        
%         %Debug:
%         figure(3501)
%         stem(selected_set, wavelet_step_energy_sum(selected_set))
        
        % Remove redundant items:
        selected_positions_indexes = unique(selected_positions_indexes);
        
        % Last) Organize in order the proposed pilot positiions:
        PILOT_VECTOR = selected_positions_indexes;
         
     otherwise
        error('Not a correct pilot design method supported! (see ''PilotDesign.METHOD'')');    
 end

end

