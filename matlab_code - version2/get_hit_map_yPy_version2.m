function Estruct = get_hit_map_yPy_version2(coords_to_use, parameters)
    
    
    
    w1 = parameters.w1;
    w2 = parameters.w2;
    Et = parameters.Et;
    Ex = parameters.Ex;
    Ey = Et - Ex;
    limit_number_of_round_hits = parameters.limit_number_of_round_hits; %boolian for limiting number of hits
    number_of_round_hits_treshold = parameters.number_of_round_hits_treshold;
    w0 = [w1 w2];
    initial_conditions_preHit = zeros(length(coords_to_use)*4, 4);
    initial_conditions_postHit = zeros(length(coords_to_use)*4, 4);
    results_preHit = zeros(length(coords_to_use)*4, 4);
    results_postHit = zeros(length(coords_to_use)*4, 4);
    results_NumberOfHits = zeros(length(coords_to_use)*4, 3); % (horizontal, circle, vertical)
    final_E_preHit = zeros(length(coords_to_use)*4, 1);
    final_E_postHit = zeros(length(coords_to_use)*4, 1);
    idx = 1;


    for w_idc = 1:length(coords_to_use)
        x_w = coords_to_use(w_idc, 1);
        y_w = coords_to_use(w_idc, 2);
        Ex_min = 0.5*(w1^2)*(x_w^2);
        EY_min = 0.5*(w2^2)*(y_w^2);
        

%         if (idx>=222) && (idx<=230)
%             disp(idx)
%         end

        if (Ex<Ex_min) || (Ey<EY_min)
            continue
        end
    %     energies = linspace(Ex_min, Et-EY_min, numberOfenergies);
        
        % the moments:
        pxw_abs = sqrt(2*Ex - (w1*x_w)^2);
        pyw_abs = sqrt(2*Ey - (w2*y_w)^2);
        
        % 4 possible directions:
        q0_right_up = [x_w; y_w; pxw_abs; pyw_abs];
        q0_right_down = [x_w; y_w; pxw_abs; -pyw_abs];
        q0_left_up = [x_w; y_w; -pxw_abs; pyw_abs];
        q0_left_down = [x_w; y_w; -pxw_abs; -pyw_abs];
        
        flag_track_right_up = checkIfInsideEuler(q0_right_up);
        flag_track_right_down = checkIfInsideEuler(q0_right_down);
        flag_track_left_up = checkIfInsideEuler(q0_left_up);
        flag_track_left_down = checkIfInsideEuler(q0_left_down);
    
        % 4 solvers, some of them will stop fast
    
        
    
    
    %     figure;
    %     hold on
    %     plot(coords_to_use(:,1),coords_to_use(:,2),'ro')
    %     plot(q0_track_right_up(:,1),q0_track_right_up(:,2),'bo')
    %     plot(q0_track_right_down(:,1),q0_track_right_down(:,2),'go')
    %     plot(q0_track_left_up(:,1),q0_track_left_up(:,2),'ko')
    %     plot(q0_track_left_down(:,1),q0_track_left_down(:,2),'mo')
    %     hold off
    
        % check if its a valid track for each of them, and add to final vector
        
    
        if ~flag_track_right_up
            [q0_track_right_up, q0_hits_round_right_up, ...
                q0_hits_horizontal_right_up, q0_hits_vertical_right_up] ...
                = calc_track_until_px_is_zero(q0_right_up, w0, 0);
            final_PreHit_right_up = q0_track_right_up(end, :);
            p_postHit_right_up = return_momentum(-[q0_right_up(3) q0_right_up(4)], coords_to_use(w_idc, 3));
            q1_right_up = [x_w; y_w; p_postHit_right_up(1); p_postHit_right_up(2)];
            [q1_track_right_up, q1_hits_round_right_up, ...
                q1_hits_horizontal_right_up, q1_hits_vertical_right_up] ...
                = calc_track_until_px_is_zero(q1_right_up, w0, 1);
            final_PostHit_right_up = q1_track_right_up(end, :);
            
%             figure;
%             hold on
%             plot(coords_to_use(:,1),coords_to_use(:,2),'ro')
% %             plot(coords_circ(:,1),coords_circ(:,2),'ro')
%             plot(q0_track_right_up(:,1),q0_track_right_up(:,2),'bo')
%             plot(q1_track_right_up(:,1),q1_track_right_up(:,2),'go')
%             legend('step', 'step', 'preHit', 'postHit')
%             xlim([-10,10])
%             ylim([-10,10])
%             hold off

            total_number_of_round_hits_rightUp = ...
                        q0_hits_round_right_up + q1_hits_round_right_up;
            total_number_of_horizontal_hits_rightUp = ...
                        q0_hits_horizontal_right_up + q1_hits_horizontal_right_up;
            total_number_of_vertical_hits_rightUp = ...
                        q0_hits_vertical_right_up + q1_hits_vertical_right_up;

            if (~limit_number_of_round_hits) || ...
            (total_number_of_round_hits_rightUp <= number_of_round_hits_treshold)
                % locate results
                final_E_preHit(idx) = calcTotalEnergy(final_PreHit_right_up);
                final_E_postHit(idx) = calcTotalEnergy(final_PostHit_right_up);
                results_preHit(idx, :) = final_PreHit_right_up;
                results_postHit(idx, :) = final_PostHit_right_up;
                initial_conditions_preHit(idx, :) = q0_right_up.';
                initial_conditions_postHit(idx, :) = q1_right_up.';
                results_NumberOfHits(idx, :) = ...
                    [total_number_of_horizontal_hits_rightUp,...
                     total_number_of_round_hits_rightUp,...
                     total_number_of_vertical_hits_rightUp];
                idx = idx + 1;
            end
        end
        
        
        if ~flag_track_right_down
            [q0_track_right_down, q0_hits_round_right_down, ...
                q0_hits_horizontal_right_down, q0_hits_vertical_right_down] ...
                = calc_track_until_px_is_zero(q0_right_down, w0, 0);
            final_PreHit_right_down = q0_track_right_down(end, :);

%             if idx==366
%                 disp(idx)
%             end
            p_postHit_right_down = return_momentum(-[q0_right_down(3) q0_right_down(4)], coords_to_use(w_idc, 3));
            q1_right_down = [x_w; y_w; p_postHit_right_down(1); p_postHit_right_down(2)];
    %         disp(idx)
            [q1_track_right_down, q1_hits_round_right_down, ...
                q1_hits_horizontal_right_down, q1_hits_vertical_right_down] = ...
                calc_track_until_px_is_zero(q1_right_down, w0, 1);
            
            final_PostHit_right_down = q1_track_right_down(end, :);
            
%             figure;
%             hold on
%             plot(coords_to_use(:,1),coords_to_use(:,2),'ro')
%             plot(q0_track_right_down(:,1),q0_track_right_down(:,2),'bo')
%             plot(q1_track_right_down(:,1),q1_track_right_down(:,2),'go')
%             hold off
            
            % locate results
            total_number_of_round_hits_rightDown = ...
                        q0_hits_round_right_down + q1_hits_round_right_down;
            total_number_of_horizontal_hits_rightDown = ...
                        q0_hits_horizontal_right_down + q1_hits_horizontal_right_down;
            total_number_of_vertical_hits_rightDown = ...
                        q0_hits_vertical_right_down + q1_hits_vertical_right_down;

            if (~limit_number_of_round_hits) || ...
            (total_number_of_round_hits_rightDown <= number_of_round_hits_treshold)

                final_E_preHit(idx) = calcTotalEnergy(final_PreHit_right_down);
                final_E_postHit(idx) = calcTotalEnergy(final_PostHit_right_down);
                results_preHit(idx, :) = final_PreHit_right_down;
                results_postHit(idx, :) = final_PostHit_right_down;
                initial_conditions_preHit(idx, :) = q0_right_down.';
                initial_conditions_postHit(idx, :) = q1_right_down.';
                results_NumberOfHits(idx, :) = ...
                    [total_number_of_horizontal_hits_rightDown,...
                     total_number_of_round_hits_rightDown,...
                     total_number_of_vertical_hits_rightDown];
                idx = idx + 1;
            end
        end
    
    
        
        if ~flag_track_left_up
            [q0_track_left_up, q0_hits_round_left_up, ...
                q0_hits_horizontal_left_up, q0_hits_vertical_left_up] ...
                = calc_track_until_px_is_zero(q0_left_up, w0, 0);
            final_PreHit_left_up = q0_track_left_up(end, :);
            p_postHit_left_up = return_momentum(-[q0_left_up(3) q0_left_up(4)], coords_to_use(w_idc, 3));
            q1_left_up = [x_w; y_w; p_postHit_left_up(1); p_postHit_left_up(2)];
            [q1_track_left_up, q1_hits_round_left_up, ...
                q1_hits_horizontal_left_up, q1_hits_vertical_left_up] ...
                = calc_track_until_px_is_zero(q1_left_up, w0, 1);
            final_PostHit_left_up = q1_track_left_up(end, :);

%             figure;
%             hold on
%             plot(coords_to_use(:,1),coords_to_use(:,2),'ro')
%             plot(coords_circ(:,1),coords_circ(:,2),'ro')
%             plot(q0_track_left_up(:,1),q0_track_left_up(:,2),'bo')
%             plot(q1_track_left_up(:,1),q1_track_left_up(:,2),'go')
%             legend('step', 'step', 'preHit', 'postHit')
%             xlim([-4,4])
%             ylim([-4,4])
%             hold off

            % locate results
            total_number_of_round_hits_leftUp = ...
                        q0_hits_round_left_up + q1_hits_round_left_up;
            total_number_of_horizontal_hits_leftUp = ...
                        q0_hits_horizontal_left_up + q1_hits_horizontal_left_up;
            total_number_of_vertical_hits_leftUp = ...
                        q0_hits_vertical_left_up + q1_hits_vertical_left_up;

            if (~limit_number_of_round_hits) || ...
            (total_number_of_round_hits_leftUp <= number_of_round_hits_treshold)

                final_E_preHit(idx) = calcTotalEnergy(final_PreHit_left_up);
                final_E_postHit(idx) = calcTotalEnergy(final_PostHit_left_up);
                results_preHit(idx, :) = final_PreHit_left_up;
                results_postHit(idx, :) = final_PostHit_left_up;
                initial_conditions_preHit(idx, :) = q0_left_up.';
                initial_conditions_postHit(idx, :) = q1_left_up.';
                results_NumberOfHits(idx, :) = ...
                    [total_number_of_horizontal_hits_leftUp,...
                     total_number_of_round_hits_leftUp,...
                     total_number_of_vertical_hits_leftUp];
                idx = idx + 1;
            end
        end
        
        
        if ~flag_track_left_down
    
            [q0_track_left_down, q0_hits_round_left_down, ...
                q0_hits_horizontal_left_down, q0_hits_vertical_left_down]...
                = calc_track_until_px_is_zero(q0_left_down, w0, 0);
            final_PreHit_left_down = q0_track_left_down(end, :);
            p_postHit_left_down = return_momentum(-[q0_left_down(3) q0_left_down(4)], coords_to_use(w_idc, 3));
            q1_left_down = [x_w; y_w; p_postHit_left_down(1); p_postHit_left_down(2)];
            [q1_track_left_down, q1_hits_round_left_down, ...
                q1_hits_horizontal_left_down, q1_hits_vertical_left_down]...
                = calc_track_until_px_is_zero(q1_left_down, w0, 1);
            
            final_PostHit_left_down = q1_track_left_down(end, :);
    
            % locate results
            total_number_of_round_hits_leftDown = ...
                        q0_hits_round_left_down + q1_hits_round_left_down;
            total_number_of_horizontal_hits_leftDown = ...
                        q0_hits_horizontal_left_down + q1_hits_horizontal_left_down;
            total_number_of_vertical_hits_leftDown = ...
                        q0_hits_vertical_left_down + q1_hits_vertical_left_down;

            if (~limit_number_of_round_hits) || ...
            (total_number_of_round_hits_leftDown <= number_of_round_hits_treshold)
                
                final_E_preHit(idx) = calcTotalEnergy(final_PreHit_left_down);
                final_E_postHit(idx) = calcTotalEnergy(final_PostHit_left_down);
                initial_conditions_preHit(idx, :) = q0_left_down.';
                initial_conditions_postHit(idx, :) = q1_left_down.';
                results_preHit(idx, :) = final_PreHit_left_down;
                results_postHit(idx, :) = final_PostHit_left_down;
                results_NumberOfHits(idx, :) = ...
                    [total_number_of_horizontal_hits_leftDown,...
                     total_number_of_round_hits_leftDown,...
                     total_number_of_vertical_hits_leftDown];
                idx = idx + 1;
            end
        end

    end

    results_preHit = results_preHit(1:idx-1, :);
    results_postHit = results_postHit(1:idx-1, :);
    final_E_preHit = final_E_preHit(1:idx-1);
    final_E_postHit = final_E_postHit(1:idx-1);
    results_NumberOfHits = results_NumberOfHits(1:idx-1, :);
    initial_conditions_preHit = initial_conditions_preHit(1:idx-1, :);
    initial_conditions_postHit = initial_conditions_postHit(1:idx-1, :);
    
    %% calc the Angle-Action (AA) coords:

    AA_preHit_results = convertToAngleActionCoordsOnHarmonicHamiltonian(results_preHit, w0);
    AA_postHit_results = convertToAngleActionCoordsOnHarmonicHamiltonian(results_postHit, w0);

    
    %% set the return struct
    Estruct = struct();
    Estruct.results_preHit = results_preHit;
    Estruct.results_postHit = results_postHit;
    Estruct.final_E_preHit = final_E_preHit;
    Estruct.final_E_postHit = final_E_postHit;
    Estruct.results_NumberOfHits = results_NumberOfHits;
    Estruct.initial_conditions_preHit = initial_conditions_preHit;
    Estruct.initial_conditions_postHit = initial_conditions_postHit;
    Estruct.AA_preHit_results = AA_preHit_results;
    Estruct.AA_postHit_results = AA_postHit_results;
end