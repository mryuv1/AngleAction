function Bstruct = get_circle_border_coords(y_energy_min, y_energy_max, parameters, numberOfEnergies)
    
    if nargin<4
        numberOfEnergies = 1e2;
    end

    circle_center = parameters.circle_center;
    w1 = parameters.w1;
    w2 = parameters.w2;
    w0 = [w1 w2];
    Et = parameters.Et;
    
    % the two relevant points:
    point_right = [circle_center(1) + parameters.radius, circle_center(2)];
    point_up = [circle_center(1), circle_center(2) + parameters.radius];
    
    % the energies for the iteration
    energies = linspace(y_energy_min, y_energy_max, numberOfEnergies);

    results_preHit_up = zeros(length(energies)*2, 4);
    results_preHit_right = zeros(length(energies)*2, 4);
    results_postHit_up = zeros(length(energies)*2, 4);
    results_postHit_right = zeros(length(energies)*2, 4);
    final_E_preHit_right = zeros(length(energies)*2, 1);
    final_E_preHit_up = zeros(length(energies)*2, 1);
    
    results_postHit_right_p = zeros(length(energies)*2, 4);
    results_postHit_right_m = zeros(length(energies)*2, 4);

    idx_up = 1;
    idx_right = 1;
    
    % check
    idx_right_p = 1;
    idx_right_m = 1;
    for enidx = 1:(length(energies)-1)

        % if enidx==100
        %     display(enidx)
        % end
    

        % for the right point
        Ey_right = energies(enidx);
        Ex_right = Et - Ey_right;

        Vx_right = 0.5*(w1^2)*(point_right(1)^2);
        Vy_right = 0.5*(w2^2)*(point_right(2)^2);
        
        px_right = sqrt(2*(Ex_right - Vx_right));
        py_right = sqrt(2*(Ey_right - Vy_right));
            
        q0_right_plus = [point_right(1); point_right(2); px_right; py_right];
        q0_right_minus = [point_right(1); point_right(2); px_right; -py_right];
        
        if isreal(q0_right_plus)
             [q1_track_right_plus,~, ~, ~] ...
                = calc_track_until_px_is_zero_V2(q0_right_plus, w0, 0, parameters);
            final_PreHit_right_plus = q1_track_right_plus(end, :);
            final_PreHit_right_plus(3:4) = -final_PreHit_right_plus(3:4);

            p_postHit_right_plus = return_momentum(-[q0_right_plus(3) q0_right_plus(4)], -inf);
            q2_right_plus = [point_right(1); point_right(2); p_postHit_right_plus(1); p_postHit_right_plus(2)];
            [q2_track_right_plus, ~, ~, ~] ...
                = calc_track_until_px_is_zero_V2(q2_right_plus, w0, 1, parameters);
            final_PostHit_right_plus = q2_track_right_plus(end, :);
            
            
            % for the first option:
            final_E_preHit_right(idx_right) = calcTotalEnergy_V2(final_PreHit_right_plus, parameters);
            results_preHit_right(idx_right, :) = final_PreHit_right_plus;
            results_postHit_right(idx_right, :) = final_PostHit_right_plus;
            % check
            results_postHit_right_p(idx_right_p, :) = final_PostHit_right_plus;
            idx_right_p = idx_right_p + 1;
            idx_right = idx_right + 1;

            % figure;
            % hold on
            % plot(q1_track_right_plus(:,1),q1_track_right_plus(:,2),'bo')
            % plot(q2_track_right_plus(:,1),q2_track_right_plus(:,2),'ro')
            % hold off
        end

        % for the second option:
        if isreal(q0_right_minus)
            [q1_track_right_minus, ~, ~, ~] ...
                = calc_track_until_px_is_zero_V2(q0_right_minus, w0, 0, parameters);
            final_PreHit_right_minus = q1_track_right_minus(end, :);
            final_PreHit_right_minus(3:4) = -final_PreHit_right_minus(3:4);

            p_postHit_right_plus = return_momentum(-[q0_right_minus(3) q0_right_minus(4)], -inf);
            q2_right_minus = [point_right(1); point_right(2); p_postHit_right_plus(1); p_postHit_right_plus(2)];
            [q2_track_right_minus, ~, ~, ~] ...
                = calc_track_until_px_is_zero_V2(q2_right_minus, w0, 1, parameters);
            final_PostHit_right_minus = q2_track_right_minus(end, :);
                
            final_E_preHit_right(idx_right) = calcTotalEnergy_V2(final_PreHit_right_minus, parameters);
            results_preHit_right(idx_right, :) = final_PreHit_right_minus;
            results_postHit_right(idx_right, :) = final_PostHit_right_minus;
            % check
            results_postHit_right_m(idx_right, :) = final_PostHit_right_minus;
            idx_right_m = idx_right_m + 1;

            idx_right = idx_right + 1;

            % figure;
            % hold on
            % plot(q1_track_right_minus(:,1),q1_track_right_minus(:,2),'bo')
            % plot(q2_track_right_minus(:,1),q2_track_right_minus(:,2),'ro')
            % hold off
        end

        % for the up point
        Ey_up = energies(enidx);
        Ex_up = Et - Ey_up;

        Vx_up = 0.5*(w1^2)*(point_up(1)^2);
        Vy_up = 0.5*(w2^2)*(point_up(2)^2);
        
        px_up = sqrt(2*(Ex_up - Vx_up));
        py_up = sqrt(2*(Ey_up - Vy_up));
        
        q0_up_plus = [point_up(1); point_up(2); -px_up; py_up];
        q0_up_minus = [point_up(1); point_up(2); px_up; py_up];
        
        if isreal(q0_up_plus)
            [q1_track_up_plus,~, ~, ~] ...
                = calc_track_until_px_is_zero_V2(q0_up_plus, w0, 0, parameters);
            final_PreHit_up_plus = q1_track_up_plus(end, :);
            final_PreHit_up_plus(3:4) = -final_PreHit_up_plus(3:4);

            p_postHit_up_plus = return_momentum(-[q0_up_plus(3) q0_up_plus(4)], 0);
            q2_right_up = [point_up(1); point_up(2); p_postHit_up_plus(1); p_postHit_up_plus(2)];
            [q2_track_up_plus, ~, ~, ~] ...
                = calc_track_until_px_is_zero_V2(q2_right_up, w0, 1, parameters);
            final_PostHit_up_plus = q2_track_up_plus(end, :);

            % for the first option:
            final_E_preHit_up(idx_up) = calcTotalEnergy_V2(final_PreHit_up_plus, parameters);
            results_preHit_up(idx_up, :) = final_PreHit_up_plus;
            results_postHit_up(idx_up, :) = final_PostHit_up_plus;
            idx_up = idx_up + 1;

            % figure;
            % hold on
            % plot(q1_track_up_plus(:,1),q1_track_up_plus(:,2),'bo')
            % plot(q2_track_up_plus(:,1),q2_track_up_plus(:,2),'ro')
            % hold off
        end
        
        if isreal(q0_up_minus)
        % for the second option:
            [q1_track_up_minus, ~, ~, ~] ...
                = calc_track_until_px_is_zero_V2(q0_up_minus, w0, 0, parameters);
            final_PreHit_up_minus = q1_track_up_minus(end, :);
            final_PreHit_up_minus(3:4) = -final_PreHit_up_minus(3:4);

            p_postHit_up_minus = return_momentum(-[q0_up_minus(3) q0_up_minus(4)], 0);
            q2_right_minus = [point_up(1); point_up(2); p_postHit_up_minus(1); p_postHit_up_minus(2)];
            [q2_track_up_minus, ~, ~, ~] ...
                = calc_track_until_px_is_zero_V2(q2_right_minus, w0, 1, parameters);
            final_PostHit_up_minus = q2_track_up_minus(end, :);

            final_E_preHit_up(idx_up) = calcTotalEnergy_V2(final_PreHit_up_minus, parameters);
            results_preHit_up(idx_up, :) = final_PreHit_up_minus;
            results_postHit_up(idx_up, :) = final_PostHit_up_minus;
            idx_up = idx_up + 1;

            % figure;
            % hold on
            % plot(q1_track_up_minus(:,1),q1_track_up_minus(:,2),'bo')
            % plot(q2_track_up_minus(:,1),q2_track_up_minus(:,2),'ro')
            % hold off
        end      
    end

% make results in the right side:
results_preHit_up = results_preHit_up(1:idx_up-1, :);
results_postHit_up = results_postHit_up(1:idx_up-1, :);
final_E_preHit_up = final_E_preHit_up(1:idx_up-1);
results_preHit_right = results_preHit_right(1:idx_right-1, :);
results_postHit_right = results_postHit_right(1:idx_right-1, :);
final_E_preHit_right = final_E_preHit_right(1:idx_right-1);

% check
results_postHit_right_p = results_postHit_right_m((1:idx_right_p-1), :);
results_postHit_right_m = results_postHit_right_m((1:idx_right_m-1), :);

AA_postHit_results_right_p = convertToAngleActionCoordsOnHarmonicHamiltonian(results_postHit_right_p, w0);
AA_postHit_results_right_m = convertToAngleActionCoordsOnHarmonicHamiltonian(results_postHit_right_m, w0);
% calc the Angle-Action (AA) coords:

AA_preHit_results_up = convertToAngleActionCoordsOnHarmonicHamiltonian(results_preHit_up, w0);
AA_preHit_results_right = convertToAngleActionCoordsOnHarmonicHamiltonian(results_preHit_right, w0);

AA_postHit_results_up = convertToAngleActionCoordsOnHarmonicHamiltonian(results_postHit_up, w0);
AA_postHit_results_right = convertToAngleActionCoordsOnHarmonicHamiltonian(results_postHit_right, w0);

% set the return struct
Bstruct = struct();
Bstruct.results_preHit_up = results_preHit_up;
Bstruct.results_postHit_up = results_postHit_up;
Bstruct.final_E_preHit_up = final_E_preHit_up;
Bstruct.AA_preHit_results_up = AA_preHit_results_up;
Bstruct.AA_postHit_results_up = AA_postHit_results_up;
Bstruct.results_preHit_right = results_preHit_right;
Bstruct.results_postHit_right = results_postHit_right;
Bstruct.final_E_preHit_right = final_E_preHit_right;
Bstruct.AA_preHit_results_right = AA_preHit_results_right;
Bstruct.AA_postHit_results_right = AA_postHit_results_right;
end