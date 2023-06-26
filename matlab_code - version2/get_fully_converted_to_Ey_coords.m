function Cstruct = get_fully_converted_to_Ey_coords(circle_coords, parameters)

    circle_center = parameters.circle_center;
    w1 = parameters.w1;
    w2 = parameters.w2;
    w0 = [w1 w2];
    Et = parameters.Et;
    
    results_preHit = zeros(length(circle_coords), 4);
    results_postHit = zeros(length(circle_coords), 4);
    final_E_preHit = zeros(length(circle_coords), 1);
    final_E_postHit = zeros(length(circle_coords), 1);   
    
    idx = 1;
    
    for w_idc = 1:length(circle_coords)
        x_w = circle_coords(w_idc, 1);
        y_w = circle_coords(w_idc, 2);

        Vx = 0.5*(w1^2)*(x_w^2);
        Vy = 0.5*(w2^2)*(y_w^2);

        Vtot = Vx + Vy;
        Ktot = Et - Vtot;
        
        % if w_idc==87
        %     display(w_idc)
        % end
        % 
        % if w_idc==86
        %     display(w_idc)
        % end
        
        % calc the postHit coords:
        Ex_after = Vx;
        Ey_after = Et - Ex_after;
        py_after = sqrt(2*(Ey_after - Vy));

        q0_after = [x_w; y_w; 0; py_after];
        p_postHit_right_up = return_momentum(-[q0_after(3) q0_after(4)], circle_coords(w_idc, 3));

        q0_before = [x_w; y_w; p_postHit_right_up(1); p_postHit_right_up(2)];
        
        [q0_track_before,~, ~, ~] ...
            = calc_track_until_px_is_zero_V2(q0_before, w0, 0, parameters);
        final_PreHit_before = q0_track_before(end, :);
        final_PreHit_before(3:4) = -final_PreHit_before(3:4);

        [q1_track_after, ~, ~, ~] ...
            = calc_track_until_px_is_zero_V2(q0_after, w0, 0, parameters);
        final_postHit_after = q1_track_after(end, :);
                
        % % plot results
        %     figure;
        %     hold on
        %     plot(circle_coords(:,1),circle_coords(:,2),'ro')
        %     plot(q0_track_before(:,1),q0_track_before(:,2),'bo')
        %     plot(q1_track_after(:,1),q1_track_after(:,2),'go')
        %     legend('step', 'preHit', 'postHit')
        %     xlim([-10,10])
        %     ylim([-10,10])
        %     hold off

        % for the first option:
        final_E_preHit(idx) = calcTotalEnergy_V2(final_PreHit_before, parameters);
        final_E_postHit(idx) = calcTotalEnergy_V2(final_postHit_after, parameters);
        results_preHit(idx, :) = final_PreHit_before;
        results_postHit(idx, :) = final_postHit_after;

        idx = idx + 1;
   end

% make results in the right side:
results_preHit = results_preHit(1:idx-1, :);
final_E_preHit = final_E_preHit(1:idx-1);
results_postHit = results_postHit(1:idx-1, :);
final_E_postHit = final_E_postHit(1:idx-1);

% calc the Angle-Action (AA) coords:

AA_preHit_results = convertToAngleActionCoordsOnHarmonicHamiltonian(results_preHit, w0);
AA_postHit_results = convertToAngleActionCoordsOnHarmonicHamiltonian(results_postHit, w0);

% set the return struct
Cstruct = struct();
Cstruct.results_preHit = results_preHit;
Cstruct.final_E_preHit = final_E_preHit;
Cstruct.results_postHit = results_postHit;
Cstruct.final_E_postHit = final_E_postHit;
Cstruct.AA_preHit_results = AA_preHit_results;
Cstruct.AA_postHit_results = AA_postHit_results;

end