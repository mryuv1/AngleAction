function Bstruct = get_tangent_coords(circle_coords, parameters)

    circle_center = parameters.circle_center;
    w1 = parameters.w1;
    w2 = parameters.w2;
    w0 = [w1 w2];
    Et = parameters.Et;
    
    results_preHit = zeros(length(circle_coords)*2, 4);
    results_postHit = zeros(length(circle_coords)*2, 4);
    final_E_preHit = zeros(length(circle_coords)*2, 1);

    
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
        % if idx==24 || idx==23
        %     display(w_idc)
        % end

        if (y_w - circle_center(2))==0
            m_traj = inf;
            px = 0;
            py = sqrt(2*Ktot);
        else
            m_traj = (circle_center(1) - x_w)/(y_w - circle_center(2));
            px = sqrt((2*Ktot)/(1 + m_traj^2));
            py = m_traj*px;
        end

        

        q0_plus = [x_w; y_w; px; py];
        q0_minus = [x_w; y_w; -px; -py];

        [q0_track_plus,~, ~, ~] ...
            = calc_track_until_px_is_zero_V2(q0_plus, w0, 0, parameters);
        final_PreHit_plus = q0_track_plus(end, :);
        final_PostHit_minus = final_PreHit_plus;
        final_PreHit_plus(3:4) = -final_PreHit_plus(3:4);


        [q1_track_minus, ~, ~, ~] ...
            = calc_track_until_px_is_zero_V2(q0_minus, w0, 0, parameters);
        final_PreHit_minus = q1_track_minus(end, :);
        final_PostHit_plus = final_PreHit_minus;
        final_PreHit_minus(3:4) = -final_PreHit_minus(3:4);

        % for the first option:
        final_E_preHit(idx) = calcTotalEnergy_V2(final_PreHit_plus, parameters);
        results_preHit(idx, :) = final_PreHit_plus;
        results_postHit(idx, :) = final_PostHit_plus;

        idx = idx + 1;

        % for the second option:
        final_E_preHit(idx) = calcTotalEnergy_V2(final_PreHit_minus, parameters);
        results_preHit(idx, :) = final_PreHit_minus;
        results_postHit(idx, :) = final_PostHit_minus;
        idx = idx + 1;
        
        % if mod(w_idc, 10)==0
            % figure;
            % hold on
            % plot(circle_coords(:,1),circle_coords(:,2),'ro')
            % plot(q0_track_plus(:,1),q0_track_plus(:,2),'bo')
            % plot(q1_track_minus(:,1),q1_track_minus(:,2),'go')
            % legend('step', 'preHit', 'postHit')
            % xlim([-1,1])
            % ylim([-1,1])
            % hold off
        % end

   end

% make results in the right side:
results_preHit = results_preHit(1:idx-1, :);
results_postHit = results_postHit(1:idx-1, :);
final_E_preHit = final_E_preHit(1:idx-1);


% calc the Angle-Action (AA) coords:

AA_preHit_results = convertToAngleActionCoordsOnHarmonicHamiltonian(results_preHit, w0);
AA_postHit_results = convertToAngleActionCoordsOnHarmonicHamiltonian(results_postHit, w0);

% set the return struct
Bstruct = struct();
Bstruct.results_preHit = results_preHit;
Bstruct.results_postHit = results_postHit;
Bstruct.final_E_preHit = final_E_preHit;
Bstruct.AA_preHit_results = AA_preHit_results;
Bstruct.AA_postHit_results = AA_postHit_results;
end