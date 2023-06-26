function [track, numOfRoundHit, numOfHorizontalHit, numOfVerticalHit] = calc_track_until_px_is_zero_V2(initial_conds, w, isSecondHit, parameters)
    
    if nargin<3
        isSecondHit = 0;
    end

    if (~isSecondHit)
            [numOfHorizontalHit, numOfRoundHit, numOfVerticalHit] = ...
                                checkIfOnBorderline_V2(initial_conds, parameters);
    else
        numOfHorizontalHit = 0;
        numOfRoundHit = 0;
        numOfVerticalHit = 0;
    end

    tsapn = [0 50];
    stop_integration_Function = @(t, q) EventsFunction_V2(t, q, parameters);
    opt = odeset('AbsTol',5e-11,'RelTol',5e-11, 'Events',stop_integration_Function);
    [~,track] = ode45(@(t,q)moveMentEq(t, q, w),tsapn, initial_conds, opt);
    final_position = track(end, :);
    
    [final_position, hitOnCircle, horizontalHit, verticalHit] = ...
                                        correct_coords_V2(final_position, parameters);
    numOfRoundHit = numOfRoundHit + hitOnCircle;
    numOfHorizontalHit = numOfHorizontalHit + horizontalHit;
    numOfVerticalHit = numOfVerticalHit + verticalHit;
    
    while abs(final_position(3))>10e-13
        slope = get_step_slope_V2(final_position, parameters);
        p_new = return_momentum([final_position(3) final_position(4)], slope);
        new_initials = [final_position(1) final_position(2) p_new(1) p_new(2)];
        [~,track_tmp] = ode45(@(t,q)moveMentEq(t, q, w),tsapn, new_initials, opt);
        track = [track; track_tmp];
        final_position = track_tmp(end, :);
        [final_position, hitOnCircle, horizontalHit, verticalHit] = correct_coords_V2(final_position, parameters);
        numOfRoundHit = numOfRoundHit + hitOnCircle;
        numOfHorizontalHit = numOfHorizontalHit + horizontalHit;
        numOfVerticalHit = numOfVerticalHit + verticalHit;

%         figure;
%         hold on
%         plot(step_coords.coords_circ(:,1),step_coords.coords_circ(:,2),'ro')
% %         plot(track_tmp(:,1),track_tmp(:,2),'ro')
%         plot(track(:,1),track(:,2),'bo')
%         hold off
    end
end