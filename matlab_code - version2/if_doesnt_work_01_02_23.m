
clc;
clear;
close all;

%% variables

parameters = params();
w1 = parameters.w1;
w2 = parameters.w2;
Et = parameters.Et;
Ex = parameters.Ex;
Ey = Et - Ex;
circle_center = parameters.circle_center;
radius = parameters.radius;
limit_number_of_round_hits = 1; %boolian for limiting number of hits
number_of_round_hits_treshold = 1;

%% check for one point

% r = 0.2;
% x_w = -0.4;
% y_w = -0.6;

% make line

% ys = linspace(-1,y_w, 1000);
% x = ones(1,1000)*(-0.4);
% slope = ones(1,1000)*(inf);
% w = [x.' ys.' slope.'];

% Ex_min = 0.5*(w1^2)*(x_w^2);
% EY_min = 0.5*(w2^2)*(y_w^2);
numberOfenergies = 1e2;
% E_x = Ex_min:0.05:Et;


vertical_step = create_vertical_line(circle_center(2) + radius, [-1 circle_center(1)]);
horizontal_step = create_horizontal_line(circle_center(1) + radius, [-1 circle_center(2)]);
coords_circ = create_cuarter_circle(circle_center, radius, 5e2);

united_coords = [horizontal_step; coords_circ; vertical_step];

w0 = [w1, w2];
results_preHit = zeros(length(coords_circ)*4, 4);
results_postHit = zeros(length(coords_circ)*4, 4);
final_E_preHit = zeros(length(coords_circ)*4, 1);
final_E_postHit = zeros(length(coords_circ)*4, 1);
idx = 1;

coords_to_use = united_coords;

for w_idc = 1:length(coords_to_use)
    x_w = coords_to_use(w_idc, 1);
    y_w = coords_to_use(w_idc, 2);
    Ex_min = 0.5*(w1^2)*(x_w^2);
    EY_min = 0.5*(w2^2)*(y_w^2);

    if Et-EY_min<Ex_min
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

    tsapn = [0 50];
    opt = odeset('RelTol',1e-9, 'Events',@EventsFunction);
    
    flag_track_right_up = checkIfInside(q0_right_up);
    flag_track_right_down = checkIfInside(q0_right_down);
    flag_track_left_up = checkIfInside(q0_left_up);
    flag_track_left_down = checkIfInside(q0_left_down);

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
        [q0_track_right_up, q0_hits_round_right_up] ...
            = calc_track_until_px_is_zero(q0_right_up, w0);
        final_PreHit_right_up = q0_track_right_up(end, :);
        p_postHit_right_up = return_momentum(-[q0_right_up(3) q0_right_up(4)], coords_to_use(w_idc, 3));
        q1_right_up = [x_w; y_w; p_postHit_right_up(1); p_postHit_right_up(2)];
        [q1_track_right_up, q1_hits_round_right_up] ...
            = calc_track_until_px_is_zero(q1_right_up, w0);
        final_PostHit_right_up = q1_track_right_up(end, :);
        
%         figure;
%         hold on
%         plot(coords_to_use(:,1),coords_to_use(:,2),'ro')
%         plot(q0_track_right_up(:,1),q0_track_right_up(:,2),'bo')
%         plot(q1_track_right_up(:,1),q1_track_right_up(:,2),'go')
%         hold off
        if (~limit_number_of_round_hits) || ...
        ((q0_hits_round_right_up < number_of_round_hits_treshold) && ...
        (q1_hits_round_right_up < number_of_round_hits_treshold))
            % locate results
            final_E_preHit(idx) = calcTotalEnergy(final_PreHit_right_up);
            final_E_postHit(idx) = calcTotalEnergy(final_PostHit_right_up);
            results_preHit(idx, :) = final_PreHit_right_up;
            results_postHit(idx, :) = final_PostHit_right_up;
            idx = idx + 1;
        end
    end
    
    
    if ~flag_track_right_down
        [q0_track_right_down, q0_hits_right_down] ...
            = calc_track_until_px_is_zero(q0_right_down, w0);
        final_PreHit_right_down = q0_track_right_down(end, :);
        
%         if idx==232
%             disp(idx)
%         end
% 
%         if idx==366
%             disp(idx)
%         end
        p_postHit_right_down = return_momentum(-[q0_right_down(3) q0_right_down(4)], coords_to_use(w_idc, 3));
        q1_right_down = [x_w; y_w; p_postHit_right_down(1); p_postHit_right_down(2)];
%         disp(idx)
        [q1_track_right_down, q1_hits_right_down] = ...
            calc_track_until_px_is_zero(q1_right_down, w0);
        
        final_PostHit_right_down = q1_track_right_down(end, :);
%         
%         figure;
%         hold on
%         plot(coords_to_use(:,1),coords_to_use(:,2),'ro')
%         plot(q0_track_right_down(:,1),q0_track_right_down(:,2),'bo')
%         plot(q1_track_right_down(:,1),q1_track_right_down(:,2),'go')
%         hold off
        
        % locate results
        if (~limit_number_of_round_hits) || ...
        ((q0_hits_right_down < number_of_round_hits_treshold) && ...
        (q1_hits_right_down < number_of_round_hits_treshold))
            final_E_preHit(idx) = calcTotalEnergy(final_PreHit_right_down);
            final_E_postHit(idx) = calcTotalEnergy(final_PostHit_right_down);
            results_preHit(idx, :) = final_PreHit_right_down;
            results_postHit(idx, :) = final_PostHit_right_down;
            idx = idx + 1;
        end
    end


    
    if ~flag_track_left_up
        [q0_track_left_up, q0_hits_left_up] ...
            = calc_track_until_px_is_zero(q0_left_up, w0);
        final_PreHit_left_up = q0_track_left_up(end, :);
        p_postHit_left_up = return_momentum(-[q0_left_up(3) q0_left_up(4)], coords_to_use(w_idc, 3));
        q1_left_up = [x_w; y_w; p_postHit_left_up(1); p_postHit_left_up(2)];
        [q1_track_left_up, q1_hits_left_up] ...
            = calc_track_until_px_is_zero(q1_left_up, w0);
        
        final_PostHit_left_up = q1_track_left_up(end, :);
        % locate results
        if (~limit_number_of_round_hits) || ...
        ((q0_hits_left_up < number_of_round_hits_treshold) && ...
            (q1_hits_left_up < number_of_round_hits_treshold))
            final_E_preHit(idx) = calcTotalEnergy(final_PreHit_left_up);
            final_E_postHit(idx) = calcTotalEnergy(final_PostHit_left_up);
            results_preHit(idx, :) = final_PreHit_left_up;
            results_postHit(idx, :) = final_PostHit_left_up;
            idx = idx + 1;
        end
    end
    
    
    if ~flag_track_left_down

        [q0_track_left_down, q0_hits_left_down] = calc_track_until_px_is_zero(q0_left_down, w0);
        final_PreHit_left_down = q0_track_left_down(end, :);
        p_postHit_left_down = return_momentum(-[q0_left_down(3) q0_left_down(4)], coords_to_use(w_idc, 3));
        q1_left_down = [x_w; y_w; p_postHit_left_down(1); p_postHit_left_down(2)];
        [q1_track_left_down, q1_hits_left_down] = calc_track_until_px_is_zero(q1_left_down, w0);
        
        final_PostHit_left_down = q1_track_left_down(end, :);

        % locate results
        if (~limit_number_of_round_hits) || ...
        ((q0_hits_left_up < number_of_round_hits_treshold) && ...
        (q1_hits_left_up < number_of_round_hits_treshold))
            
            final_E_preHit(idx) = calcTotalEnergy(final_PreHit_left_down);
            final_E_postHit(idx) = calcTotalEnergy(final_PostHit_left_down);
            results_preHit(idx, :) = final_PreHit_left_down;
            results_postHit(idx, :) = final_PostHit_left_down;
            idx = idx + 1;
        end
    end

end

results_preHit = results_preHit(1:idx-1, :);
results_postHit = results_postHit(1:idx-1, :);
final_E_preHit = final_E_preHit(1:idx-1);
final_E_postHit = final_E_postHit(1:idx-1);

% if energy isn't conserved
assert(var((final_E_preHit - final_E_postHit)./Et)<1e-5, "The Energy is NOT conserved");

Energy_Y_line_coords = calcY_Energy_level();
figure;
hold on
plot(Energy_Y_line_coords(:,1),Energy_Y_line_coords(:,2), 'k')
plot(results_preHit(:,2),results_preHit(:,4), 'rx')
plot(results_postHit(:,2),results_postHit(:,4), 'bo')
grid on
legend('Y preHit energy levels','preHit', 'PostHit')
xlabel("y")
ylabel('P_y')
hold off
xlim([-3.5, 3.5]);
ylim([-3.5, 3.5]);

% figure;
% sanity_check = results_preHit(:,2).^2 + results_preHit(:,4).^2;
% plot(1:length(sanity_check),sanity_check, '+');
%% functions:


% function [position,isterminal,direction] = EventsFunction(t,q) % q = [x, y, px, py];
%     
%     position = q(3);             % The value that we want to be zero
%     isterminal = 1;   % Stop the integration
%     direction = -1;   % Negative direction only
% end













