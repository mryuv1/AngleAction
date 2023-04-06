function [position,isterminal,direction] = EventsFunction(t,q) % q = [x, y, px, py];
    
    myParams = params();
    circle_center = myParams.circle_center;
    radius = myParams.radius;
    % check that the track isn't inside the step, with euiler aproximation
    inside_circ = ((q(1) - circle_center(1))^2 + (q(2) - circle_center(2))^2) +eps < radius^2;
    inside_vertial_wall = (q(1) < (circle_center(1) + radius)) && (q(2) < (circle_center(2)));
    inside_horizontal_wall = (q(2) < (circle_center(2) + radius)) && (q(1) < circle_center(1));
%     if (inside_circ) % || inside_vertial_wall || inside_horizontal_wall)
%         disp(inside_circ | inside_vertial_wall | inside_horizontal_wall)
%     end
    position = [q(3); ~inside_circ; ~inside_vertial_wall; ~inside_horizontal_wall];             % The value that we want to be zero
    isterminal = [1; 1; 1; 1];   % Stop the integration
    direction = [-1; 0; 0; 0];   % Negative direction only
end