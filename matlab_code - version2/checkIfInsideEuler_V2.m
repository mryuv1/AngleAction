function inside = checkIfInsideEuler_V2(q, parameters) % q = [x, y, px, py];
    epsilon = 1e-4;
    myParams = parameters;
    circle_center = myParams.circle_center;
    radius = myParams.radius;
    % check that the track isn't inside the step, with euiler aproximation
    inside_circ = ((q(1) + q(3)*epsilon) - circle_center(1))^2 + ((q(2) + q(4)*epsilon) - circle_center(2))^2 < radius^2;
    inside_vertial_wall = ((q(1) + q(3)*epsilon) < (circle_center(1) + radius)) && ((q(2) + q(4)*epsilon) < (circle_center(2)));
    inside_horizontal_wall = ((q(2) + q(4)*epsilon) < (circle_center(2) + radius)) && ((q(1) + q(3)*epsilon) < (circle_center(1)));
    
    if inside_horizontal_wall
        inside = 1;
    elseif inside_circ
        inside = 2;
    elseif inside_vertial_wall
        inside = 3;
    else
        inside = 0;
    end
    
end