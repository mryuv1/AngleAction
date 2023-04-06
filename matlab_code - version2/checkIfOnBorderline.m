function [horizontalHit, hitOnCircle, verticalHit] = checkIfOnBorderline(q) % q = [x, y, px, py];

    myParams = params();
    circle_center = myParams.circle_center;
    radius = myParams.radius;
    % check that the track isn't inside the step, with euiler aproximation
    inside_circ = abs((q(1) - circle_center(1))^2 + (q(2) - circle_center(2))^2 - radius^2) < eps;
    inside_vertial_wall = (q(1) == (circle_center(1) + radius)) && (q(2) <= (circle_center(2)));
    inside_horizontal_wall = (q(2) == (circle_center(2) + radius)) && (q(1) < (circle_center(1)));
    
    hitOnCircle = 0;
    horizontalHit = 0;
    verticalHit = 0;

    if inside_horizontal_wall
        horizontalHit = 1;
    elseif inside_circ
        hitOnCircle = 1;
    elseif inside_vertial_wall
        verticalHit = 1;
    end
    
end