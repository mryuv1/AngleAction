function p_new = return_momentum(p, slope) % q = [px, py];

    theta = atan(slope);
    p_new = zeros(2,1);           
    p_new(1) = p(1)*cos(2*theta) + p(2)* sin(2*theta);
    p_new(2) = p(1)*sin(2*theta) - p(2)* cos(2*theta);
end