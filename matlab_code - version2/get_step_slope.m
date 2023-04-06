function slope_ = get_step_slope(q)
    
    inside = checkIfInsideEuler(q);
    myParams = params();

    if inside==1
        slope_ = 0;
    elseif inside==2
        slope_ = (myParams.circle_center(1) - q(1))./ ...
                (q(2) - myParams.circle_center(2) + eps);
    elseif inside==3
        slope_ = -inf;
    end
end
