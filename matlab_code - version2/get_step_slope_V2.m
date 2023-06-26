function slope_ = get_step_slope_V2(q, parameters)
    
    inside = checkIfInsideEuler_V2(q, parameters);
    myParams = parameters;

    if inside==1
        slope_ = 0;
    elseif inside==2
        slope_ = (myParams.circle_center(1) - q(1))./ ...
                (q(2) - myParams.circle_center(2) + eps);
    elseif inside==3
        slope_ = -inf;
    end
end
