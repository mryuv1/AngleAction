function Energy_line_coords = calcY_Energy_level(parameters) 
    % calc the y, py coords
    myParams = parameters;
    y_max = sqrt(2*myParams.Ey/(myParams.w2^2));
    py_max = sqrt(2*myParams.Ey);
    thetas = linspace(0, 2*pi, 1e2);

    xs = y_max*cos(thetas);
    ys = py_max*sin(thetas);

    Energy_line_coords = [xs.' ys.'];
end