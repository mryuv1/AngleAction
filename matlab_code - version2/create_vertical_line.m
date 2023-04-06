function coords = create_vertical_line(x_cord, ylims, numberOfPoints) % ylims = [y_min, y_max];  
    if nargin<3
        numberOfPoints = 1e2;
    end

    xs = x_cord.*ones(numberOfPoints, 1);
    ys = linspace(ylims(1), ylims(2), numberOfPoints);
    m = -inf.*ones(numberOfPoints, 1);

    coords = [xs ys.' m];
end