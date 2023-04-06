function coords = create_horizontal_line(y_cord, xlims, numberOfPoints) % ylims = [y_min, y_max];  
    if nargin<3
        numberOfPoints = 1e2;
    end

    ys = y_cord.*ones(numberOfPoints, 1);
    xs = linspace(xlims(1), xlims(2), numberOfPoints);
    m = zeros(numberOfPoints, 1);

    coords = [xs.' ys m];
    
end