function coords = create_cuarter_circle(center, radius, numberOfPoints) % center = [x, y]; coords = [x, y, slope]  
    if nargin<3
        numberOfPoints = 1e2;
    end
    
    thetas = linspace(0, pi/2, numberOfPoints);
    xs = center(1) + radius*cos(thetas);
    ys = center(2) + radius*sin(thetas);
    m = (center(1) - xs)./(ys - center(2) + eps);

    coords = [xs.' ys.' m.'];
    
end