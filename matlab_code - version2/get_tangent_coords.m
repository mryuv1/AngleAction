function Bstruct = get_tangent_coords(circle_coords, parameters)

    circle_center = parameters.circle_center;
    w1 = parameters.w1;
    w2 = parameters.w2;
    Et = parameters.Et;

    for w_idc = 1:length(circle_coords)
        x_w = coords_to_use(w_idc, 1);
        y_w = coords_to_use(w_idc, 2);
        m_traj = (circle_center(1) - x_w)/(y_w - circle_center(2) + eps);
        
        Vx = 0.5*(w1^2)*(x_w^2);
        Vy = 0.5*(w2^2)*(y_w^2);

        Vtot = Vx + Vy;
        Ktot = Et - Vtot;

    end

end