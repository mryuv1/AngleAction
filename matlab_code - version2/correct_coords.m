function [qc, roundPartHit, horizontalHit, verticalHit] = correct_coords(position)
    myParams = params();
    qc = position;
    horizontalHit = 0;
    roundPartHit = 0;
    verticalHit = 0;
    switch checkIfInsideEuler(position)
        case 1      % inside_horizontal_wall
            qc(2) = myParams.circle_center(1) + myParams.radius;
            horizontalHit = 1;
        case 2      % inside_circ
            ang = atan((qc(2) - myParams.circle_center(2))./ ...
                        ((qc(1) - myParams.circle_center(1))));
            qc(1) = myParams.circle_center(1) + myParams.radius*cos(ang);
            qc(2) = myParams.circle_center(2) + myParams.radius*sin(ang);
            roundPartHit = 1;
        case 3      % inside_vertial_wall
            qc(1) = myParams.circle_center(2) + myParams.radius;
            verticalHit = 1;
        otherwise
    %         disp(final_position);
    end

end