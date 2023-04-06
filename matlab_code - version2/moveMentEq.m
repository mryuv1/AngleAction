function dqdt = moveMentEq(t, q, w) % q = [x, y, px, py]; w = [w1, w2]
    
    x = q(1);
    y = q(2);
    px = q(3);
    py = q(4);

    dqdt = zeros(4,1);
    dqdt(1) = px; % dx/dt
    dqdt(2) = py; % dy/dt
    dqdt (3) = -(w(1)^2)*x; % dPx/dt
    dqdt (4) = -(w(2)^2)*y; % dPy/dt
end