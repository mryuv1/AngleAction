function totalE = calcTotalEnergy_V2(q, parameters) % q = [x, y, px, py];
    myParams = parameters;
    totalE = 0.5*(q(3)^2 + q(4)^2) + 0.5*((myParams.w1*q(1))^2 + (myParams.w2*q(2))^2);
end
