function totalE = calcTotalEnergy(q) % q = [x, y, px, py];
    myParams = params();
    totalE = 0.5*(q(3)^2 + q(4)^2) + 0.5*((myParams.w1*q(1))^2 + (myParams.w2*q(2))^2);
end
