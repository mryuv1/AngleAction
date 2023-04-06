function myParams = params()
    myParams = struct;
    myParams.w1 = 1;
    myParams.w2 = sqrt(2);
    myParams.Et = 2*(myParams.w1^2 + myParams.w2^2)*2;
    myParams.Ex = 1.5;
    myParams.Ey = myParams.Et - myParams.Ex;
    myParams.circle_center = [-0.4 -0.4];
    myParams.radius = 0.1;
    myParams.limit_number_of_round_hits = 1; % max number of hits
    myParams.number_of_round_hits_treshold = 1;
end