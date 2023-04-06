
clc;
clear;
close all;

%% Debug variables

plot_indexes = 0;
color = ["r" "g" "b" "c" "m" "y"];
marker = ["o" "+" "*" "x" "square" "diamond" "pentagram" "^" "hexagram"];

%% variables

parameters = params();
w1 = parameters.w1;
w2 = parameters.w2;
Et = parameters.Et;
Ex = parameters.Ex;
Ey = Et - Ex;
circle_center = parameters.circle_center;
radius = parameters.radius;
limit_number_of_round_hits = parameters.limit_number_of_round_hits;
%boolian for limiting number of hits
number_of_round_hits_treshold = parameters.number_of_round_hits_treshold;

min_possible_Ex = 0.5*(parameters.w1*(parameters.circle_center(1) + parameters.radius))^2;
min_possible_Ey = 0.5*(parameters.w2*(parameters.circle_center(2) + parameters.radius))^2;
max_possible_Ex = Et - min_possible_Ey;
max_possible_Ey = Et - min_possible_Ex;

max_py = sqrt(2*Et);
max_y = sqrt(2*Et)/w2;
%% check for one point
numberOfenergies = 1e2;

% make the step coords
step_coords = struct();
step_coords.coords_ver = ...
    create_vertical_line(circle_center(2) + radius, [-3 circle_center(1)], 5e2);
step_coords.coords_horz = ...
    create_horizontal_line(circle_center(1) + radius, [-3 circle_center(2)], 5e2);
step_coords.coords_circ = ...
    create_cuarter_circle(circle_center, radius, 8e2);

% set all the parameters
w0 = [w1, w2];
eParams = struct();
Ex_max_ratio = [7/8, 6/8, 4/8, 2/8, 1/8 1/16 1/24 1/32 1/40];
side = ["horz", "circ", "ver"];
results = struct();
energyLevels = struct();

% calculate the energy results
for index=1:length(Ex_max_ratio)
    for jj=1:length(side)
        parameters.Ex = max_possible_Ex*Ex_max_ratio(index);
        parameters.Ey = parameters.Et - parameters.Ex;
        tmp_result = ...
            get_hit_map_yPy_version2(step_coords.('coords_' + side(jj)), parameters);
        % if energy isn't conserved
        assert(var((tmp_result.final_E_preHit - tmp_result.final_E_postHit)./Et)<1e-5,...
            "The " + 'e' + string(index) + ' ' + side(jj) + ' coords' +  " is NOT conserved");
        % insert the results into the struct
        results.('e' + string(index)).(side(jj)) = tmp_result;

        energyLevels.('e' + string(index)) = calcY_Energy_level(parameters);
    end
end

parameters.Ex = min_possible_Ex;
parameters.Ey = parameters.Et - parameters.Ex;
Max_Energy_Y_line_coords = calcY_Energy_level(parameters);

parameters.Ey = min_possible_Ey;
parameters.Ex = parameters.Et - parameters.Ey;
Min_Energy_Y_line_coords = calcY_Energy_level(parameters);

%% plot the preHit results:

legend_info = [];
fields_eLevels = fieldnames(energyLevels);


figure;
hold on
for f=1:numel(fields_eLevels)
    plot(energyLevels.(fields_eLevels{f})(:,1),energyLevels.(fields_eLevels{f})(:,2), 'k')
    legend_info = [legend_info string(fields_eLevels{f})];
end
plot(Max_Energy_Y_line_coords(:,1),Max_Energy_Y_line_coords(:,2), '--k')
plot(Min_Energy_Y_line_coords(:,1),Min_Energy_Y_line_coords(:,2), '--k')
legend_info = [legend_info "Max_Energy_Y" "Min_Energy_Y"];


for f=1:numel(fields_eLevels)
    step_part = fieldnames(results.(fields_eLevels{f}));
    for g=1:numel(step_part)
%         energy_atributs = fieldnames(results.(fields_eLevels{f}).(step_part{g}));
        plotStyle = marker(f) + marker(mod(g,length(marker)));
        plot(results.(fields_eLevels{f}).(step_part{g}).("results_preHit")(:,2),...
            results.(fields_eLevels{f}).(step_part{g}).("results_preHit")(:,4),...
            marker(f)+ color(g))
        legend_info = [legend_info,  [fields_eLevels{f}  ' '  step_part{g}  ' results_preHit' ]];
    end
end

xlabel("y")
ylabel('P_y')
title("preHit results")
legend(legend_info);
xlim([-max_y, max_y]);
ylim([-max_py, max_py]);
hold off

%% plot the postHit results:

legend_info = [];
fields_eLevels = fieldnames(energyLevels);

figure;
hold on
for f=1:numel(fields_eLevels)
    plot(energyLevels.(fields_eLevels{f})(:,1),energyLevels.(fields_eLevels{f})(:,2), 'k')
    legend_info = [legend_info string(fields_eLevels{f})];
end
plot(Max_Energy_Y_line_coords(:,1),Max_Energy_Y_line_coords(:,2), '--k')
plot(Min_Energy_Y_line_coords(:,1),Min_Energy_Y_line_coords(:,2), '--k')
legend_info = [legend_info "Max_Energy_Y" "Min_Energy_Y"];


for f=1:numel(fields_eLevels)
    step_part = fieldnames(results.(fields_eLevels{f}));
    for g=1:numel(step_part)
%         energy_atributs = fieldnames(results.(fields_eLevels{f}).(step_part{g}));
        plotStyle = marker(f) + marker(g);
        plot(results.(fields_eLevels{f}).(step_part{g}).("results_postHit")(:,2),...
            results.(fields_eLevels{f}).(step_part{g}).("results_postHit")(:,4),...
            marker(f)+ color(g))
        legend_info = [legend_info,  [fields_eLevels{f}  ' '  step_part{g}  ' results_postHit' ]];
    end
end

xlabel("y")
ylabel('P_y')
title("preHit results")
legend(legend_info);
xlim([-max_y, max_y]);
ylim([-max_py, max_py]);
hold off

%% find the index of a result:

if plot_indexes

    % for the preHit
    legend_info = [];
    fields_eLevels = fieldnames(energyLevels);
    color = ["r" "g" "b" "c" "m" "y"];
    marker = ["o" "+" "*" "x" "square" "diamond" "pentagram"];
    
    figure;
    hold on
    for f=1:numel(fields_eLevels)
        plot3(energyLevels.(fields_eLevels{f})(:,1),energyLevels.(fields_eLevels{f})(:,2), 'k')
        legend_info = [legend_info string(fields_eLevels{f})];
    end
    plot3(Max_Energy_Y_line_coords(:,1), Max_Energy_Y_line_coords(:,2), '--k')
    plot3(Min_Energy_Y_line_coords(:,1), Min_Energy_Y_line_coords(:,2), '--k')
    legend_info = [legend_info "Max_Energy_Y" "Min_Energy_Y"];
    
    
    for f=1:numel(fields_eLevels)
        step_part = fieldnames(results.(fields_eLevels{f}));
        for g=1:numel(step_part)
    %         energy_atributs = fieldnames(results.(fields_eLevels{f}).(step_part{g}));
            plotStyle = marker(f) + marker(g);
            plot3(results.(fields_eLevels{f}).(step_part{g}).("results_preHit")(:,2),...
                results.(fields_eLevels{f}).(step_part{g}).("results_preHit")(:,4),...
                marker(f)+ color(g))
            legend_info = [legend_info,  [fields_eLevels{f}  ' '  step_part{g}  ' results_preHit' ]];
        end
    end
    
    xlabel("y")
    ylabel('P_y')
    title("preHit results")
    legend(legend_info);
    xlim([-max_y, max_y]);
    ylim([-max_py, max_py]);
    hold off

    % for the postHit
    legend_info = [];
    fields_eLevels = fieldnames(energyLevels);
    color = ["r" "g" "b" "c" "m" "y"];
    marker = ["o" "+" "*" "x" "square" "diamond" "pentagram"];
    
    figure;
    hold on
    for f=1:numel(fields_eLevels)
        plot(energyLevels.(fields_eLevels{f})(:,1),energyLevels.(fields_eLevels{f})(:,2), 'k')
        legend_info = [legend_info string(fields_eLevels{f})];
    end
    plot(Max_Energy_Y_line_coords(:,1),Max_Energy_Y_line_coords(:,2), '--k')
    plot(Min_Energy_Y_line_coords(:,1),Min_Energy_Y_line_coords(:,2), '--k')
    legend_info = [legend_info "Max_Energy_Y" "Min_Energy_Y"];
    
    
    for f=1:numel(fields_eLevels)
        step_part = fieldnames(results.(fields_eLevels{f}));
        for g=1:numel(step_part)
    %         energy_atributs = fieldnames(results.(fields_eLevels{f}).(step_part{g}));
            plotStyle = marker(f) + marker(g);
            plot(results.(fields_eLevels{f}).(step_part{g}).("results_postHit")(:,2),...
                results.(fields_eLevels{f}).(step_part{g}).("results_postHit")(:,4),...
                marker(f)+ color(g))
            legend_info = [legend_info,  [fields_eLevels{f}  ' '  step_part{g}  ' results_postHit' ]];
        end
    end
    
    xlabel("y")
    ylabel('P_y')
    title("preHit results")
    legend(legend_info);
    xlim([-max_y, max_y]);
    ylim([-max_py, max_py]);
    hold off

end
%% check a specific trajactory

check_conditions_preHit = results.e2.circ.initial_conditions_preHit(265, :);
[track_check_preHit, numOfRoundHit_check_preHit, numOfHorizontalHit_check_preHit, numOfVerticalHit_check_preHit] = ...
                                        calc_track_until_px_is_zero(check_conditions_preHit, w0);

check_conditions_postHit = results.e2.circ.initial_conditions_postHit(265, :);
[track_check_postHit, numOfRoundHit_check_postHit, numOfHorizontalHit_check_postHit, numOfVerticalHit_check_postHit] = ...
                                        calc_track_until_px_is_zero(check_conditions_postHit, w0);

figure;
hold on
plot(step_coords.coords_ver(:,1),step_coords.coords_ver(:,2),'ro')
plot(step_coords.coords_circ(:,1),step_coords.coords_circ(:,2),'ro')
plot(step_coords.coords_horz(:,1),step_coords.coords_horz(:,2),'ro')
plot(track_check_preHit(:,1),track_check_preHit(:,2),'bo')
plot(track_check_postHit(:,1),track_check_postHit(:,2),'go')
legend('step', 'step', 'step', 'track_check_preHit', 'track_check_postHit')
xlim([-10,10])
ylim([-10,10])
hold off

%% plot pre to Angle-Action coordinates:

figure;
hold on

legend_info = ["max Iy"];
plot(linspace(-pi, pi, 100),(max_possible_Ey/w0(2)*ones(1, 100)),'--k')

for f=1:numel(fields_eLevels)
    step_part = fieldnames(results.(fields_eLevels{f}));
    for g=1:numel(step_part)
%         energy_atributs = fieldnames(results.(fields_eLevels{f}).(step_part{g}));
        plotStyle = marker(f) + marker(g);
        plot(results.(fields_eLevels{f}).(step_part{g}).AA_preHit_results.Thetay./pi,...
            results.(fields_eLevels{f}).(step_part{g}).AA_preHit_results.Jy,...
            marker(f)+ color(g))
        legend_info = [legend_info  string([fields_eLevels{f}  ' '  step_part{g}  ' AA_preHit' ])];
    end
end

xlabel('\theta_y');
ylabel('I_y');
title("AA preHit results")
legend(legend_info);
ylim([0, (1 + max_possible_Ey/w0(2))]);
xlim([-pi, pi]);
hold off


%% plot post to Angle-Action coordinates:

figure;
hold on

legend_info = ["max Iy"];
plot(linspace(-pi, pi, 100),(max_possible_Ey/w0(2)*ones(1, 100)),'--k')

for f=1:numel(fields_eLevels)
    step_part = fieldnames(results.(fields_eLevels{f}));
    for g=1:numel(step_part)
%         energy_atributs = fieldnames(results.(fields_eLevels{f}).(step_part{g}));
        plotStyle = marker(f) + marker(g);
        plot(results.(fields_eLevels{f}).(step_part{g}).AA_postHit_results.Thetay./pi,...
            results.(fields_eLevels{f}).(step_part{g}).AA_postHit_results.Jy,...
            marker(f)+ color(g))
        legend_info = [legend_info  string([fields_eLevels{f}  ' '  step_part{g}  ' AA_postHit' ])];
    end
end

xlabel('\theta_y');
ylabel('I_y');
title("AA postHit results")
legend(legend_info);
ylim([0, (1 + max_possible_Ey/w0(2))]);
xlim([-pi, pi]);
hold off

%% check

% syms x y px py
% 
% eqn1 = y*py - circle_center(2)*py - circle_center(1)*px - x*px == 0;
% eqn2 = 0.5*px^2 + 0.5*(w1^2)*x^2 - parameters.Ex == 0;
% eqn3 = 0.5*py^2 + 0.5*(w2^2)*y^2 - parameters.Ey == 0;
% eqn4 = x^2 - 2*x*circle_center(1) + circle_center(1)^2 ...
%             + y^2 - 2*y*circle_center(2) + circle_center(2)^2 - radius^2 == 0;
% 
% [X,Y, PX, PY] = solve([eqn1, eqn2, eqn3, eqn4], [x y px py])
% 
% 




%%

% 
% % eqn1 = y*py - circle_center(2)*py - circle_center(1)*px - x*px == 0;
% % eqn2 = 0.5*px^2 + 0.5*(w1^2)*x^2 - parameters.Ex == 0;
% % eqn3 = 0.5*py^2 + 0.5*(w2^2)*y^2 - parameters.Ey == 0;
% % eqn4 = x^2 - 2*x*circle_center(1) + circle_center(1)^2 ...
% %             + y^2 - 2*y*circle_center(2) + circle_center(2)^2 - radius^2 == 0;
% 
% 
% % q = (x y px py)
% F = @(q) [q(2)*q(4) - circle_center(2)*q(4) - circle_center(1)*q(3) - q(1)*q(3);
%          0.5*q(3)^2 + 0.5*(w1^2)*q(1)^2 - parameters.Ex;
%          0.5*q(4)^2 + 0.5*(w2^2)*q(2)^2 - parameters.Ey;
%          q(1)^2 - 2*q(1)*circle_center(1) + circle_center(1)^2 ...
%             + q(2)^2 - 2*q(2)*circle_center(2) + circle_center(2)^2 - radius^2];
% q0 = [circle_center(1) + radius*cos(pi/4);
%       circle_center(2) + radius/2*sin(pi/4);
%       0;
%       0];
% 
% options = optimoptions('fsolve','Display','iter');
% [x,fval] = fsolve(F,q0, options);
% 
% 
% %%
x = optimvar('x',4);
eqn1 = x(2)*x(4) - circle_center(2)*x(4) - circle_center(1)*x(3) + x(1)*x(3) == 0;
eqn2 = 0.5*x(3)^2 + 0.5*(w1^2)*x(1)^2 - parameters.Ex == 0;
eqn3 = 0.5*x(4)^2 + 0.5*(w2^2)*x(2)^2 - parameters.Ey == 0;
eqn4 = x(1)^2 - 2*x(1)*circle_center(1) + circle_center(1)^2 ...
        + x(2)^2 - 2*x(2)*circle_center(2) + circle_center(2)^2 - radius^2 == 0;

prob = eqnproblem;
prob.Equations.eq1 = eqn1;
prob.Equations.eq2 = eqn2;
prob.Equations.eq3 = eqn3;
prob.Equations.eq4 = eqn4;

show(prob)

x0.x = [circle_center(1) + radius*cos(pi/4);
      circle_center(2) + radius/2*sin(pi/4);
      0;
      0];

[sol,fval,exitflag] = solve(prob,x0);


%% another option

x = optimvar('x',2);
eqn1 = 0.5*(w1^2)*(x(1)^4) - circle_center(1)*(w1^2)*(x(1)^3) + ...
    (-parameters.Ex + 0.5*(circle_center(1)^2)*(w1^2))*(x(1)^2) + ...
    2*circle_center(1)*parameters.Ex*x(1) - parameters.Ex*circle_center(1)^2 ...
    - 0.5*(w2^2)*(x(2)^4) + circle_center(2)*(w2^2)*(x(2)^3) +...
    (parameters.Ey - 0.5*(circle_center(2)^2)*(w2^2))*(x(2)^2) -...
    2*circle_center(2)*parameters.Ey*x(2) + parameters.Ey*circle_center(2)^2 == 0;
eqn2 = x(1)^2 - 2*x(1)*circle_center(1) + circle_center(1)^2 ...
        + x(2)^2 - 2*x(2)*circle_center(2) + circle_center(2)^2 - radius^2 == 0;

prob = eqnproblem;
prob.Equations.eq1 = eqn1;
prob.Equations.eq2 = eqn2;

show(prob)

x0.x = [circle_center(1) + radius*cos(pi/4);
      circle_center(2) + radius*sin(pi/4);];

[sol,fval,exitflag] = solve(prob,x0);
