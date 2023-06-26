
clc;
clear;
close all;

%% Debug variables

plot_indexes = 0;
numberOfenergies = 3e2;
color = ["r" "g" "b" "c" "m" "y" "k" ];
groups_colors = ["r" "g" "b" "#808000" "y" "#FF99FF"];
marker = ["o" "+" "*" "x" "square" "diamond" "pentagram" "^" "hexagram"];
hexa_decimal_color_array = ["#B2FF66", "#80FF00", "#66CC00", "#4C9900", "#336600", ...
    "#193300", "#66FF66", "#33FF33", "#00CC00", "#009900", "#006600", "#66FFB2", ...
    "#33FF99", "#00FF80", "#00CC66", "#00994C", "#006633"];
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


%% calc all the points


% make the step coords
step_coords = struct();
step_coords.coords_ver = ...
    create_vertical_line(circle_center(2) + radius, [-3 circle_center(1)], numberOfenergies); % 5e2
step_coords.coords_horz = ...
    create_horizontal_line(circle_center(1) + radius, [-3 circle_center(2)], numberOfenergies); % 5e2
step_coords.coords_circ = ... 
    create_cuarter_circle(circle_center, radius, numberOfenergies); % 8e2

% set all the parameters
w0 = [w1, w2];
eParams = struct();
Ex_max_ratio = [24/25 15/16 7/8, 6/8, 4/8, 2/8, 1/8 1/16 1/24 1/32 1/40];
side = ["horz", "circ", "ver"];
results = struct();
energyLevels = struct();

% calculate the energy results
for index=1:length(Ex_max_ratio)
    for jj=1:length(side)
        parameters.Ex = max_possible_Ex*Ex_max_ratio(index);
        parameters.Ey = parameters.Et - parameters.Ex;
        % if strcmp(side(jj),"circ") && index==8
        %     disp(side(jj))
        % end
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

%% testbox

Cstruct = get_fully_converted_to_Ey_coords(step_coords.coords_circ, parameters);
%% the border lines

Bstruct = get_tangent_coords(step_coords.coords_circ, parameters);
circle_border_results = get_circle_border_coords(min_possible_Ey, max_possible_Ey, parameters, 2e2);
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


plot(circle_border_results.results_preHit_up(:,2),circle_border_results.results_preHit_up(:,4), '.m')
plot(circle_border_results.results_preHit_right(:,2),circle_border_results.results_preHit_right(:,4), '.m')
plot(Bstruct.results_preHit(:,2),Bstruct.results_preHit(:,4), '.c')
legend_info = [legend_info "circle_border_up" "circle_border_right" "tangent traj"];


for f=1:numel(fields_eLevels)
    step_part = fieldnames(results.(fields_eLevels{f}));
    for g=1:numel(step_part)
%         energy_atributs = fieldnames(results.(fields_eLevels{f}).(step_part{g}));
        % plotStyle = marker(f) + marker(mod(g,length(marker)));
        if strcmp(step_part{g}, 'circ')
            plot(results.(fields_eLevels{f}).(step_part{g}).("results_preHit")(:,2),...
                results.(fields_eLevels{f}).(step_part{g}).("results_preHit")(:,4),...
                'marker',marker(mod(f,length(marker))+1),'color',...
                hexa_decimal_color_array(mod(f,length(hexa_decimal_color_array))),...
                'LineStyle', 'none')
        else
            plot(results.(fields_eLevels{f}).(step_part{g}).("results_preHit")(:,2),...
                results.(fields_eLevels{f}).(step_part{g}).("results_preHit")(:,4),...
                'marker',marker(mod(f,length(marker))+1),...
                'color', color(mod(g,length(color))),'LineStyle', 'none')
        end
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


plot(circle_border_results.results_preHit_up(:,2),circle_border_results.results_preHit_up(:,4), '.m')
plot(circle_border_results.results_preHit_right(:,2),circle_border_results.results_preHit_right(:,4), '.m')
plot(Bstruct.results_preHit(:,2),Bstruct.results_preHit(:,4), '.c')
legend_info = [legend_info "circle_border_up" "circle_border_right" "tangent traj"];

for f=1:numel(fields_eLevels)
    step_part = fieldnames(results.(fields_eLevels{f}));
    for g=1:numel(step_part)
        if strcmp(step_part{g}, 'circ')
            plot(results.(fields_eLevels{f}).(step_part{g}).("results_postHit")(:,2),...
                results.(fields_eLevels{f}).(step_part{g}).("results_postHit")(:,4),...
                                'marker',marker(mod(f,length(marker))+1),'color',...
                hexa_decimal_color_array(mod(f,length(hexa_decimal_color_array))),...
                'LineStyle', 'none')           
        else
            plot(results.(fields_eLevels{f}).(step_part{g}).("results_postHit")(:,2),...
                results.(fields_eLevels{f}).(step_part{g}).("results_postHit")(:,4),...
                marker(mod(f,length(marker))+1)+ color(mod(g,length(color))))
        end
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

    figure;
    hold on
    for f=1:numel(fields_eLevels)
        plot3(energyLevels.(fields_eLevels{f})(:,1),energyLevels.(fields_eLevels{f})(:,2),...
            1:length(energyLevels.(fields_eLevels{f})(:,2)), 'k')
        legend_info = [legend_info string(fields_eLevels{f})];
    end
    plot3(Max_Energy_Y_line_coords(:,1), Max_Energy_Y_line_coords(:,2),1:length(Max_Energy_Y_line_coords(:,1)) , '--k')
    plot3(Min_Energy_Y_line_coords(:,1), Min_Energy_Y_line_coords(:,2),...
        1:length(Min_Energy_Y_line_coords(:,2)), '--k')
    legend_info = [legend_info "Max_Energy_Y" "Min_Energy_Y"];
    
    
    for f=1:numel(fields_eLevels)
        step_part = fieldnames(results.(fields_eLevels{f}));
        for g=1:numel(step_part)
             if strcmp(step_part{g}, 'circ')
                plot3(results.(fields_eLevels{f}).(step_part{g}).("results_preHit")(:,2),...
                    results.(fields_eLevels{f}).(step_part{g}).("results_preHit")(:,4),...
                    1:length(results.(fields_eLevels{f}).(step_part{g}).("results_preHit")(:,4)),...
                    'marker',marker(mod(f,length(marker))+1),'color',...
                hexa_decimal_color_array(mod(f,length(hexa_decimal_color_array))),...
                'LineStyle', 'none')
             else
                plot3(results.(fields_eLevels{f}).(step_part{g}).("results_preHit")(:,2),...
                    results.(fields_eLevels{f}).(step_part{g}).("results_preHit")(:,4),...
                    1:length(results.(fields_eLevels{f}).(step_part{g}).("results_preHit")(:,4)),...
                    marker(mod(f,length(marker))+1)+ color(mod(g,length(color))))
             end
                legend_info = [legend_info,  [fields_eLevels{f}  ' '  step_part{g}  ' results_preHit' ]];
        end
    end
    
    xlabel("y")
    ylabel('P_y')
    title("preHit results with indexes")
    legend(legend_info);
    xlim([-max_y, max_y]);
    ylim([-max_py, max_py]);
    hold off

    % for the postHit
    legend_info = [];
    fields_eLevels = fieldnames(energyLevels);
    
    figure;
    hold on
    for f=1:numel(fields_eLevels)
        plot3(energyLevels.(fields_eLevels{f})(:,1),energyLevels.(fields_eLevels{f})(:,2),...
            1:length(energyLevels.(fields_eLevels{f})(:,1)),'k')
        legend_info = [legend_info string(fields_eLevels{f})];
    end
    plot3(Max_Energy_Y_line_coords(:,1),Max_Energy_Y_line_coords(:,2),...
        1:length(Max_Energy_Y_line_coords(:,1)) ,'--k')
    plot3(Min_Energy_Y_line_coords(:,1),Min_Energy_Y_line_coords(:,2),...
        1:length(Min_Energy_Y_line_coords(:,1)),'--k')
    legend_info = [legend_info "Max_Energy_Y" "Min_Energy_Y"];
    
    
    for f=1:numel(fields_eLevels)
        step_part = fieldnames(results.(fields_eLevels{f}));
        for g=1:numel(step_part)
    %         energy_atributs = fieldnames(results.(fields_eLevels{f}).(step_part{g}));
        if strcmp(step_part{g}, 'circ')
            plot3(results.(fields_eLevels{f}).(step_part{g}).("results_postHit")(:,2),...
                results.(fields_eLevels{f}).(step_part{g}).("results_postHit")(:,4),...
                1:length(results.(fields_eLevels{f}).(step_part{g}).("results_postHit")(:,2)),...
                'marker',marker(mod(f,length(marker))+1),'color',...
                hexa_decimal_color_array(mod(f,length(hexa_decimal_color_array))),...
                'LineStyle', 'none')
        else
            plot3(results.(fields_eLevels{f}).(step_part{g}).("results_postHit")(:,2),...
                results.(fields_eLevels{f}).(step_part{g}).("results_postHit")(:,4),...
                1:length(results.(fields_eLevels{f}).(step_part{g}).("results_postHit")(:,2)),...
                marker(mod(f,length(marker))+1)+ color(mod(g,length(color))))
        end
            legend_info = [legend_info,  [fields_eLevels{f}  ' '  step_part{g}  ' results_postHit' ]];
        end
    end
    
    xlabel("y")
    ylabel('P_y')
    title("postHit results with indexes")
    legend(legend_info);
    xlim([-max_y, max_y]);
    ylim([-max_py, max_py]);
    hold off

end
%% check a specific trajactory

% check_conditions_preHit = results.e2.circ.initial_conditions_preHit(265, :);
% [track_check_preHit, numOfRoundHit_check_preHit, numOfHorizontalHit_check_preHit, numOfVerticalHit_check_preHit] = ...
%                                         calc_track_until_px_is_zero(check_conditions_preHit, w0);
% 
% check_conditions_postHit = results.e2.circ.initial_conditions_postHit(265, :);
% [track_check_postHit, numOfRoundHit_check_postHit, numOfHorizontalHit_check_postHit, numOfVerticalHit_check_postHit] = ...
%                                         calc_track_until_px_is_zero(check_conditions_postHit, w0);
% 
% figure;
% hold on
% plot(step_coords.coords_ver(:,1),step_coords.coords_ver(:,2),'ro')
% plot(step_coords.coords_circ(:,1),step_coords.coords_circ(:,2),'ro')
% plot(step_coords.coords_horz(:,1),step_coords.coords_horz(:,2),'ro')
% plot(track_check_preHit(:,1),track_check_preHit(:,2),'bo')
% plot(track_check_postHit(:,1),track_check_postHit(:,2),'go')
% legend('step', 'step', 'step', 'track_check_preHit', 'track_check_postHit')
% xlim([-10,10])
% ylim([-10,10])
% hold off

%% plot preHit to Angle-Action coordinates:

figure;
hold on

legend_info = ["max Iy"];
plot(linspace(-pi, pi, 100),(max_possible_Ey/w0(2)*ones(1, 100)),'--k')

for f=1:numel(fields_eLevels)
    step_part = fieldnames(results.(fields_eLevels{f}));
    for g=1:numel(step_part)
        if strcmp(step_part{g}, 'circ')
            plot(results.(fields_eLevels{f}).(step_part{g}).AA_preHit_results.Thetay./pi,...
                results.(fields_eLevels{f}).(step_part{g}).AA_preHit_results.Jy,...
                'marker',marker(mod(f,length(marker))+1),'color',...
                hexa_decimal_color_array(mod(f,length(hexa_decimal_color_array))),...
                'LineStyle', 'none')        
        else
            plot(results.(fields_eLevels{f}).(step_part{g}).AA_preHit_results.Thetay./pi,...
                results.(fields_eLevels{f}).(step_part{g}).AA_preHit_results.Jy,...
                 marker(mod(f,length(marker))+1)+ color(mod(g,length(color))))
        end

        legend_info = [legend_info  string([fields_eLevels{f}  ' '  step_part{g}  ' AA_preHit' ])];
    end
end

plot(Bstruct.AA_preHit_results.Thetay./pi, Bstruct.AA_preHit_results.Jy,'.c')
plot(circle_border_results.AA_preHit_results_up.Thetay./pi, circle_border_results.AA_preHit_results_up.Jy,'.m')
plot(circle_border_results.AA_preHit_results_right.Thetay./pi, circle_border_results.AA_preHit_results_right.Jy,...
    'marker','.', 'color', "#7E2F8E", 'LineStyle', 'none')
plot(Cstruct.AA_preHit_results.Thetay./pi, Cstruct.AA_preHit_results.Jy,'.k')
legend_info = [legend_info, 'tan_borders', 'circ_upper_borders',...
                'circ_right_borders', 'max_convert'];

xlabel('\theta_y');
ylabel('I_y');
title("AA preHit results")
legend(legend_info);
ylim([0, (1 + max_possible_Ey/w0(2))]);
xlim([-pi, pi]);
hold off

%% find the index of AA preHit result:

if plot_indexes
    figure;
    hold on
    
    legend_info = ["max Iy"];
    plot3(linspace(-pi, pi, 100),(max_possible_Ey/w0(2)*ones(1, 100)),...
        1:100,'--k')
    
    for f=1:numel(fields_eLevels)
        step_part = fieldnames(results.(fields_eLevels{f}));
        for g=1:numel(step_part)
            if strcmp(step_part{g}, 'circ')
                plot3(results.(fields_eLevels{f}).(step_part{g}).AA_preHit_results.Thetay./pi,...
                    results.(fields_eLevels{f}).(step_part{g}).AA_preHit_results.Jy,...
                    1:length(results.(fields_eLevels{f}).(step_part{g}).AA_preHit_results.Jy),...
                                     'marker',marker(mod(f,length(marker))+1),'color',...
                hexa_decimal_color_array(mod(f,length(hexa_decimal_color_array))),...
                'LineStyle', 'none')
            else
                plot3(results.(fields_eLevels{f}).(step_part{g}).AA_preHit_results.Thetay./pi,...
                    results.(fields_eLevels{f}).(step_part{g}).AA_preHit_results.Jy,...
                    1:length(results.(fields_eLevels{f}).(step_part{g}).AA_preHit_results.Jy),...
                     marker(mod(f,length(marker))+1)+ color(mod(g,length(color))))
            end
            legend_info = [legend_info  string([fields_eLevels{f}  ' '  step_part{g}  ' AA_preHit' ])];
        end
    end

    plot3(Bstruct.AA_preHit_results.Thetay./pi, Bstruct.AA_preHit_results.Jy,...
        1:length(Bstruct.AA_preHit_results.Jy),'.c')

    plot3(circle_border_results.AA_preHit_results_up.Thetay./pi, circle_border_results.AA_preHit_results_up.Jy,...
        1:length(circle_border_results.AA_preHit_results_up.Jy),'.m')
    plot3(circle_border_results.AA_preHit_results_right.Thetay./pi, circle_border_results.AA_preHit_results_right.Jy,...
        1:length(circle_border_results.AA_preHit_results_right.Jy),...
        'marker','.', 'color', "#7E2F8E", 'LineStyle', 'none')
    legend_info = [legend_info, 'tan_borders', 'circ_upper_borders','circ_right_borders'];
    
    xlabel('\theta_y');
    ylabel('I_y');
    title("AA preHit results with indexes")
    legend(legend_info);
    ylim([0, (1 + max_possible_Ey/w0(2))]);
    xlim([-pi, pi]);
    hold off
end
%% plot post to Angle-Action coordinates:

figure;
hold on

legend_info = ["max Iy"];
plot(linspace(-pi, pi, 100),(max_possible_Ey/w0(2)*ones(1, 100)),'--k')

plot(Bstruct.AA_preHit_results.Thetay./pi, Bstruct.AA_preHit_results.Jy,'.c')
plot(circle_border_results.AA_preHit_results_up.Thetay./pi, circle_border_results.AA_preHit_results_up.Jy,'.m')
plot(circle_border_results.AA_preHit_results_right.Thetay./pi, circle_border_results.AA_preHit_results_right.Jy,...
    'marker','.', 'color', "#7E2F8E", 'LineStyle', 'none')
plot(Cstruct.AA_postHit_results.Thetay./pi, Cstruct.AA_postHit_results.Jy,'.k')
legend_info = [legend_info, 'tan_borders', 'circ_upper_borders',...
                'circ_right_borders', 'max_convert'];

for f=1:numel(fields_eLevels)
    step_part = fieldnames(results.(fields_eLevels{f}));
    for g=1:numel(step_part)
        if strcmp(step_part{g}, 'circ')
            plot(results.(fields_eLevels{f}).(step_part{g}).AA_postHit_results.Thetay./pi,...
                results.(fields_eLevels{f}).(step_part{g}).AA_postHit_results.Jy,...
                'marker',marker(mod(f,length(marker))+1),'color',...
                hexa_decimal_color_array(mod(f,length(hexa_decimal_color_array))),...
                'LineStyle', 'none')
        else
            plot(results.(fields_eLevels{f}).(step_part{g}).AA_postHit_results.Thetay./pi,...
                results.(fields_eLevels{f}).(step_part{g}).AA_postHit_results.Jy,...
                marker(mod(f,length(marker))+1)+ color(mod(g,length(color))))
        end
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

%% plot indexes post to Angle-Action coordinates:

if plot_indexes
    figure;
    hold on
    
    legend_info = ["max Iy"];
    plot3(linspace(-pi, pi, 100),(max_possible_Ey/w0(2)*ones(1, 100)),...
        1:100,'--k')
    
    for f=1:numel(fields_eLevels)
        step_part = fieldnames(results.(fields_eLevels{f}));
        for g=1:numel(step_part)
            if strcmp(step_part{g}, 'circ')
                                plot3(results.(fields_eLevels{f}).(step_part{g}).AA_postHit_results.Thetay./pi,...
                results.(fields_eLevels{f}).(step_part{g}).AA_postHit_results.Jy,...
                1:length(results.(fields_eLevels{f}).(step_part{g}).AA_postHit_results.Thetay),...
                 'marker',marker(mod(f,length(marker))+1),'color',...
                hexa_decimal_color_array(mod(f,length(hexa_decimal_color_array))),...
                'LineStyle', 'none')
            else
                plot3(results.(fields_eLevels{f}).(step_part{g}).AA_postHit_results.Thetay./pi,...
                results.(fields_eLevels{f}).(step_part{g}).AA_postHit_results.Jy,...
                1:length(results.(fields_eLevels{f}).(step_part{g}).AA_postHit_results.Thetay)...
                ,marker(mod(f,length(marker))+1)+ color(mod(g,length(color))))
            end
            
            legend_info = [legend_info  string([fields_eLevels{f}  ' '  step_part{g}  ' AA_postHit' ])];
        end
    end
    
    xlabel('\theta_y');
    ylabel('I_y');
    title("AA  indexes postHit results")
    legend(legend_info);
    ylim([0, (1 + max_possible_Ey/w0(2))]);
    xlim([-pi, pi]);
    hold off
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

end

%% separate to groups

results_preHit_horz1 = [];
results_preHit_horz2 = [];
results_preHit_horz3 = [];
results_preHit_ver = [];
results_preHit_only_corner = [];

results_postHit_horz1 = [];
results_postHit_horz2 = [];
results_postHit_horz3 = [];
results_psotHit_ver = [];
results_postHit_only_corner = [];

results_horz1_hits = [];
results_horz2_hits = [];
results_horz3_hits = [];
results_ver_hits = [];
results_only_c_hits = [];

for f=1:numel(fields_eLevels)
    step_part = fieldnames(results.(fields_eLevels{f}));
    for g=1:numel(step_part)
        hits_horz = results.(fields_eLevels{f}).(step_part{g}).("results_NumberOfHits")(:,1);
        hits_ver = results.(fields_eLevels{f}).(step_part{g}).("results_NumberOfHits")(:,3);
        hits_corner = results.(fields_eLevels{f}).(step_part{g}).("results_NumberOfHits")(:,2);

        [horz_hits_idxes_0, ~] = find(hits_horz==0);
        [horz_hits_idxes_1, ~] = find(hits_horz==1);
        [horz_hits_idxes_2, ~] = find(hits_horz==2);
        [horz_hits_idxes_3, ~] = find(hits_horz==3);
        [ver_hits_idxes, ~] = find(hits_ver==1);
        [hits_corner, ~] = find(hits_corner>=1);

        [only_corner_indexes, ~] = intersect(horz_hits_idxes_0,hits_corner);

        % the results
        results_preHit_horz1_tmp = results.(fields_eLevels{f}).(step_part{g}).("results_preHit")(horz_hits_idxes_1, :);
        results_preHit_horz2_tmp = results.(fields_eLevels{f}).(step_part{g}).("results_preHit")(horz_hits_idxes_2, :);
        results_preHit_horz3_tmp = results.(fields_eLevels{f}).(step_part{g}).("results_preHit")(horz_hits_idxes_3, :);
        results_preHit_ver_tmp = results.(fields_eLevels{f}).(step_part{g}).("results_preHit")(ver_hits_idxes, :);
        results_preHit_only_corner_tmp = results.(fields_eLevels{f}).(step_part{g}).("results_preHit")(only_corner_indexes, :);

        results_postHit_horz1_tmp = results.(fields_eLevels{f}).(step_part{g}).("results_postHit")(horz_hits_idxes_1, :);
        results_postHit_horz2_tmp = results.(fields_eLevels{f}).(step_part{g}).("results_postHit")(horz_hits_idxes_2, :);
        results_postHit_horz3_tmp = results.(fields_eLevels{f}).(step_part{g}).("results_postHit")(horz_hits_idxes_3, :);
        results_postHit_ver_tmp = results.(fields_eLevels{f}).(step_part{g}).("results_postHit")(ver_hits_idxes, :);
        results_postHit_only_corner_tmp = results.(fields_eLevels{f}).(step_part{g}).("results_postHit")(only_corner_indexes, :);

        % hits data
        results_horz1_hits_tmp = results.(fields_eLevels{f}).(step_part{g}).("results_NumberOfHits")(horz_hits_idxes_1, :);
        results_horz2_hits_tmp = results.(fields_eLevels{f}).(step_part{g}).("results_NumberOfHits")(horz_hits_idxes_2, :);
        results_horz3_hits_tmp = results.(fields_eLevels{f}).(step_part{g}).("results_NumberOfHits")(horz_hits_idxes_3, :);
        results_ver_hits_tmp = results.(fields_eLevels{f}).(step_part{g}).("results_NumberOfHits")(ver_hits_idxes, :);
        results_only_c_hits_tmp = results.(fields_eLevels{f}).(step_part{g}).("results_NumberOfHits")(only_corner_indexes, :);

        % set into arrays
        results_horz1_hits = [results_horz1_hits; results_horz1_hits_tmp];
        results_horz2_hits = [results_horz2_hits; results_horz2_hits_tmp];
        results_horz3_hits = [results_horz3_hits; results_horz3_hits_tmp];
        results_ver_hits = [results_ver_hits; results_ver_hits_tmp];
        results_only_c_hits = [results_only_c_hits; results_only_c_hits_tmp];

        results_preHit_horz1 = [results_preHit_horz1; results_preHit_horz1_tmp];
        results_preHit_horz2 = [results_preHit_horz2; results_preHit_horz2_tmp];
        results_preHit_horz3 = [results_preHit_horz3; results_preHit_horz3_tmp];
        results_preHit_ver = [results_preHit_ver; results_preHit_ver_tmp];
        results_preHit_only_corner = [results_preHit_only_corner; results_preHit_only_corner_tmp];


        results_postHit_horz1 = [results_postHit_horz1; results_postHit_horz1_tmp];
        results_postHit_horz2 = [results_postHit_horz2; results_postHit_horz2_tmp];
        results_postHit_horz3 = [results_postHit_horz3; results_postHit_horz3_tmp];
        results_psotHit_ver = [results_psotHit_ver; results_postHit_ver_tmp];
        results_postHit_only_corner = [results_postHit_only_corner; results_postHit_only_corner_tmp];

    end
end


AA_results_preHit_horz1 = convertToAngleActionCoordsOnHarmonicHamiltonian(results_preHit_horz1, w0);
AA_results_preHit_horz2 = convertToAngleActionCoordsOnHarmonicHamiltonian(results_preHit_horz2, w0);
AA_results_preHit_horz3 = convertToAngleActionCoordsOnHarmonicHamiltonian(results_preHit_horz3, w0);
AA_results_preHit_ver = convertToAngleActionCoordsOnHarmonicHamiltonian(results_preHit_ver, w0);
AA_results_preHit_OC = convertToAngleActionCoordsOnHarmonicHamiltonian(results_preHit_only_corner, w0);

AA_results_postHit_horz1 = convertToAngleActionCoordsOnHarmonicHamiltonian(results_postHit_horz1, w0);
AA_results_postHit_horz2 = convertToAngleActionCoordsOnHarmonicHamiltonian(results_postHit_horz2, w0);
AA_results_postHit_horz3 = convertToAngleActionCoordsOnHarmonicHamiltonian(results_postHit_horz3, w0);
AA_results_psotHit_ver = convertToAngleActionCoordsOnHarmonicHamiltonian(results_psotHit_ver, w0);
AA_results_postHit_OC = convertToAngleActionCoordsOnHarmonicHamiltonian(results_postHit_only_corner, w0);

J = struct();

J.preHit.results.jr = results_preHit_ver;
J.preHit.results.j1 = results_preHit_horz1;
J.preHit.results.j2 = results_preHit_horz2;
J.preHit.results.j3 = results_preHit_horz3;
J.preHit.results.jOC = results_preHit_only_corner;

J.postHit.results.jr = results_psotHit_ver;
J.postHit.results.j1 = results_postHit_horz1;
J.postHit.results.j2 = results_postHit_horz2;
J.postHit.results.j3 = results_postHit_horz3;
J.postHit.results.jOC = results_postHit_only_corner;

J.preHit.AA.jr = AA_results_preHit_ver;
J.preHit.AA.j1 = AA_results_preHit_horz1;
J.preHit.AA.j2 = AA_results_preHit_horz2;
J.preHit.AA.j3 = AA_results_preHit_horz3;
J.preHit.AA.jOC = AA_results_preHit_OC;

J.postHit.AA.jr = AA_results_psotHit_ver;
J.postHit.AA.j1 = AA_results_postHit_horz1;
J.postHit.AA.j2 = AA_results_postHit_horz2;
J.postHit.AA.j3 = AA_results_postHit_horz3;
J.postHit.AA.jOC = AA_results_postHit_OC;

J.hits_data.jr = results_ver_hits;
J.hits_data.j1 = results_horz1_hits;
J.hits_data.j2 = results_horz2_hits;
J.hits_data.j3 = results_horz3_hits;
J.hits_data.jOC = results_only_c_hits;


%% preHit results j's groups:

legend_info = [];
fields_jGroups = fields(J.preHit.results);

figure;
hold on
for f=1:numel(fields_eLevels)
    plot(energyLevels.(fields_eLevels{f})(:,1),energyLevels.(fields_eLevels{f})(:,2), 'k')
    legend_info = [legend_info string(fields_eLevels{f})];
end
plot(Max_Energy_Y_line_coords(:,1),Max_Energy_Y_line_coords(:,2), '--k')
plot(Min_Energy_Y_line_coords(:,1),Min_Energy_Y_line_coords(:,2), '--k')
legend_info = [legend_info "Max_Energy_Y" "Min_Energy_Y"];

plot(circle_border_results.results_preHit_up(:,2),circle_border_results.results_preHit_up(:,4), '.m')
plot(circle_border_results.results_preHit_right(:,2),circle_border_results.results_preHit_right(:,4), '.m')
plot(Bstruct.results_preHit(:,2),Bstruct.results_preHit(:,4), '.c')
legend_info = [legend_info "circle_border_up" "circle_border_right" "tangent traj"];

for f=1:numel(fields_jGroups)
    if ~isempty(J.preHit.results.(fields_jGroups{f}))
        plot(J.preHit.results.(fields_jGroups{f})(:,2),...
                J.preHit.results.(fields_jGroups{f})(:,4),...
                'marker',marker(mod(f,length(marker))+1),'color',...
                groups_colors(mod(f,length(groups_colors))),...
                'LineStyle', 'none')
    
        legend_info = [legend_info  string([fields_jGroups{f}  ' result preHit' ])];
    end
end

xlabel("y")
ylabel('P_y')
title("groups preHit results")
legend(legend_info);
xlim([-max_y, max_y]);
ylim([-max_py, max_py]);
hold off

%% postHit results j's groups:

legend_info = [];
fields_jGroups = fields(J.preHit.results);

figure;
hold on
for f=1:numel(fields_eLevels)
    plot(energyLevels.(fields_eLevels{f})(:,1),energyLevels.(fields_eLevels{f})(:,2), 'k')
    legend_info = [legend_info string(fields_eLevels{f})];
end
plot(Max_Energy_Y_line_coords(:,1),Max_Energy_Y_line_coords(:,2), '--k')
plot(Min_Energy_Y_line_coords(:,1),Min_Energy_Y_line_coords(:,2), '--k')
legend_info = [legend_info "Max_Energy_Y" "Min_Energy_Y"];


plot(circle_border_results.results_preHit_up(:,2),circle_border_results.results_preHit_up(:,4), '.m')
plot(circle_border_results.results_preHit_right(:,2),circle_border_results.results_preHit_right(:,4), '.m')
plot(Bstruct.results_preHit(:,2),Bstruct.results_preHit(:,4), '.c')
legend_info = [legend_info "circle_border_up" "circle_border_right" "tangent traj"];

for f=1:numel(fields_jGroups)
    if ~isempty(J.postHit.results.(fields_jGroups{f}))
        plot(J.postHit.results.(fields_jGroups{f})(:,2),...
                J.postHit.results.(fields_jGroups{f})(:,4),...
                'marker',marker(mod(f,length(marker))+1),'color',...
                groups_colors(mod(f,length(groups_colors))),...
                'LineStyle', 'none')
    
        legend_info = [legend_info  string([fields_jGroups{f}  ' result postHit' ])];
    end
end

xlabel("y")
ylabel('P_y')
title("groups postHit results")
legend(legend_info);
xlim([-max_y, max_y]);
ylim([-max_py, max_py]);
hold off


%% preHit AA j's groups:

fields_jGroups = fields(J.preHit.AA);
figure;
hold on

legend_info = ["max Iy"];
plot(linspace(-pi, pi, 100),(max_possible_Ey/w0(2)*ones(1, 100)),'--k')

plot(Bstruct.AA_preHit_results.Thetay./pi, Bstruct.AA_preHit_results.Jy,'.c')
plot(circle_border_results.AA_preHit_results_up.Thetay./pi, circle_border_results.AA_preHit_results_up.Jy,'.m')
plot(circle_border_results.AA_preHit_results_right.Thetay./pi, circle_border_results.AA_preHit_results_right.Jy,...
    'marker','.', 'color', "#7E2F8E", 'LineStyle', 'none')
plot(Cstruct.AA_postHit_results.Thetay./pi, Cstruct.AA_postHit_results.Jy,'.k')
legend_info = [legend_info, 'tan_borders', 'circ_upper_borders',...
                'circ_right_borders', 'max_convert'];


for f=1:numel(fields_jGroups)
    if ~isempty(J.preHit.AA.(fields_jGroups{f}).Thetay)
        plot(J.preHit.AA.(fields_jGroups{f}).Thetay./pi,...
                J.preHit.AA.(fields_jGroups{f}).Jy,...
                'marker',marker(mod(f,length(marker))+1),'color',...
                groups_colors(mod(f,length(groups_colors))),...
                'LineStyle', 'none')
    
        legend_info = [legend_info  string([fields_jGroups{f}  ' AA preHit' ])];
    end
end

xlabel('\theta_y');
ylabel('I_y');
title("Groups - AA preHit results")
legend(legend_info);
ylim([0, (1 + max_possible_Ey/w0(2))]);
xlim([-pi, pi]);
hold off

%% postHit AA j's groups:

fields_jGroups = fields(J.postHit.AA);
figure;
hold on

legend_info = ["max Iy"];
plot(linspace(-pi, pi, 100),(max_possible_Ey/w0(2)*ones(1, 100)),'--k')

plot(Bstruct.AA_preHit_results.Thetay./pi, Bstruct.AA_preHit_results.Jy,'.c')
plot(circle_border_results.AA_preHit_results_up.Thetay./pi, circle_border_results.AA_preHit_results_up.Jy,'.m')
plot(circle_border_results.AA_preHit_results_right.Thetay./pi, circle_border_results.AA_preHit_results_right.Jy,...
    'marker','.', 'color', "#7E2F8E", 'LineStyle', 'none')
plot(Cstruct.AA_postHit_results.Thetay./pi, Cstruct.AA_postHit_results.Jy,'.k')
legend_info = [legend_info, 'tan_borders', 'circ_upper_borders',...
                'circ_right_borders', 'max_convert'];


for f=1:numel(fields_jGroups)
    if ~isempty(J.preHit.AA.(fields_jGroups{f}).Thetay)
        plot(J.postHit.AA.(fields_jGroups{f}).Thetay./pi,...
                J.postHit.AA.(fields_jGroups{f}).Jy,...
                'marker',marker(mod(f,length(marker))+1),'color',...
                groups_colors(mod(f,length(groups_colors))),...
                'LineStyle', 'none')
    
        legend_info = [legend_info  string([fields_jGroups{f}  ' AA postHit' ])];
    end
end

xlabel('\theta_y');
ylabel('I_y');
title("Groups - AA postHit results")
legend(legend_info);
ylim([0, (1 + max_possible_Ey/w0(2))]);
xlim([-pi, pi]);
hold off
%%

% % 
% % % eqn1 = y*py - circle_center(2)*py - circle_center(1)*px - x*px == 0;
% % % eqn2 = 0.5*px^2 + 0.5*(w1^2)*x^2 - parameters.Ex == 0;
% % % eqn3 = 0.5*py^2 + 0.5*(w2^2)*y^2 - parameters.Ey == 0;
% % % eqn4 = x^2 - 2*x*circle_center(1) + circle_center(1)^2 ...
% % %             + y^2 - 2*y*circle_center(2) + circle_center(2)^2 - radius^2 == 0;
% % 
% % 
% % % q = (x y px py)
% % F = @(q) [q(2)*q(4) - circle_center(2)*q(4) - circle_center(1)*q(3) - q(1)*q(3);
% %          0.5*q(3)^2 + 0.5*(w1^2)*q(1)^2 - parameters.Ex;
% %          0.5*q(4)^2 + 0.5*(w2^2)*q(2)^2 - parameters.Ey;
% %          q(1)^2 - 2*q(1)*circle_center(1) + circle_center(1)^2 ...
% %             + q(2)^2 - 2*q(2)*circle_center(2) + circle_center(2)^2 - radius^2];
% % q0 = [circle_center(1) + radius*cos(pi/4);
% %       circle_center(2) + radius/2*sin(pi/4);
% %       0;
% %       0];
% % 
% % options = optimoptions('fsolve','Display','iter');
% % [x,fval] = fsolve(F,q0, options);
% % 
% % 
% % %%
% x = optimvar('x',4);
% eqn1 = x(2)*x(4) - circle_center(2)*x(4) - circle_center(1)*x(3) + x(1)*x(3) == 0;
% eqn2 = 0.5*x(3)^2 + 0.5*(w1^2)*x(1)^2 - parameters.Ex == 0;
% eqn3 = 0.5*x(4)^2 + 0.5*(w2^2)*x(2)^2 - parameters.Ey == 0;
% eqn4 = x(1)^2 - 2*x(1)*circle_center(1) + circle_center(1)^2 ...
%         + x(2)^2 - 2*x(2)*circle_center(2) + circle_center(2)^2 - radius^2 == 0;
% 
% prob = eqnproblem;
% prob.Equations.eq1 = eqn1;
% prob.Equations.eq2 = eqn2;
% prob.Equations.eq3 = eqn3;
% prob.Equations.eq4 = eqn4;
% 
% show(prob)
% 
% x0.x = [circle_center(1) + radius*cos(pi/4);
%       circle_center(2) + radius/2*sin(pi/4);
%       0;
%       0];
% 
% [sol,fval,exitflag] = solve(prob,x0);
% 
% 
% %% another option
% 
% x = optimvar('x',2);
% eqn1 = 0.5*(w1^2)*(x(1)^4) - circle_center(1)*(w1^2)*(x(1)^3) + ...
%     (-parameters.Ex + 0.5*(circle_center(1)^2)*(w1^2))*(x(1)^2) + ...
%     2*circle_center(1)*parameters.Ex*x(1) - parameters.Ex*circle_center(1)^2 ...
%     - 0.5*(w2^2)*(x(2)^4) + circle_center(2)*(w2^2)*(x(2)^3) +...
%     (parameters.Ey - 0.5*(circle_center(2)^2)*(w2^2))*(x(2)^2) -...
%     2*circle_center(2)*parameters.Ey*x(2) + parameters.Ey*circle_center(2)^2 == 0;
% eqn2 = x(1)^2 - 2*x(1)*circle_center(1) + circle_center(1)^2 ...
%         + x(2)^2 - 2*x(2)*circle_center(2) + circle_center(2)^2 - radius^2 == 0;
% 
% prob = eqnproblem;
% prob.Equations.eq1 = eqn1;
% prob.Equations.eq2 = eqn2;
% 
% show(prob)
% 
% x0.x = [circle_center(1) + radius*cos(pi/4);
%       circle_center(2) + radius*sin(pi/4);];
% 
% [sol,fval,exitflag] = solve(prob,x0);
