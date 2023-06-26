clc;
clear;
close all;
%%
set(0,'units','pixels')  ;
%Obtains this pixel information
Pix_SS = get(0,'screensize');

parameters = params();
% parameters.radius = 0.05;
% parameters.w2 = sqrt(2);
numberOfenergies = 2e2;
debug = 1;
plot_indexes = 0;


%% effort to save time:
use_saved_res = 0;

if use_saved_res
    results_struct = load("saved_results\numberOfenergies_2e2_radius_0p1_w1_1_w2_sqrt2.mat");
else
    results_struct = get_setup_results(parameters, numberOfenergies, debug, plot_indexes);
    save('saved_results\numberOfenergies_2e2_radius_0p1_w1_1_w2_sqrt2.mat', '-struct','results_struct')
end
% saved_results = dir('saved_results\');

%% show a energy spread CIRC:
set(groot,'defaultLineMarkerSize',4);

figure;
hold on

AA_coords_up = results_struct.results_by_energy.e6.circ;
AA_coords_down = results_struct.results_by_energy.e7.circ;
% 
% for the upper part
relevant_indices_up = 1:length(AA_coords_up.AA_preHit_results.Thetay);
% relevant_indices_up = find(abs(AA_coords_up.AA_preHit_results.Thetay./pi - 0.23) < 0.1);
% relevant_indices_up = find(abs(AA_coords_up.AA_preHit_results.Thetay./pi - 0.45) < 0.1);
% relevant_indices_up = find(abs(AA_coords_up.AA_preHit_results.Thetay./pi - (-0.7)) < 0.1);


% for the lower part
relevant_indices_low = 1:length(AA_coords_down.AA_preHit_results.Thetay);
% relevant_indices_low = find(abs(AA_coords_down.AA_preHit_results.Thetay./pi - 0.1) < 0.07);
% relevant_indices_low = find(abs(AA_coords_down.AA_preHit_results.Thetay./pi - 0.275) < 0.1);
% relevant_indices_low = find(abs(AA_coords_down.AA_preHit_results.Thetay./pi - (-0.625)) < 0.1);

legend_info = ["max Iy"];
plot(linspace(-pi, pi, 100),((results_struct.max_possible_Ey/parameters.w2)*ones(1, 100)),'--k')

% for the upper energy
plot(AA_coords_up.AA_preHit_results.Thetay(relevant_indices_up)./pi,...
    AA_coords_up.AA_preHit_results.Jy(relevant_indices_up),...
    'go')
plot(AA_coords_up.AA_postHit_results.Thetay(relevant_indices_up)./pi,...
    AA_coords_up.AA_postHit_results.Jy(relevant_indices_up),...
    'r*')
legend_info = [legend_info  , 'circ side AA-preHit upper-energy', 'circ side AA-postHit upper-energy'];

% for the lower energy
plot(AA_coords_down.AA_preHit_results.Thetay(relevant_indices_low)./pi,...
    AA_coords_down.AA_preHit_results.Jy(relevant_indices_low),...
    'yo')
plot(AA_coords_down.AA_postHit_results.Thetay(relevant_indices_low)./pi,...
    AA_coords_down.AA_postHit_results.Jy(relevant_indices_low),...
    'b*')
legend_info = [legend_info  , 'circ side AA-preHit  lower energy', 'circ side AA-postHit lower energy'];

% the borders
plot(results_struct.tangent_coords.AA_preHit_results.Thetay./pi, results_struct.tangent_coords.AA_preHit_results.Jy,'.c')
plot(results_struct.circle_border_results.AA_preHit_results_up.Thetay./pi, results_struct.circle_border_results.AA_preHit_results_up.Jy,'.m')
plot(results_struct.circle_border_results.AA_preHit_results_right.Thetay./pi, results_struct.circle_border_results.AA_preHit_results_right.Jy,...
    'marker','.', 'color', "#66B2FF", 'LineStyle', 'none')
plot(results_struct.fully_converted.AA_preHit_results.Thetay./pi, results_struct.fully_converted.AA_preHit_results.Jy,'.k')

plot(results_struct.fully_converted.AA_postHit_results.Thetay./pi, results_struct.fully_converted.AA_postHit_results.Jy, 'Dc')
plot(results_struct.circle_border_results.AA_postHit_results_right.Thetay./pi, results_struct.circle_border_results.AA_postHit_results_right.Jy,...
    'marker','diamond', 'color', "#66B2FF", 'LineStyle', 'none')
plot(results_struct.circle_border_results.AA_postHit_results_up.Thetay./pi, results_struct.circle_border_results.AA_postHit_results_up.Jy, 'Dm')




legend_info = [legend_info, 'tan_borders', 'circ_upper_borders',...
                'circ_right_borders', 'max_convert', 'post_fullyconverted'...
                ,'post_circ_border_right', 'post_circ_border_up'];

plot(results_struct.circle_border_results.AA_postHit_results_up.Thetay./pi, results_struct.circle_border_results.AA_postHit_results_up.Jy,'.b')
plot(results_struct.circle_border_results.AA_postHit_results_right.Thetay./pi, results_struct.circle_border_results.AA_postHit_results_right.Jy,...
    'marker','.', 'color', "#9933FF", 'LineStyle', 'none')
legend_info = [legend_info, 'circ_post_upper_borders', 'circ_post_right_borders'];


xlabel('\theta_y');
ylabel('I_y');
title("AA map of two energies - right group")
legend(legend_info);
ylim([0, (1 + results_struct.max_possible_Ey/parameters.w2)]);
xlim([-pi, pi]);
hold off

%% plot all map graph:

intervals_centers = [-0.455, -0.23, 0.696];
title_names = ["left_group", "middle_group", "right_group"];
figure_size = [450 350];
graph_location_size = [50 (Pix_SS(4) - figure_size(1)) figure_size];

figures_struct = struct();

% TODO: to add a function that can show plots on the same screen

for ii = 1:length(intervals_centers)
    figures_struct.('graph_' + title_names(ii)) = figure('Position',graph_location_size);
    hold on
    relevant_indices_up = find(abs(AA_coords_up.AA_preHit_results.Thetay./pi - intervals_centers(ii)) < 0.02);
    relevant_indices_low = find(abs(AA_coords_down.AA_preHit_results.Thetay./pi - intervals_centers(ii)) < 0.02);

    legend_info = ["max Iy"];
    plot(linspace(-pi, pi, 100),((results_struct.max_possible_Ey/parameters.w2)*ones(1, 100)),'--k')
    
    % for the upper energy
    plot(AA_coords_up.AA_preHit_results.Thetay(relevant_indices_up)./pi,...
        AA_coords_up.AA_preHit_results.Jy(relevant_indices_up),...
        'go')
    plot(AA_coords_up.AA_postHit_results.Thetay(relevant_indices_up)./pi,...
        AA_coords_up.AA_postHit_results.Jy(relevant_indices_up),...
        'r*')
    legend_info = [legend_info  , 'circ side AA-preHit upper-energy', 'circ side AA-postHit upper-energy'];
    
    % for the lower energy
    plot(AA_coords_down.AA_preHit_results.Thetay(relevant_indices_low)./pi,...
        AA_coords_down.AA_preHit_results.Jy(relevant_indices_low),...
        'yo')
    plot(AA_coords_down.AA_postHit_results.Thetay(relevant_indices_low)./pi,...
        AA_coords_down.AA_postHit_results.Jy(relevant_indices_low),...
        'b*')
    legend_info = [legend_info  , 'circ side AA-preHit  lower energy', 'circ side AA-postHit lower energy'];
    
    % the borders
    plot(results_struct.tangent_coords.AA_preHit_results.Thetay./pi, results_struct.tangent_coords.AA_preHit_results.Jy,'.c')
    plot(results_struct.circle_border_results.AA_preHit_results_up.Thetay./pi, results_struct.circle_border_results.AA_preHit_results_up.Jy,'.m')
    plot(results_struct.circle_border_results.AA_preHit_results_right.Thetay./pi, results_struct.circle_border_results.AA_preHit_results_right.Jy,...
        'marker','.', 'color', "#66B2FF", 'LineStyle', 'none')
    plot(results_struct.fully_converted.AA_preHit_results.Thetay./pi, results_struct.fully_converted.AA_preHit_results.Jy,'.k')
    
    plot(results_struct.Xfully_converted_toX.AA_preHit_results.Thetay./pi, results_struct.Xfully_converted_toX.AA_preHit_results.Jy, 'Dc')
    
    legend_info = [legend_info, 'tan_borders', 'circ_upper_borders',...
                    'circ_right_borders', 'max_convert to Y', 'max_convert to X'];
    
    plot(results_struct.circle_border_results.AA_postHit_results_up.Thetay./pi, results_struct.circle_border_results.AA_postHit_results_up.Jy,'.b')
    plot(results_struct.circle_border_results.AA_postHit_results_right.Thetay./pi, results_struct.circle_border_results.AA_postHit_results_right.Jy,...
        'marker','.', 'color', "#9933FF", 'LineStyle', 'none')
    legend_info = [legend_info, 'circ_post_upper_borders', 'circ_post_right_borders'];
    
    
    xlabel('\theta_y');
    ylabel('I_y');
    title(title_names(ii))
    legend(legend_info);
    ylim([0, (1 + results_struct.max_possible_Ey/parameters.w2)]);
    xlim([-pi, pi]);
    hold off
    
    movegui(figures_struct.('graph_' + title_names(ii)));
    graph_location_size(1) = graph_location_size(1) + 1.2*figure_size(2);
end

%% show a energy spread VRTICAL:

figure;
hold on

AA_coords_up = results_struct.results_by_energy.e6.ver;
% relevant_indices = 1:length(AA_coords.AA_preHit_results.Thetay);
relevant_indices_up = find(abs(AA_coords_up.AA_preHit_results.Thetay./pi - 0.23) < 0.1);
% relevant_indices = find(abs(AA_coords.AA_preHit_results.Thetay./pi - 0.45) < 0.1);
% relevant_indices = find(abs(AA_coords.AA_preHit_results.Thetay./pi - (-0.7)) < 0.1);

legend_info = ["max Iy"];
plot(linspace(-pi, pi, 100),((results_struct.max_possible_Ey/parameters.w2)*ones(1, 100)),'--k')
plot(AA_coords_up.AA_preHit_results.Thetay(relevant_indices_up)./pi,...
    AA_coords_up.AA_preHit_results.Jy(relevant_indices_up),...
    'go')
plot(AA_coords_up.AA_postHit_results.Thetay(relevant_indices_up)./pi,...
    AA_coords_up.AA_postHit_results.Jy(relevant_indices_up),...
    'r*')
legend_info = [legend_info  , 'circ side AA-preHit', 'circ side AA-postHit'];

plot(results_struct.tangent_coords.AA_preHit_results.Thetay./pi, results_struct.tangent_coords.AA_preHit_results.Jy,'.c')
plot(results_struct.circle_border_results.AA_preHit_results_up.Thetay./pi, results_struct.circle_border_results.AA_preHit_results_up.Jy,'.m')
plot(results_struct.circle_border_results.AA_preHit_results_right.Thetay./pi, results_struct.circle_border_results.AA_preHit_results_right.Jy,...
    'marker','.', 'color', "#66B2FF", 'LineStyle', 'none')
plot(results_struct.fully_converted.AA_preHit_results.Thetay./pi, results_struct.fully_converted.AA_preHit_results.Jy,'.k')
legend_info = [legend_info, 'tan_borders', 'circ_upper_borders',...
                'circ_right_borders', 'max_convert'];

% plot(results_struct.circle_border_results.AA_postHit_results_up.Thetay./pi, results_struct.circle_border_results.AA_postHit_results_up.Jy,'.b')
% plot(results_struct.circle_border_results.AA_postHit_results_right.Thetay./pi, results_struct.circle_border_results.AA_postHit_results_right.Jy,...
%     'marker','.', 'color', "#9933FF", 'LineStyle', 'none')
% legend_info = [legend_info, 'circ_post_upper_borders', 'circ_post_right_borders'];


xlabel('\theta_y');
ylabel('I_y');
title("AA preHit results")
legend(legend_info);
ylim([0, (1 + results_struct.max_possible_Ey/parameters.w2)]);
xlim([-pi, pi]);
hold off

%% show a energy spread HORIZONTAL:

figure;
hold on

AA_coords_up = results_struct.results_by_energy.e6.horz;
% relevant_indices = 1:length(AA_coords.AA_preHit_results.Thetay);
% relevant_indices = find(abs(AA_coords.AA_preHit_results.Thetay./pi - 0.23) < 0.1);
% relevant_indices = find(abs(AA_coords.AA_preHit_results.Thetay./pi - 0.45) < 0.1);
relevant_indices_up = find(abs(AA_coords_up.AA_preHit_results.Thetay./pi - (-0.7)) < 0.1);

legend_info = ["max Iy"];
plot(linspace(-pi, pi, 100),((results_struct.max_possible_Ey/parameters.w2)*ones(1, 100)),'--k')
plot(AA_coords_up.AA_preHit_results.Thetay(relevant_indices_up)./pi,...
    AA_coords_up.AA_preHit_results.Jy(relevant_indices_up),...
    'g+')
plot(AA_coords_up.AA_postHit_results.Thetay(relevant_indices_up)./pi,...
    AA_coords_up.AA_postHit_results.Jy(relevant_indices_up),...
    'ro')
legend_info = [legend_info  , 'circ side AA_preHit', 'circ side AA_postHit'];

plot(results_struct.tangent_coords.AA_preHit_results.Thetay./pi, results_struct.tangent_coords.AA_preHit_results.Jy,'.c')
plot(results_struct.circle_border_results.AA_preHit_results_up.Thetay./pi, results_struct.circle_border_results.AA_preHit_results_up.Jy,'.m')
plot(results_struct.circle_border_results.AA_preHit_results_right.Thetay./pi, results_struct.circle_border_results.AA_preHit_results_right.Jy,...
    'marker','.', 'color', "#66B2FF", 'LineStyle', 'none')
plot(results_struct.fully_converted.AA_preHit_results.Thetay./pi, results_struct.fully_converted.AA_preHit_results.Jy,'.k')
legend_info = [legend_info, 'tan_borders', 'circ_upper_borders',...
                'circ_right_borders', 'max_convert'];

% plot(results_struct.circle_border_results.AA_postHit_results_up.Thetay./pi, results_struct.circle_border_results.AA_postHit_results_up.Jy,'.b')
% plot(results_struct.circle_border_results.AA_postHit_results_right.Thetay./pi, results_struct.circle_border_results.AA_postHit_results_right.Jy,...
%     'marker','.', 'color', "#9933FF", 'LineStyle', 'none')
% legend_info = [legend_info, 'circ_post_upper_borders', 'circ_post_right_borders'];


xlabel('\theta_y');
ylabel('I_y');
title("AA preHit results")
legend(legend_info);
ylim([0, (1 + results_struct.max_possible_Ey/parameters.w2)]);
xlim([-pi, pi]);
hold off
%% show a energy spread:

figure;
hold on

AA_coords_up = results_struct.results_by_energy.e3.circ;
% relevant_indices = 1:length(AA_coords.AA_preHit_results.Thetay);
% relevant_indices = find(abs(AA_coords.AA_preHit_results.Thetay./pi - 0.15) < 0.1);
% relevant_indices = find(abs(AA_coords.AA_preHit_results.Thetay./pi - 0.35) < 0.1);
relevant_indices_up = find(abs(AA_coords_up.AA_preHit_results.Thetay./pi - (-0.65)) < 0.1);

legend_info = ["max Iy"];
plot(linspace(-pi, pi, 100),((results_struct.max_possible_Ey/parameters.w2)*ones(1, 100)),'--k')
plot(AA_coords_up.AA_preHit_results.Thetay(relevant_indices_up)./pi,...
    AA_coords_up.AA_preHit_results.Jy(relevant_indices_up),...
    'go')
plot(AA_coords_up.AA_postHit_results.Thetay(relevant_indices_up)./pi,...
    AA_coords_up.AA_postHit_results.Jy(relevant_indices_up),...
    'ro')
legend_info = [legend_info  , 'circ side AA_preHit', 'circ side AA_postHit'];

plot(results_struct.tangent_coords.AA_preHit_results.Thetay./pi, results_struct.tangent_coords.AA_preHit_results.Jy,'.c')
plot(results_struct.circle_border_results.AA_preHit_results_up.Thetay./pi, results_struct.circle_border_results.AA_preHit_results_up.Jy,'.m')
plot(results_struct.circle_border_results.AA_preHit_results_right.Thetay./pi, results_struct.circle_border_results.AA_preHit_results_right.Jy,...
    'marker','.', 'color', "#9933FF", 'LineStyle', 'none')
plot(results_struct.fully_converted.AA_preHit_results.Thetay./pi, results_struct.fully_converted.AA_preHit_results.Jy,'.k')
legend_info = [legend_info, 'tan_borders', 'circ_upper_borders',...
                'circ_right_borders', 'max_convert'];

xlabel('\theta_y');
ylabel('I_y');
title("AA preHit results")
legend(legend_info);
ylim([0, (1 + results_struct.max_possible_Ey/parameters.w2)]);
xlim([-pi, pi]);
hold off

%%

figure;
hold on;
% plot(results_struct.fully_converted.AA_preHit_results.Thetay./pi, results_struct.fully_converted.AA_preHit_results.Jy, 'ro')
plot(results_struct.fully_converted.AA_postHit_results.Thetay./pi, results_struct.fully_converted.AA_postHit_results.Jy, 'r+')
% plot(results_struct.circle_border_results.AA_preHit_results_right.Thetay./pi, results_struct.circle_border_results.AA_preHit_results_right.Jy, 'go')
plot(results_struct.circle_border_results.AA_postHit_results_right.Thetay./pi, results_struct.circle_border_results.AA_postHit_results_right.Jy, 'g.')
% plot(results_struct.circle_border_results.AA_preHit_results_up.Thetay./pi, results_struct.circle_border_results.AA_preHit_results_up.Jy, 'co')
plot(results_struct.circle_border_results.AA_postHit_results_up.Thetay./pi, results_struct.circle_border_results.AA_postHit_results_up.Jy, 'c.')
% plot(results_struct.tangent_coords.AA_preHit_results.Thetay./pi, results_struct.tangent_coords.AA_preHit_results.Jy, 'mo')
plot(results_struct.tangent_coords.AA_postHit_results.Thetay./pi, results_struct.tangent_coords.AA_postHit_results.Jy, 'm+')
% plot(results_struct.Xfully_converted_toX.AA_preHit_results.Thetay./pi, results_struct.Xfully_converted_toX.AA_preHit_results.Jy, 'ko')
plot(results_struct.Xfully_converted_toX.AA_postHit_results.Thetay./pi, results_struct.Xfully_converted_toX.AA_postHit_results.Jy, 'k+')
% plot(results_struct.tangent_coords.AA_postHit_results.Thetay./pi, results_struct.tangent_coords.AA_postHit_results.Jy, 'm+')

%% try so sort the data:
