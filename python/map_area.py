import numpy as np
import matplotlib.pyplot as plt
from straight_space_animation import twoD_Harmonic_separable_hamiltonian, dynamics, Step, Trimmed_step, circled_edge
from utils import save_results_to_date_file

h_initial_conditions = {
    'current_q1': 400,
    'current_q2': 400,
    'w1': 1,
    'w2': 2,
    'res': (800, 800)
}


def calc_map_of_area(initial_conditions, step, q1_left, q1_right, q2_down, q2_up, **kwargs):
    debug = kwargs.setdefault('debug', 0)
    show_verticle_coordinates = kwargs.setdefault('show_verticle_coordinates', 0)
    show_slope_coordinates = kwargs.setdefault('show_slope_coordinates', 0)
    show_horizontal_coordinates = kwargs.setdefault('show_horizontal_coordinates', 0)
    total_iterations = (q1_right - q1_left) * (q2_up - q2_down)
    iterator = 0
    coords_firstHit = list()
    coords_SecondtHit = list()
    for i in range(q1_left, q1_right):
        for j in range(q2_down, q2_up):
            if step.inside(i, j):
                continue
            stop_iteration = 0
            num_of_hits = 0
            # print(i)
            initial_conditions['current_q1'] = i
            initial_conditions['current_q2'] = j
            h = twoD_Harmonic_separable_hamiltonian(**initial_conditions)
            iterator += 1

            if iterator % 1000 == 0:
                print(f'precentage of process {(iterator / total_iterations) * 100:.2f}%')

            while True:
                dynamics(hamiltonian=h, obstacle=step)

                if stop_iteration > 5000:
                    break
                if step.hit_flag:
                    if num_of_hits == 0:
                        J_2 = 0.5 * h.w2 * (initial_conditions['current_q2'] ** 2)
                        if initial_conditions['current_q2'] == h.previous_q2:
                            theta_wall = np.pi
                        else:
                            theta_wall = (np.pi / 2) - np.arctan(
                                h.previous_q2 / np.sqrt(initial_conditions['current_q2'] ** 2 - h.previous_q2 ** 2))
                            if h.previous_p2 < 0:
                                theta_wall = 2 * np.pi - theta_wall

                        coords_firstHit.append((J_2, theta_wall, h.previous_q1, h.previous_q2, (i, j), iterator))
                        num_of_hits += 1
                    else:
                        E_1 = 0.5 * (h.current_p1 ** 2) + 0.5 * ((h.w1 * h.current_q1) ** 2)
                        E_2 = 0.5 * (h.current_p2 ** 2) + 0.5 * ((h.w2 * h.current_q2) ** 2)
                        J_2 = E_2 / h.w2
                        q1_max = np.sign(h.current_p1) * np.sqrt(2 * E_1) / h.w1
                        q2_max = np.sign(h.current_p2) * np.sqrt(2 * E_2) / h.w2
                        if q2_max == h.previous_q2:
                            theta_wall = np.pi
                        else:
                            theta_wall = (np.pi / 2) - np.arctan(
                                h.previous_q2 / np.sqrt(q2_max ** 2 - h.previous_q2 ** 2))
                            if h.previous_p2 < 0:
                                theta_wall = 2 * np.pi - theta_wall
                        coords_SecondtHit.append(
                            (J_2, theta_wall, h.previous_q1, h.previous_q2, (q1_max, q2_max), iterator))
                        break
                stop_iteration += 1

    if debug:
        x1 = [val[2] for val in coords_firstHit]
        y1 = [val[3] for val in coords_firstHit]

        x2 = [val[2] for val in coords_SecondtHit]
        y2 = [val[3] for val in coords_SecondtHit]

        fig = plt.figure(figsize=(16, 5))  # create a figure
        ax1 = fig.add_subplot(1, 2, 1)  # create a subplot of certain size
        ax1.scatter(x1, y1, c='r')
        # for i, txt in enumerate(x1):
        #     ax1.annotate(str(coords_firstHit[i][4]), (x1[i], y1[i]))
        ax1.set_title("first hit")

        ax2 = fig.add_subplot(1, 2, 2)
        ax2.scatter(x2, y2, c='b')
        # for i, txt in enumerate(x2):
        #     ax2.annotate(str(coords_SecondtHit[i][4]), (x2[i], y2[i]))
        ax2.set_title("second hit")
        plt.show()

    J_vertical1 = [val[0] for val in coords_firstHit if val[2] == step.x_right]
    theta_vertical1 = [val[1] for val in coords_firstHit if val[2] == step.x_right]

    J_horizontal1 = [val[0] for val in coords_firstHit if val[3] == step.y_up]
    theta_horizontal1 = [val[1] for val in coords_firstHit if val[3] == step.y_up]

    J_vertical2 = [val[0] for val in coords_SecondtHit if val[2] == step.x_right]
    theta_vertical2 = [val[1] for val in coords_SecondtHit if val[2] == step.x_right]

    J_horizontal2 = [val[0] for val in coords_SecondtHit if val[3] == step.y_up]
    theta_horizontal2 = [val[1] for val in coords_SecondtHit if val[3] == step.y_up]

    J_middle_rigion1 = [val[0] for val in coords_firstHit if ((val[3] != step.y_up) and (val[2] != step.x_right))]
    theta_middle_rigion1 = [val[1] for val in coords_firstHit if ((val[3] != step.y_up) and (val[2] != step.x_right))]

    J_middle_rigion2 = [val[0] for val in coords_SecondtHit if ((val[3] != step.y_up) and (val[2] != step.x_right))]
    theta_middle_rigion2 = [val[1] for val in coords_SecondtHit if ((val[3] != step.y_up) and (val[2] != step.x_right))]

    fig = plt.figure(figsize=(16, 5))  # create a figure
    ax1 = fig.add_subplot(1, 2, 1)  # create a subplot of certain size

    ax1.scatter(theta_vertical1, J_vertical1, c='r')
    if show_verticle_coordinates:
        txt = [str(val[4]) for val in coords_firstHit if val[2] == step.x_right]
        for i, txt_i in enumerate(txt):
            ax1.annotate(txt_i, (theta_vertical1[i], J_vertical1[i]))

    ax1.scatter(theta_horizontal1, J_horizontal1, c='b')
    if show_horizontal_coordinates:
        txt = [str(val[4]) for val in coords_firstHit if val[3] == step.y_up]
        for i, txt_i in enumerate(txt):
            ax1.annotate(txt_i, (theta_horizontal1[i], J_horizontal1[i]))

    ax1.scatter(theta_middle_rigion1, J_middle_rigion1, c='g')
    if show_slope_coordinates:
        txt = [str(val[4]) for val in coords_firstHit if ((val[3] != step.y_up) and (val[2] != step.x_right))]
        for i, txt_i in enumerate(txt):
            ax1.annotate(txt_i, (theta_middle_rigion1[i], J_middle_rigion1[i]))

    ax1.set_title("first hit")
    ax1.legend(['vertical side', 'horizontal side', 'slope side'])
    ax1.set_xlabel('theta1')
    ax1.set_ylabel('J1')

    ax2 = fig.add_subplot(1, 2, 2)

    ax2.scatter(theta_vertical2, J_vertical2, c='r')
    if show_verticle_coordinates:
        txt = [str(val[4]) for val in coords_SecondtHit if val[2] == step.x_right]
        for i, txt_i in enumerate(txt):
            ax2.annotate(txt_i, (theta_vertical2[i], J_vertical2[i]))

    ax2.scatter(theta_horizontal2, J_horizontal2, c='b')
    if show_horizontal_coordinates:
        txt = [str(val[4]) for val in coords_SecondtHit if val[3] == step.y_up]
        for i, txt_i in enumerate(txt):
            ax2.annotate(txt_i, (theta_horizontal2[i], J_horizontal2[i]))

    ax2.scatter(theta_middle_rigion2, J_middle_rigion2, c='g')
    if show_slope_coordinates:
        txt = [str(val[4]) for val in coords_SecondtHit if ((val[3] != step.y_up) and (val[2] != step.x_right))]
        for i, txt_i in enumerate(txt):
            ax2.annotate(txt_i, (theta_middle_rigion2[i], J_middle_rigion2[i]))

    ax2.set_title("second hit")
    ax2.legend(['vertical side', 'horizontal side', 'slope side'])
    ax2.set_xlabel('theta2')
    ax2.set_ylabel('J2')
    plt.show()

    return coords_firstHit, coords_SecondtHit


# calc_map_of_area(initial_conditions=h_initial_conditions, step=Step(), q1_left=250, q1_right=300, q2_down=250,
#                  q2_up=300, debug=0)
# calc_map_of_area(initial_conditions=h_initial_conditions, step=Step(), q1_left=250, q1_right=300, q2_down=-400,
#                  q2_up=-300, debug=0)
# Trimmed step
# calc_map_of_area(initial_conditions=h_initial_conditions, step=Trimmed_step(), q1_left=250, q1_right=300, q2_down=-400,
#                  q2_up=-350, debug=0)


show_lables = {
    'show_verticle_coordinates': 0,
    'show_slope_coordinates': 0,
    'show_horizontal_coordinates': 0,
    'debug': 0
}
# ROUNDED step
# calc_map_of_area(initial_conditions=h_initial_conditions, step=circled_edge(top_right_radius=80), q1_left=200,
#                  q1_right=250, q2_down=-300,
#                  q2_up=-200, **show_lables)

coords_firstHit_rounded_r40, coords_SecondHit_rounded_r40 = calc_map_of_area(initial_conditions=h_initial_conditions,
                                                                             step=circled_edge(top_right_radius=40),
                                                                             q1_left=-400,
                                                                             q1_right=400, q2_down=-400,
                                                                             q2_up=400, **show_lables)
save_results_to_date_file({'first_hit': coords_firstHit_rounded_r40, 'second_hit': coords_SecondHit_rounded_r40},
                          file_name='all_board_r40')

#
# coords_firstHit_rounded_r40, coords_SecondHit_rounded_r40 = calc_map_of_area(initial_conditions=h_initial_conditions,
#                                                                              step=circled_edge(top_right_radius=40),
#                                                                              q1_left=50,
#                                                                              q1_right=400, q2_down=-400,
#                                                                              q2_up=400, **show_lables)
# save_results_to_date_file({'first_hit': coords_firstHit_rounded_r40, 'second_hit': coords_SecondHit_rounded_r40},
#                           file_name='right_side_r40')
#
# coords_firstHit_rounded_r40, coords_SecondHit_rounded_r40 = calc_map_of_area(initial_conditions=h_initial_conditions,
#                                                                              step=circled_edge(top_right_radius=40),
#                                                                              q1_left=-400,
#                                                                              q1_right=-50, q2_down=50,
#                                                                              q2_up=400, **show_lables)
# save_results_to_date_file({'first_hit': coords_firstHit_rounded_r40, 'second_hit': coords_SecondHit_rounded_r40},
#                           file_name='left_upper_side_r40')
