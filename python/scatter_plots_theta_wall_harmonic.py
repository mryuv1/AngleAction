import numpy as np
import matplotlib.pyplot as plt
from straight_space_animation import twoD_Harmonic_separable_hamiltonian, dynamics, Step, Trimmed_step, circled_edge

h_initial_conditions = {
    'current_q1': 55,
    'current_q2': 300,
    'w1': 1,
    'w2': 2,
    'res': (800, 800)
}

clac_step = 1
calc_trim = 1
calc_circ = 1

if clac_step:
    coords_firstHit = list()
    coords_SecondtHit = list()
    step = Step()
    for i in range(-400, 401):
        stop_iteration = 0
        num_of_hits = 0
        # print(i)
        h_initial_conditions['current_q2'] = i
        h = twoD_Harmonic_separable_hamiltonian(**h_initial_conditions)

        while True:
            dynamics(hamiltonian=h, obstacle=step)

            if stop_iteration > 5000:
                break
            if step.hit_flag:
                if num_of_hits == 0:
                    J_2 = 0.5 * h.w2 * (h_initial_conditions['current_q2'] ** 2)
                    if h_initial_conditions['current_q2'] == h.previous_q2:
                        theta_wall = np.pi
                    else:
                        theta_wall = (np.pi / 2) - np.arctan(
                            h.previous_q2 / np.sqrt(h_initial_conditions['current_q2'] ** 2 - h.previous_q2 ** 2))
                        if h.previous_p2 < 0:
                            theta_wall = 2 * np.pi - theta_wall
                    coords_firstHit.append((J_2, theta_wall, h.previous_q1, h.previous_q2, i))
                    num_of_hits += 1
                else:
                    J_2 = 0.5 * h.w2 * (h_initial_conditions['current_q2'] ** 2)
                    if h_initial_conditions['current_q2'] == h.previous_q2:
                        theta_wall = np.pi
                    else:
                        theta_wall = (np.pi / 2) - np.arctan(
                            h.previous_q2 / np.sqrt(h_initial_conditions['current_q2'] ** 2 - h.previous_q2 ** 2))
                        if h.previous_p2 < 0:
                            theta_wall = 2 * np.pi - theta_wall
                    coords_SecondtHit.append((J_2, theta_wall, h.previous_q1, h.previous_q2, i))
                    break
            stop_iteration += 1

    J_vertical1 = [val[0] for val in coords_firstHit if val[2] == step.x_right]
    theta_vertical1 = [val[1] for val in coords_firstHit if val[2] == step.x_right]
    text_J_vertical1 = [f'({val[2]:.1f},{val[3]:.1f})' for val in coords_firstHit if val[2] == step.x_right]
    text_J_initial_vertical1 = [f'({val[4]:.1f})' for val in coords_firstHit if val[2] == step.x_right]

    J_horizontal1 = [val[0] for val in coords_firstHit if val[3] == step.y_up]
    theta_horizontal1 = [val[1] for val in coords_firstHit if val[3] == step.y_up]

    J_vertical2 = [val[0] for val in coords_SecondtHit if val[2] == step.x_right]
    theta_vertical2 = [val[1] for val in coords_SecondtHit if val[2] == step.x_right]
    text_J_vertical2 = [f'({val[2]:.1f},{val[3]:.1f})' for val in coords_SecondtHit if val[2] == step.x_right]
    text_J_initial_vertical2 = [f'({val[4]:.1f})' for val in coords_SecondtHit if val[2] == step.x_right]

    J_horizontal2 = [val[0] for val in coords_SecondtHit if val[3] == step.y_up]
    theta_horizontal2 = [val[1] for val in coords_SecondtHit if val[3] == step.y_up]

    #sanity check
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
    for i, txt in enumerate(x2):
        ax2.annotate(str(coords_SecondtHit[i][4]), (x2[i], y2[i]))
    ax2.set_title("second hit")
    plt.show()


    fig = plt.figure(figsize=(16, 5))  # create a figure
    ax1 = fig.add_subplot(1, 2, 1)  # create a subplot of certain size
    ax1.scatter(theta_vertical1, J_vertical1, c='r')
    # for i, txt in enumerate(text_J_vertical1):
    #     ax1.annotate(txt, (theta_vertical1[i], J_vertical1[i]))

    ax1.scatter(theta_horizontal1, J_horizontal1, c='b')
    ax1.set_title("first hit")

    ax2 = fig.add_subplot(1, 2, 2)
    ax2.scatter(theta_vertical2, J_vertical2, c='r')
    # for i, txt in enumerate(text_J_initial_vertical2):
    #     ax2.annotate(txt, (theta_vertical2[i], J_vertical2[i]))
    ax2.scatter(theta_horizontal2, J_horizontal2, c='b')
    ax2.set_title("second hit")
    plt.show()


if calc_trim:
    coords_trim_firstHit = list()
    coords_trim_SecondtHit = list()

    step = Trimmed_step()
    for i in range(-400, 401):
        stop_iteration = 0
        num_of_hits = 0
        # print(i)
        h_initial_conditions['current_q2'] = i
        h = twoD_Harmonic_separable_hamiltonian(**h_initial_conditions)

        while True:
            dynamics(hamiltonian=h, obstacle=step)

            if stop_iteration > 5000:
                break
            if step.hit_flag:
                if num_of_hits == 0:
                    J_2 = 0.5 * h.w2 * (h_initial_conditions['current_q2'] ** 2)
                    if h_initial_conditions['current_q2'] == h.previous_q2:
                        theta_wall = np.pi
                    else:
                        theta_wall = (np.pi / 2) - np.arctan(
                            h.previous_q2 / np.sqrt(h_initial_conditions['current_q2'] ** 2 - h.previous_q2 ** 2))
                        if h.previous_p2 < 0:
                            theta_wall = 2 * np.pi - theta_wall
                    coords_trim_firstHit.append((J_2, theta_wall, h.previous_q1, h.previous_q2, i))
                    num_of_hits += 1
                else:
                    E_2 = 0.5*(h.current_p2**2)+0.5*((h.w2*h.current_q2)**2)
                    J_2 = E_2/h.w2
                    q2_max = np.sqrt(2*E_2/h.w2)
                    if q2_max == h.previous_q2:
                        theta_wall = np.pi
                    else:
                        theta_wall = (np.pi / 2) - np.arctan(
                            h.previous_q2 / np.sqrt(q2_max ** 2 - h.previous_q2 ** 2))
                        if h.previous_p2 < 0:
                            theta_wall = 2 * np.pi - theta_wall
                    coords_trim_SecondtHit.append((J_2, theta_wall, h.previous_q1, h.previous_q2, i))
                    break

            stop_iteration += 1

    x1_trim = [val[2] for val in coords_trim_firstHit]
    y1_trim = [val[3] for val in coords_trim_firstHit]

    x2_trim = [val[2] for val in coords_trim_SecondtHit]
    y2_trim = [val[3] for val in coords_trim_SecondtHit]

    text_first_hit_trim = [f'({val[4]:.1f})' for val in coords_trim_firstHit]
    text_second_hit_trim = [f'({val[4]:.1f})' for val in coords_trim_SecondtHit]


    fig = plt.figure(figsize=(16, 5))  # create a figure
    ax1 = fig.add_subplot(1, 2, 1)  # create a subplot of certain size
    ax1.scatter(x1_trim, y1_trim, c='r')
    ax1.legend(['hits locations'])
    ax1.set_xlabel('x1')
    ax1.set_ylabel('y1')
    # for i, txt in enumerate(text_first_hit_trim):
    #     ax1.annotate(txt, (x1_trim[i], y1_trim[i]))

    ax2 = fig.add_subplot(1, 2, 2)  # create a subplot of certain size
    ax2.scatter(x2_trim, y2_trim, c='r')
    for i, txt in enumerate(text_second_hit_trim):
        ax2.annotate(txt, (x2_trim[i], y2_trim[i]))
    ax2.legend(['hits locations'])
    ax2.set_xlabel('x2')
    ax2.set_ylabel('y2')

    plt.show()

    J_vertical1 = [val[0] for val in coords_trim_firstHit if val[2] == step.x_right]
    theta_vertical1 = [val[1] for val in coords_trim_firstHit if val[2] == step.x_right]

    J_horizontal1 = [val[0] for val in coords_trim_firstHit if val[3] == step.y_up]
    theta_horizontal1 = [val[1] for val in coords_trim_firstHit if val[3] == step.y_up]

    J_slope1 = [val[0] for val in coords_trim_firstHit if ((val[3] != step.y_up) and (val[2] != step.x_right))]
    theta_slope1 = [val[1] for val in coords_trim_firstHit if ((val[3] != step.y_up) and (val[2] != step.x_right))]

    J_vertical2 = [val[0] for val in coords_trim_SecondtHit if val[2] == step.x_right]
    theta_vertical2 = [val[1] for val in coords_trim_SecondtHit if val[2] == step.x_right]

    J_horizontal2 = [val[0] for val in coords_trim_SecondtHit if val[3] == step.y_up]
    theta_horizontal2 = [val[1] for val in coords_trim_SecondtHit if val[3] == step.y_up]

    J_slope2 = [val[0] for val in coords_trim_SecondtHit if ((val[3] != step.y_up) and (val[2] != step.x_right))]
    theta_slope2 = [val[1] for val in coords_trim_SecondtHit if ((val[3] != step.y_up) and (val[2] != step.x_right))]

    fig = plt.figure(figsize=(16, 5))
    ax1 = fig.add_subplot(1, 2, 1)  # create a subplot of certain size
    ax1.scatter(theta_vertical1, J_vertical1, c='r')
    for i, txt in enumerate(theta_vertical1):
        ax1.annotate(str(i), (theta_vertical1[i], J_vertical1[i]))
    ax1.scatter(theta_horizontal1, J_horizontal1, c='b')
    ax1.scatter(theta_slope1, J_slope1, c='g')
    ax1.legend(['vertical side', 'horizontal side', 'slope side'])
    ax1.set_xlabel('theta1')
    ax1.set_ylabel('J1')


    ax2 = fig.add_subplot(1, 2, 2)  # create a subplot of certain size
    ax2.scatter(theta_vertical2, J_vertical2, c='r')
    for i, txt in enumerate(theta_vertical2):
        ax2.annotate(str(i), (theta_vertical2[i], J_vertical2[i]))
    ax2.scatter(theta_horizontal2, J_horizontal2, c='b')
    ax2.scatter(theta_slope2, J_slope2, c='g')
    ax2.legend(['vertical side', 'horizontal side', 'slope side'])
    ax2.set_xlabel('theta2')
    ax2.set_ylabel('J2')

    plt.show()
print("new_map")


if calc_circ:
    circ_trim_firstHit = list()
    circ_trim_SecondtHit = list()

    step = circled_edge(top_right_radius=80)

    for i in range(-400, 401):
        stop_iteration = 0
        num_of_hits = 0
        # print(i)
        h_initial_conditions['current_q2'] = i
        h = twoD_Harmonic_separable_hamiltonian(**h_initial_conditions)

        while True:
            dynamics(hamiltonian=h, obstacle=step)

            if stop_iteration > 5000:
                break
            if step.hit_flag:
                if num_of_hits == 0:
                    J_2 = 0.5 * h.w2 * (h_initial_conditions['current_q2'] ** 2)
                    if h_initial_conditions['current_q2'] == h.previous_q2:
                        theta_wall = np.pi
                    else:
                        theta_wall = (np.pi / 2) - np.arctan(
                            h.previous_q2 / np.sqrt(h_initial_conditions['current_q2'] ** 2 - h.previous_q2 ** 2))
                        if h.previous_p2 < 0:
                            theta_wall = 2 * np.pi - theta_wall
                    circ_trim_firstHit.append((J_2, theta_wall, h.previous_q1, h.previous_q2, i))
                    num_of_hits += 1
                else:
                    E_2 = 0.5*(h.current_p2**2)+0.5*((h.w2*h.current_q2)**2)
                    J_2 = E_2/h.w2
                    q2_max = np.sqrt(2*E_2/h.w2)
                    if q2_max == h.previous_q2:
                        theta_wall = np.pi
                    else:
                        theta_wall = (np.pi / 2) - np.arctan(
                            h.previous_q2 / np.sqrt(q2_max ** 2 - h.previous_q2 ** 2))
                        if h.previous_p2 < 0:
                            theta_wall = 2 * np.pi - theta_wall
                    circ_trim_SecondtHit.append((J_2, theta_wall, h.previous_q1, h.previous_q2, i))
                    break

            stop_iteration += 1

    x1_circ = [val[2] for val in circ_trim_firstHit]
    y1_circ = [val[3] for val in circ_trim_firstHit]

    x2_circ = [val[2] for val in circ_trim_SecondtHit]
    y2_circ = [val[3] for val in circ_trim_SecondtHit]

    fig = plt.figure(figsize=(16, 5))  # create a figure
    ax1 = fig.add_subplot(1, 2, 1)  # create a subplot of certain size
    ax1.scatter(x1_circ, y1_circ, c='r')
    ax1.legend(['hits locations'])
    ax1.set_xlabel('x1')
    ax1.set_ylabel('y1')
    for i, txt in enumerate(x1_circ):
        ax1.annotate(circ_trim_firstHit[i][4], (x1_circ[i], y1_circ[i]))

    ax2 = fig.add_subplot(1, 2, 2)  # create a subplot of certain size
    ax2.scatter(x2_circ, y2_circ, c='r')
    # for i, txt in enumerate(text_second_hit_trim):
    #     ax2.annotate(txt, (x2_trim[i], y2_trim[i]))
    ax2.legend(['hits locations'])
    ax2.set_xlabel('x2')
    ax2.set_ylabel('y2')

    plt.show()

    J_vertical1 = [val[0] for val in circ_trim_firstHit if val[2] == step.x_right]
    theta_vertical1 = [val[1] for val in circ_trim_firstHit if val[2] == step.x_right]

    J_horizontal1 = [val[0] for val in circ_trim_firstHit if val[3] == step.y_up]
    theta_horizontal1 = [val[1] for val in circ_trim_firstHit if val[3] == step.y_up]

    J_circ1 = [val[0] for val in circ_trim_firstHit if ((val[3] != step.y_up) and (val[2] != step.x_right))]
    theta_circ1 = [val[1] for val in circ_trim_firstHit if ((val[3] != step.y_up) and (val[2] != step.x_right))]

    J_vertical2 = [val[0] for val in circ_trim_SecondtHit if val[2] == step.x_right]
    theta_vertical2 = [val[1] for val in circ_trim_SecondtHit if val[2] == step.x_right]

    J_horizontal2 = [val[0] for val in circ_trim_SecondtHit if val[3] == step.y_up]
    theta_horizontal2 = [val[1] for val in circ_trim_SecondtHit if val[3] == step.y_up]

    J_circ2 = [val[0] for val in circ_trim_SecondtHit if ((val[3] != step.y_up) and (val[2] != step.x_right))]
    theta_circ2 = [val[1] for val in circ_trim_SecondtHit if ((val[3] != step.y_up) and (val[2] != step.x_right))]

    fig = plt.figure(figsize=(16, 5))
    ax1 = fig.add_subplot(1, 2, 1)  # create a subplot of certain size
    ax1.scatter(theta_vertical1, J_vertical1, c='r')
    # for i, txt in enumerate(theta_vertical1):
    #     ax1.annotate(str(i), (theta_vertical1[i], J_vertical1[i]))
    ax1.scatter(theta_horizontal1, J_horizontal1, c='b')
    ax1.scatter(theta_circ1, J_circ1, c='g')
    ax1.legend(['vertical side', 'horizontal side', 'slope side'])
    ax1.set_xlabel('theta1')
    ax1.set_ylabel('J1')


    ax2 = fig.add_subplot(1, 2, 2)  # create a subplot of certain size
    ax2.scatter(theta_vertical2, J_vertical2, c='r')
    # for i, txt in enumerate(theta_vertical2):
    #     ax2.annotate(str(i), (theta_vertical2[i], J_vertical2[i]))
    ax2.scatter(theta_horizontal2, J_horizontal2, c='b')
    ax2.scatter(theta_circ2, J_circ2, c='g')
    ax2.legend(['vertical side', 'horizontal side', 'slope side'])
    ax2.set_xlabel('theta2')
    ax2.set_ylabel('J2')

    plt.show()