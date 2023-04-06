import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

show_all_initial_points = 0
path_allPoints = os.path.join('angleAction', 'saved_results', '16_01_2023', 'all_board_r40.pickle')
obj_allPoints = pd.read_pickle(path_allPoints)

# HUSH MAPS
ids_first_hit = [val[5] for val in obj_allPoints['first_hit']]
indexed_first_hit = [i for i in range(len(ids_first_hit))]
hushMap_first_hit = dict(zip(ids_first_hit, indexed_first_hit))

ids_second_hit = [val[5] for val in obj_allPoints['second_hit']]
indexed_second_hit = [i for i in range(len(ids_second_hit))]
hushMap_second_hit = dict(zip(ids_second_hit, indexed_second_hit))

ids_round_fh = [val[5] for val in obj_allPoints['first_hit'] if (val[2] != -50) and (val[3] != -50)]
indexed_round_fh = [i for i in range(len(ids_round_fh))]
hushMap_round_fh = dict(zip(ids_round_fh, indexed_round_fh))

ids_round_sh = [val[5] for val in obj_allPoints['second_hit'] if (val[2] != -50) and (val[3] != -50)]
indexed_round_sh = [i for i in range(len(ids_round_sh))]
hushMap_round_sh = dict(zip(ids_round_sh, indexed_round_sh))

ids_sh_max_vals = [(np.round(val[4][0]).astype(int), np.rint(val[4][1]).astype(int)) for val in
                   obj_allPoints['second_hit']]
hushMap_sh_max_vals = dict(zip(ids_sh_max_vals, indexed_second_hit))

ids_fh_max_vals = [(val[4][0], val[4][1]) for val in obj_allPoints['first_hit']]
hushMap_fh_max_vals = dict(zip(ids_fh_max_vals, indexed_first_hit))

# Box 3

x_fh_round_box3_fh = [obj_allPoints['first_hit'][hushMap_round_fh[idx]][4][0] for idx in ids_round_fh]
y_fh_round_box3_fh = [obj_allPoints['first_hit'][hushMap_round_fh[idx]][4][1] for idx in ids_round_fh]

x_fh_round_box3_sh = [obj_allPoints['second_hit'][hushMap_round_fh[idx]][4][0] for idx in ids_round_fh]
y_fh_round_box3_sh = [obj_allPoints['second_hit'][hushMap_round_fh[idx]][4][1] for idx in ids_round_fh]

x_sh_round_box3_fh = [obj_allPoints['first_hit'][hushMap_round_sh[idx]][4][0] for idx in ids_round_sh]
y_sh_round_box3_fh = [obj_allPoints['first_hit'][hushMap_round_sh[idx]][4][1] for idx in ids_round_sh]

x_sh_round_box3_sh = [obj_allPoints['second_hit'][hushMap_round_sh[idx]][4][0] for idx in ids_round_sh]
y_sh_round_box3_sh = [obj_allPoints['second_hit'][hushMap_round_sh[idx]][4][1] for idx in ids_round_sh]

fig = plt.figure(figsize=(16, 5))  # create a figure
ax1 = fig.add_subplot(1, 2, 1)  # create a subplot of certain size
ax1.scatter(x_fh_round_box3_fh, y_fh_round_box3_fh, c='r', marker='^')
ax1.scatter(x_fh_round_box3_sh, y_fh_round_box3_sh, c='b', marker='v')
ax1.set_xlim([-400, 400])
ax1.set_ylim([-400, 400])
ax1.set_title("initial_right")
ax1.legend(['FH - FH', 'FH - SH'])

ax2 = fig.add_subplot(1, 2, 2)  # create a subplot of certain size
ax2.scatter(x_sh_round_box3_fh, y_sh_round_box3_fh, c='r', marker='^')
ax2.scatter(x_sh_round_box3_sh, y_sh_round_box3_sh, c='b', marker='v')
ax2.set_xlim([-400, 400])
ax2.set_ylim([-400, 400])
ax2.set_title("initial_right")
ax2.legend(['SH - FH', 'SH - SH'])
plt.show()

# Box 2
box21_ids = [indx[5] for indx in obj_allPoints['first_hit'] if
             (indx[4][0] >= 250 and indx[4][0] <= 300) and (
                     indx[4][1] >= -100 and indx[4][1] <= -50)]

box22_ids = [indx[5] for indx in obj_allPoints['first_hit'] if
             (indx[4][0] >= 210 and indx[4][0] <= 230) and (
                     indx[4][1] >= 110 and indx[4][1] <= 130)]

x_initial_box21_fh = [obj_allPoints['first_hit'][hushMap_first_hit[idx]][4][0] for idx in box21_ids]
y_initial_box21_fh = [obj_allPoints['first_hit'][hushMap_first_hit[idx]][4][1] for idx in box21_ids]

x_initial_box21_sh = [obj_allPoints['second_hit'][hushMap_second_hit[idx]][4][0] for idx in box21_ids]
y_initial_box21_sh = [obj_allPoints['second_hit'][hushMap_second_hit[idx]][4][1] for idx in box21_ids]

x_initial_box22_fh = [obj_allPoints['first_hit'][hushMap_first_hit[idx]][4][0] for idx in box22_ids]
y_initial_box22_fh = [obj_allPoints['first_hit'][hushMap_first_hit[idx]][4][1] for idx in box22_ids]

x_initial_box22_sh = [obj_allPoints['second_hit'][hushMap_second_hit[idx]][4][0] for idx in box22_ids]
y_initial_box22_sh = [obj_allPoints['second_hit'][hushMap_second_hit[idx]][4][1] for idx in box22_ids]

fig = plt.figure(figsize=(16, 5))  # create a figure
ax1 = fig.add_subplot(1, 2, 1)  # create a subplot of certain size
ax1.scatter(x_initial_box21_fh, y_initial_box21_fh, c='r', marker='^')
ax1.scatter(x_initial_box21_sh, y_initial_box21_sh, c='b', marker='v')
ax1.set_xlim([-400, 400])
ax1.set_ylim([-400, 400])
ax1.set_title("initial_right")
ax1.legend(['initial points', 'max_after_hit'])

ax2 = fig.add_subplot(1, 2, 2)  # create a subplot of certain size
ax2.scatter(x_initial_box22_fh, y_initial_box22_fh, c='r', marker='^')
ax2.scatter(x_initial_box22_sh, y_initial_box22_sh, c='b', marker='v')
ax2.set_xlim([-400, 400])
ax2.set_ylim([-400, 400])
ax2.set_title("initial_right")
ax2.legend(['initial points', 'max_after_hit'])
plt.show()

J2_box2_fs = [obj_allPoints['first_hit'][hushMap_first_hit[idx]][0] for idx in box2_ids]
theta_box2_fs = [obj_allPoints['first_hit'][hushMap_first_hit[idx]][1] for idx in box2_ids]

J2_box2_sh = [obj_allPoints['second_hit'][hushMap_second_hit[idx]][0] for idx in box2_ids]
theta_box2_sh = [obj_allPoints['second_hit'][hushMap_second_hit[idx]][1] for idx in box2_ids]

fig = plt.figure(figsize=(8, 5))  # create a figure
ax1 = fig.add_subplot(1, 1, 1)  # create a subplot of certain size
ax1.scatter(theta_box2_fs, J2_box2_fs, c='r', marker='^')
ax1.scatter(theta_box2_sh, J2_box2_sh, c='b', marker='v')
# ax1.set_xlim([-400, 400])
# ax1.set_ylim([-400, 400])
ax1.set_title("initial_right")
ax1.legend(['initial points', 'max_after_hit'])
plt.show()

# Box 1
box1_ids = [indx[5] for indx in obj_allPoints['first_hit'] if
            (indx[4][0] >= 350 and indx[4][0] <= 400) and (
                    indx[4][1] >= -400 and indx[4][1] <= -350)]

x_initial_box1_fh = [obj_allPoints['first_hit'][hushMap_first_hit[idx]][4][0] for idx in box1_ids]
y_initial_box1_fh = [obj_allPoints['first_hit'][hushMap_first_hit[idx]][4][1] for idx in box1_ids]

x_initial_box1_sh = [obj_allPoints['second_hit'][hushMap_second_hit[idx]][4][0] for idx in box1_ids]
y_initial_box1_sh = [obj_allPoints['second_hit'][hushMap_second_hit[idx]][4][1] for idx in box1_ids]

fig = plt.figure(figsize=(8, 5))  # create a figure
ax1 = fig.add_subplot(1, 1, 1)  # create a subplot of certain size
ax1.scatter(x_initial_box1_fh, y_initial_box1_fh, c='r', marker='^')
ax1.scatter(x_initial_box1_sh, y_initial_box1_sh, c='b', marker='v')
ax1.set_xlim([-400, 400])
ax1.set_ylim([-400, 400])
ax1.set_title("initial_right")
ax1.legend(['initial points', 'max_after_hit'])
plt.show()

# show the initial conditions:


if show_all_initial_points:
    x_initial = [val[4][0] for val in obj_allPoints['first_hit']]
    y_initial = [val[4][1] for val in obj_allPoints['first_hit']]
    fig = plt.figure(figsize=(16, 5))  # create a figure
    ax1 = fig.add_subplot(1, 1, 1)  # create a subplot of certain size
    ax1.scatter(x_initial, y_initial, c='r')
    # ax1.scatter(x_initial_left, y_initial_left, c='b')
    ax1.set_xlim([-400, 400])
    ax1.set_ylim([-400, 400])
    ax1.set_title("initial_right")
    ax1.legend(['initial points'])
    plt.show()

#
# ax2 = fig.add_subplot(1, 2, 2)
# ax2.scatter(x_second_hit_max_right, y_second_hit_max_right, c='r')
# ax2.scatter(x_second_hit_max_left, y_second_hit_max_left, c='b')
# ax2.set_xlim([-400, 400])
# ax2.set_ylim([-400, 400])
# ax2.set_title("second hit max values")
# ax2.legend(['initial right', 'initial left'])


print('check')

# filter indexes
