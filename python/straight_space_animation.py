import pygame
import numpy as np

# to reduce calculations:
global sqrt2
sqrt2 = np.sqrt(2)

colors = {'BLACK': (0, 0, 0),
          'BLACK_lighter': (20, 20, 20),
          'GRAY': (127, 127, 127),
          'WHITE': (255, 255, 255),
          'YELLOW': (255, 255, 0),
          'CYAN': (0, 255, 255),
          'MAGENTA': (255, 0, 255),
          'RED': (255, 0, 0),
          'GREEN': (0, 255, 0)
          }

h_initial_conditions = {
    'current_q1': 375,
    'current_q2': -375,
    'w1': 1,
    'w2': 2,
    'res': (800, 800)
}


class twoD_Harmonic_separable_hamiltonian:
    def __init__(self, **kwargs):
        self.m1 = kwargs.setdefault('m1', 1)
        self.m2 = kwargs.setdefault('m2', 1)
        self.w1 = kwargs.setdefault('w1', 1)
        self.w2 = kwargs.setdefault('w2', 1)
        self.current_q1 = kwargs.setdefault('current_q1', 0)
        self.current_p1 = kwargs.setdefault('current_p1', 0)
        self.current_q2 = kwargs.setdefault('current_q2', 0)
        self.current_p2 = kwargs.setdefault('current_p2', 0)
        self.previous_q1 = kwargs.setdefault('previous_q1', self.current_q1)
        self.previous_p1 = kwargs.setdefault('previous_p1', self.current_p1)
        self.previous_q2 = kwargs.setdefault('previous_q2', self.current_q2)
        self.previous_p2 = kwargs.setdefault('previous_p2', self.current_p2)
        self.res = kwargs.setdefault('res', (0, 0))
        self.radius = kwargs.setdefault('radius', 5)
        self.center = kwargs.setdefault('center', (0, 0))
        self.dt = kwargs.setdefault('dt', 0.01)
        # print((self.current_q1, self.current_p1, self.current_q2, self.current_p2))
        self.energy_1, self.energy_2, self.q1_max, self.q2_max = self.calc_max_q1_q2()

    def calc_max_q1_q2(self):
        # E_i = p^2/2m + 1/2*w^2*q^2
        energy_1 = (self.current_p1 ** 2) / (2 * self.m1) + 0.5 * (self.w1 ** 2) * (self.current_q1 ** 2)
        energy_2 = (self.current_p2 ** 2) / (2 * self.m2) + 0.5 * (self.w2 ** 2) * (self.current_q2 ** 2)
        # print("the total energy is" + str(energy_1+energy_2))

        q1_max = (sqrt2 * np.sqrt(energy_1)) / self.w1
        q2_max = (sqrt2 * np.sqrt(energy_2)) / self.w2
        # print(f'q1 max = {q1_max}, and q2_max = {q2_max}')
        return energy_1, energy_2, q1_max, q2_max

    def step(self):
        self.previous_q1 = self.current_q1
        self.previous_p1 = self.current_p1
        self.previous_q2 = self.current_q2
        self.previous_p2 = self.current_p2

        self.current_q1 = self.previous_q1 + self.dt * (self.previous_p1 // self.m1)
        self.current_p1 = self.previous_p1 - self.dt * (self.w1 ** 2) * self.current_q1
        self.current_q2 = self.previous_q2 + self.dt * (self.previous_p2 // self.m2)
        self.current_p2 = self.previous_p2 - self.dt * (self.w2 ** 2) * self.current_q2

        # print((self.current_q1, self.current_p1, self.current_q2, self.current_p2))
        # TODO: maybe I have a problem with the iterations here
        self.center = (self.current_q1 + self.res[0] // 2, -self.current_q2 + self.res[1] // 2)
        self.energy_1, self.energy_2, self.q1_max, self.q2_max = self.calc_max_q1_q2()


def rotate(x, y, theta):
    new_x = x * np.cos(theta) - y * np.sin(theta)
    new_y = x * np.sin(theta) + y * np.cos(theta)
    return new_x, new_y


class Step:
    def __init__(self, **kwargs):
        """
        The cordinats are NOT in the pyGame board
        :param kwargs:
        """
        self.x_right = kwargs.setdefault('x_right', -50)
        self.y_up = kwargs.setdefault('y_up', -50)
        self.x_left = kwargs.setdefault('x_left', -400)
        self.y_down = kwargs.setdefault('y_down', -400)
        self.hit_flag = 0
        self.epsilon = kwargs.setdefault('epsilon', 0.1)

    def convert_to_rect_params(self, width, height):
        left = self.x_left + width // 2
        top = -self.y_up + height // 2
        width = self.x_right - self.x_left
        height = self.y_up - self.y_down

        return left, top, width, height

    def inside(self, q1, q2):
        return (q1 <= self.x_right) and (q2 <= self.y_up)

    def return_rule(self, q1, q2, p1, p2):

        # if (q1 > self.x_right) and (q2 > self.y_up):
        #     print(f'out - ({q1}, {q2}))')
        #     self.flag = 1

        if self.inside(q1=q1, q2=q2):
            # print(f'inside - ({q1}, {q2}))')
            if not self.hit_flag:
                self.hit_flag = 1
                q1_tilde = q1 - p1 * self.epsilon
                q2_tilde = q2 + p2 * self.epsilon
                condition = self.inside(q1=q1_tilde, q2=q2_tilde)
                # if np.abs(q1 - self.x_right) < np.abs(q2 - self.y_up):
                if not condition and np.abs(q1 - self.x_right) < np.abs(q2 - self.y_up):
                    dx = self.x_right - q1
                    q1 = self.x_right
                    q2 = q2 + dx * (p2 / p1)
                    return q1, q2, -p1, p2
                else:
                    dy = self.y_up - q2
                    q2 = self.y_up
                    q1 = q1 + dy * (p1 / p2)
                    return q1, q2, p1, -p2

        else:
            self.hit_flag = 0
        # elif (q1 <= self.x_right) and (q1 >= (self.x_right - self.epsilon)):
        #     print("over the left")
        #     if q2 < self.y_up and self.flag:
        #         self.flag = 0
        #         return -p1, p2
        #
        # elif (q2 <= self.y_up) and (q2 > (self.y_up - self.epsilon)):
        #     print("over the upper")
        #     if q1 < self.x_right and self.flag:
        #         self.flag = 0
        #         return p1, -p2

        return q1, q2, p1, p2


class circled_edge(Step):
    def __init__(self, **kwargs):
        """
        The cordinats are NOT in the pyGame board
        :param kwargs:
        """
        Step.__init__(self, **kwargs)
        self.top_right_radius = kwargs.setdefault('top_right_radius', 30)
        self.x_center = self.x_right - self.top_right_radius
        self.y_center = self.y_up - self.top_right_radius

    def inside(self, q1, q2):
        if (q1 <= self.x_right) and (q2 <= self.y_up):
            if q1 < self.x_center:
                return 1
            elif q2 < self.y_center:
                return 2
            else:
                dist = np.sqrt((q1 - self.x_center) ** 2 + (q2 - self.y_center) ** 2)
                if dist < self.top_right_radius:
                    return 3

                return 0
        return 0

    def calculated_returned_momentum(self, p1, p2, slope):
        """
        calculate the return values of the momentum, by the formula of:
        |p1_return| = |cos(2*theta)      sin(2*theta)    | |p1|
        |p2_return| = |sin(2 * theta)    -cos(2 * theta) | |p2|
        :param p1:
        :param p2:
        :return: p1_return, p2_return
        """
        theta = np.arctan(slope)
        p1_return = p1 * np.cos(2 * theta) + p2 * np.sin(2 * theta)
        p2_return = p1 * np.sin(2 * theta) - p2 * np.cos(2 * theta)
        return p1_return, p2_return

    def reconstruct_x_y(self, q1, q2, p1, p2):

        m = p2 / p1
        b = q2 - m * q1
        h = self.x_center
        k = self.y_center
        r = self.top_right_radius
        A = 1 + m ** 2
        B = 2 * m * (b - k) - 2 * h
        C = (b - k) ** 2 + h ** 2 - r ** 2
        rec_x_plus = (-B + np.sqrt(B ** 2 - 4 * A * C)) / (2 * A)
        rec_x_minus = (-B - np.sqrt(B ** 2 - 4 * A * C)) / (2 * A)
        rec_x = rec_x_plus if (np.abs(rec_x_plus - q1) < np.abs(rec_x_minus - q1)) else rec_x_minus
        rec_y = m * rec_x + b

        return rec_x, rec_y

    def return_rule(self, q1, q2, p1, p2):
        switch = self.inside(q1, q2)
        if switch:
            if not self.hit_flag:
                self.hit_flag = 1
                if switch == 1:
                    dy = self.y_up - q2
                    q2 = self.y_up
                    q1 = q1 + dy * (p1 / p2)
                    return q1, q2, p1, -p2
                if switch == 2:
                    dx = self.x_right - q1
                    q1 = self.x_right
                    q2 = q2 + dx * (p2 / p1)
                    return q1, q2, -p1, p2
                if switch == 3:
                    rec_x, rec_y = self.reconstruct_x_y(q1, q2, p1, p2)
                    current_slope = (self.x_center - rec_x) / (rec_y - self.y_center)
                    # rotated_p2, rotated_p1 = rotate(p1, p2, np.arctan(current_slope))
                    # b1 = q2 - current_slope * q1
                    # b2 = q2 + q1 / current_slope
                    # perpen_q1 = (current_slope / (1 + current_slope ** 2)) * (b2 - b1)
                    # perpen_q2 = perpen_q1 * current_slope + b1
                    # dist = np.sqrt((q1 - perpen_q1) ** 2 + (q2 - perpen_q2) ** 2)
                    # rotated_fix = dist * (rotated_p1 / rotated_p2)
                    # if np.abs(rotated_fix) > 5:
                    #     rotated_fix = dist * (rotated_p2 / rotated_p1)
                    # rect_dx, rect_dy = rotate(0, rotated_fix, -np.arctan(current_slope))
                    p1_return, p2_return = self.calculated_returned_momentum(p1, p2, slope=current_slope)
                    return (rec_x), (rec_y), p1_return, p2_return
        else:
            self.hit_flag = 0

        return q1, q2, p1, p2


class Trimmed_step:
    def __init__(self, **kwargs):
        """
        The cordinats are NOT in the pyGame board
        :param kwargs:
        """
        self.x_right = kwargs.setdefault('x_right', -100)
        self.y_up = kwargs.setdefault('y_up', -100)
        self.x_left = kwargs.setdefault('x_left', -500)
        self.y_down = kwargs.setdefault('y_down', -500)
        self.x_trim = kwargs.setdefault('x_trim', -200)
        self.y_trim = kwargs.setdefault('y_trim', -200)
        self.epsilon = kwargs.setdefault('epsilon', 5)
        self.hit_flag = 1
        self.dt = kwargs.setdefault('dt', 0.01)

        if self.x_right != self.x_trim:
            self.slope = (self.y_trim - self.y_up) / (self.x_right - self.x_trim)
        else:
            self.slope = np.inf

    def convert_to_polygon_params(self, width, height):
        top_left = (self.x_left + width // 2, -self.y_up + height // 2)
        trim_up_left = (self.x_trim + width // 2, -self.y_up + height // 2)
        trim_down_right = (self.x_right + width // 2, -self.y_trim + height // 2)
        bottom_right = (self.x_right + width // 2, -self.y_down + height // 2)
        corner = (self.x_left + width // 2, -self.y_down + height // 2)

        return [top_left, trim_up_left, trim_down_right, bottom_right, corner]

    def calculated_returned_momentum(self, p1, p2):
        """
        calculate the return values of the momentum, by the formula of:
        |p1_return| = |cos(2*theta)      sin(2*theta)    | |p1|
        |p2_return| = |sin(2 * theta)    -cos(2 * theta) | |p2|
        :param p1:
        :param p2:
        :return: p1_return, p2_return
        """
        theta = np.arctan(self.slope)
        p1_return = p1 * np.cos(2 * theta) + p2 * np.sin(2 * theta)
        p2_return = p1 * np.sin(2 * theta) - p2 * np.cos(2 * theta)
        return p1_return, p2_return

    def inside(self, q1, q2):
        if q1 <= self.x_right and q2 <= self.y_up:
            if q1 <= self.x_trim:
                return 1
            elif q2 <= self.y_trim:
                return 2
            elif (q2 - self.y_up) < ((q1 - self.x_trim) * self.slope):
                # elif q2 < (self.slope*(q1-self.x_trim) + self.y_up):
                # print("this is inside the return 3")
                return 3
            else:
                # print(f'this is the slope y: {((q1 - self.x_trim) * self.slope) + self.y_up}')
                # print(f'the cordiantes are: ({q1}, {q2})')
                return 0
        return 0

    def return_rule(self, q1, q2, p1, p2):

        switch = self.inside(q1, q2)
        if switch:
            if not self.hit_flag:
                self.hit_flag = 1
                if switch == 1:
                    dy = self.y_up - q2
                    q2 = self.y_up
                    q1 = q1 + dy * (p1 / p2)
                    return q1, q2, p1, -p2
                if switch == 2:
                    dx = self.x_right - q1
                    q1 = self.x_right
                    q2 = q2 + dx * (p2 / p1)
                    return q1, q2, -p1, p2
                if switch == 3:
                    rotated_p2, rotated_p1 = rotate(p1, p2, np.arctan(self.slope))
                    b1 = self.y_up - self.slope * self.x_trim
                    b2 = q2 + q1 / self.slope
                    perpen_q1 = (self.slope / (1 + self.slope ** 2)) * (b2 - b1)
                    perpen_q2 = perpen_q1 * self.slope + b1
                    dist = np.sqrt((q1 - perpen_q1) ** 2 + (q2 - perpen_q2) ** 2)
                    rotated_fix = dist * (rotated_p1 / rotated_p2)
                    if np.abs(rotated_fix) > 5:
                        rotated_fix = dist * (rotated_p2 / rotated_p1)
                    rect_dx, rect_dy = rotate(0, rotated_fix, -np.arctan(self.slope))
                    p1_return, p2_return = self.calculated_returned_momentum(p1, p2)
                    return (perpen_q1 + rect_dy), (perpen_q2 + rect_dx), p1_return, p2_return
        else:
            self.hit_flag = 0
        # if (q1 < self.x_right) and (q1 > (self.x_right - self.epsilon)):
        #     if q2 < self.y_trim and not self.hit_flag:
        #         self.hit_flag = 0
        #         return (q1 - self.dt*p1), (q2 - self.dt*p2), -p1, p2
        #
        # if (q2 < self.y_up) and (q2 > (self.y_up - self.epsilon)):
        #     if q1 < self.x_trim and not self.hit_flag:
        #         self.hit_flag = 0
        #         return (q1 - self.dt*p1), (q2 - self.dt*p2), p1, -p2
        #
        # if (q1 < self.x_right) and (q2 < self.y_up):
        #     line_val = (q1 - self.x_trim) * self.slope
        #     if ((q2 - self.y_up) < line_val) and not self.hit_flag:
        #         self.hit_flag = 0
        #         # print("now flag is" + str(self.flag))
        #         p1_return, p2_return = self.calculated_returned_momentum(p1, p2)
        #         return (q1 - self.dt*p1), (q2 - self.dt*p2), p1_return, p2_return

        # if (q1 > self.x_right) and (q2 > self.y_up):
        #     self.hit_flag = 0

        return q1, q2, p1, p2


def dynamics(hamiltonian, screen=None, obstacle=None):
    if obstacle:
        hamiltonian.current_q1, hamiltonian.current_q2, hamiltonian.current_p1, hamiltonian.current_p2 = obstacle.return_rule(
            hamiltonian.current_q1,
            hamiltonian.current_q2,
            hamiltonian.current_p1,
            hamiltonian.current_p2)

    hamiltonian.step()
    if screen:
        pygame.draw.circle(screen, colors['WHITE'],
                           center=hamiltonian.center,
                           radius=hamiltonian.radius)


def main():
    # initializing the constructor
    pygame.init()
    pygame.font.init()
    # screen resolution
    res = h_initial_conditions['res']

    # opens up a window
    screen = pygame.display.set_mode(res)

    width, height = screen.get_width(), screen.get_height()

    # Screen Update Speed (FPS)
    clock = pygame.time.Clock()

    h = twoD_Harmonic_separable_hamiltonian(**h_initial_conditions)
    step1 = Step()
    trimmed_step = Trimmed_step()
    circled = circled_edge(top_right_radius=40)

    while True:
        # fills the screen with a color
        screen.fill(colors['BLACK'])

        # show axis
        pygame.draw.line(screen, colors['YELLOW'], (0, height // 2), (width, height // 2))  # X-axis
        pygame.draw.line(screen, colors['YELLOW'], (width // 2, 0), (width // 2, height))  # Y-axis
        pygame.draw.line(screen, colors['YELLOW'], (0, 0), (width, height))  # X=-Y-axis
        pygame.draw.line(screen, colors['YELLOW'], (width, 0), (0, height))  # X=Y-axis
        # NO OBSTACLE:
        # dynamics(screen=screen, hamiltonian=h)

        # RECT OBSTACLE:
        # pygame.draw.rect(screen, colors['RED'], step1.convert_to_rect_params(width, height))
        # dynamics(screen=screen, hamiltonian=h, obstacle=step1)

        # RECT-TRIMMED OBSTACLE:
        # pygame.draw.polygon(screen, colors['CYAN'], trimmed_step.convert_to_polygon_params(width, height))
        # dynamics(hamiltonian=h, screen=screen, obstacle=trimmed_step)

        # CIRC OBSTACLE:
        pygame.draw.rect(screen, colors['RED'], circled.convert_to_rect_params(width, height),
                         border_top_right_radius=circled.top_right_radius)
        dynamics(screen=screen, hamiltonian=h, obstacle=circled)

        pygame.display.update()

        # Setting FPS
        clock.tick(50)


if __name__ == '__main__':
    main()
