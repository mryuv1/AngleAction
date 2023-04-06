import pygame
import numpy as np
from straight_space_animation import colors, twoD_Harmonic_separable_hamiltonian, dynamics, Step


class AngleSpace:
    def __init__(self, hamiltonian, obstacle=None):
        self.hamiltonian = hamiltonian
        self.straight_obstacle = obstacle
        self.current_theta1 = 0
        self.current_theta2 = 0
        self.res = self.hamiltonian.res
        self.epsillon = 2*np.pi/90
        if obstacle:
            self.theta_wall_right, self.theta_wall_left, self.theta_wall_up, self.theta_wall_down \
                = self.calc_obstacle_params()

    def calc_obstacle_params(self):
        arg_for_q1_right = self.straight_obstacle.x_right / (
            np.sqrt(np.abs((self.hamiltonian.q1_max ** 2) - (self.straight_obstacle.x_right ** 2))))
        theta_wall_right = np.pi / 2 - np.arctan(arg_for_q1_right)
        # arg_for_q1_left = self.straight_obstacle.x_left / (
        #     np.sqrt((self.hamiltonian.q1_max ** 2) - (self.straight_obstacle.x_left ** 2)))
        theta_wall_left = 2*np.pi - theta_wall_right

        arg_for_q2_up = self.straight_obstacle.y_up / (
            np.sqrt(np.abs((self.hamiltonian.q2_max ** 2) - (self.straight_obstacle.y_up ** 2))))
        theta_wall_up = np.pi / 2 - np.arctan(arg_for_q2_up)
        theta_wall_down = 2*np.pi - theta_wall_up

        return theta_wall_right, theta_wall_left, theta_wall_up, theta_wall_down

    def step(self):
        # for axis=1
        self.current_theta1 = (self.current_theta1 + self.hamiltonian.w1 * self.hamiltonian.dt) % (2 * np.pi)
        # for axis=2
        self.current_theta2 = (self.current_theta2 + self.hamiltonian.w2 * self.hamiltonian.dt) % (2 * np.pi)

        if self.straight_obstacle:
            # if (2 * np.pi - self.current_theta2) < self.epsillon:
            #     print(f'{self.current_theta1}, {self.current_theta2}')
            if (self.current_theta1 < self.theta_wall_left and self.current_theta1 > self.theta_wall_right) and (
                    self.current_theta2 < self.theta_wall_down and self.current_theta1 > self.theta_wall_up):

                if np.abs(self.current_theta1 - self.theta_wall_right) < self.epsillon:
                    self.current_theta1 = 2 * np.pi - self.current_theta1

                elif np.abs(self.current_theta2 - self.theta_wall_up) < self.epsillon:
                    self.current_theta2 = 2 * np.pi - self.current_theta2


        # print(f'current_theta1 = {self.current_theta1/np.pi} ')
        # print(f'current_theta2 = {self.current_theta2 / np.pi} ')

    def convert_from_rad_to_pixels(self, rad):
        return np.round(((rad % (2 * np.pi)) / (2 * np.pi)) * self.res[0])

    def plot_on_screen(self, screen):
        self.step()
        convert_theta1_to_pixels = self.convert_from_rad_to_pixels(self.current_theta1)
        convert_theta2_to_pixels = self.convert_from_rad_to_pixels(self.current_theta2)

        if self.straight_obstacle:
            convert_thetaWallRight_to_pixels = self.convert_from_rad_to_pixels(self.theta_wall_right)
            convert_thetaWallLeft_to_pixels = self.convert_from_rad_to_pixels(self.theta_wall_left)
            convert_thetaWallUp_to_pixels = self.convert_from_rad_to_pixels(self.theta_wall_up)
            convert_thetaWallDown_to_pixels = self.convert_from_rad_to_pixels(self.theta_wall_down)

            point1 = (convert_thetaWallLeft_to_pixels, convert_thetaWallUp_to_pixels)
            point2 = (convert_thetaWallRight_to_pixels, convert_thetaWallUp_to_pixels)
            point3 = (convert_thetaWallRight_to_pixels, convert_thetaWallDown_to_pixels)
            point4 = (convert_thetaWallLeft_to_pixels, convert_thetaWallDown_to_pixels)

            pygame.draw.polygon(screen, colors['CYAN'], [point1, point2, point3, point4])

        # print(f'convert_theta1_to_pixels = {convert_theta1_to_pixels} ')
        # print(f'convert_theta2_to_pixels = {convert_theta2_to_pixels} ')

        pygame.draw.circle(screen, colors['RED'],
                           center=(
                               convert_theta1_to_pixels, -convert_theta2_to_pixels + self.res[1]),
                           radius=self.hamiltonian.radius)


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

    # No obstcle:
    # h = twoD_Harmonic_separable_hamiltonian(**h_initial_conditions)
    # h_angle = AngleSpace(hamiltonian=h)

    h = twoD_Harmonic_separable_hamiltonian(**h_initial_conditions)
    step1 = Step()
    h_angle = AngleSpace(hamiltonian=h, obstacle=step1)

    while True:
        # fills the screen with a color
        screen.fill(colors['BLACK_lighter'])

        # show axis
        pygame.draw.line(screen, colors['MAGENTA'], (0, height // 2), (width, height // 2))  # X-axis
        pygame.draw.line(screen, colors['MAGENTA'], (width // 2, 0), (width // 2, height))  # Y-axis
        pygame.draw.line(screen, colors['MAGENTA'], (0, 0), (width, height))  # X=-Y-axis
        pygame.draw.line(screen, colors['MAGENTA'], (width, 0), (0, height))  # X=Y-axis
        # NO OBSTACLE:
        # dynamics(screen=screen, hamiltonian=h)

        # RECT OBSTACLE:
        # pygame.draw.rect(screen, colors['RED'], step1.convert_to_rect_params(width, height))
        # dynamics(screen=screen, hamiltonian=h, obstacle=step1)

        # RECT-TRIMMED OBSTACLE:

        dynamics(hamiltonian=h)
        h_angle.plot_on_screen(screen=screen)
        pygame.display.update()

        # Setting FPS
        clock.tick(50)


if __name__ == '__main__':
    main()
