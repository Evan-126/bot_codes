import math as m
import matplotlib.pyplot as plt
import numpy as np

def linear_position(L2, L3, t2):
    asin_arg = (L2 * m.sin(m.radians(t2))) / L3
    if asin_arg < -1 or asin_arg > 1:
        return float('nan')
    t3 = 180 - (t2 + m.degrees(m.asin(asin_arg)))
    c = m.sqrt(L2**2 + L3**2 - 2*L2*L3*m.cos(m.radians(t3)))
    return c


def linear_position_reverse(L2, L3, t2):
    asin_arg = (L2 * m.sin(m.radians(abs(t2)))) / L3
    if asin_arg < -1 or asin_arg > 1:
        return float('nan')
    t3 = 180 - (abs(t2) + m.degrees(m.asin(asin_arg)))
    c = m.sqrt(L2**2 + L3**2 - 2*L2*L3*m.cos(m.radians(t3)))
    return c


def slider_velocity(L2, L3, t2, omega, delta_t=1e-4):
    pos1 = linear_position(L2, L3, t2) if t2 >= 0 else linear_position_reverse(L2, L3, t2)
    t2_next = t2 + omega * delta_t * 180 / m.pi  # degrees increment
    pos2 = linear_position(L2, L3, t2_next) if t2_next >= 0 else linear_position_reverse(L2, L3, t2_next)
    if np.isnan(pos1) or np.isnan(pos2):
        return float('nan')
    return (pos2 - pos1) / delta_t


def rod_slider_angle(L2, L3, t2):
    x_crank = L2 * m.cos(m.radians(t2))
    y_crank = L2 * m.sin(m.radians(t2))
    if t2 >= 0:
        x_slider = linear_position(L2, L3, t2)
    else:
        x_slider = linear_position_reverse(L2, L3, t2)
    y_slider = 0
    if np.isnan(x_slider):
        return float('nan')
    angle = m.degrees(m.atan2(y_slider - y_crank, x_slider - x_crank))
    return angle


def rod_angular_velocity(L2, L3, t2, omega, delta_t=1e-4):
    angle1 = rod_slider_angle(L2, L3, t2)
    angle2 = rod_slider_angle(L2, L3, t2 + omega * delta_t * 180 / m.pi)
    if np.isnan(angle1) or np.isnan(angle2):
        return float('nan')
    ang_vel_deg = (angle2 - angle1) / delta_t
    return abs(m.radians(ang_vel_deg))


def slider_acceleration(L2, L3, t2, omega, delta_t=1e-4):
    # Numerical differentiation of slider velocity to estimate acceleration
    vel1 = slider_velocity(L2, L3, t2, omega, delta_t)
    vel2 = slider_velocity(L2, L3, t2 + omega * delta_t * 180 / m.pi, omega, delta_t)
    if np.isnan(vel1) or np.isnan(vel2):
        return float('nan')
    return (vel2 - vel1) / delta_t


if __name__ == "__main__":
    L2, L3 = 0.031, 0.023
    omega = 1  # rad/s
    angles = np.linspace(140, 179, 500)

    velocities = [slider_velocity(L2, L3, angle, omega) for angle in angles]
    rod_ang_velocities = [rod_angular_velocity(L2, L3, angle, omega) for angle in angles]

    # Mechanical advantage (output force / input force)
    mechanical_advantage = [(omega * L2 / v) if not np.isnan(v) else float('nan') for v in velocities]

    accelerations = [slider_acceleration(L2, L3, angle, omega) for angle in angles]

    # Plot slider velocity
    plt.figure("Slider Velocity")
    plt.plot(angles, velocities, 'b-', label='Slider Velocity (m/s)')
    plt.xlabel("t2 (degrees)")
    plt.ylabel("Slider Velocity (m/s)")
    plt.title(f"Slider Velocity vs Crank Angle (ω={omega:.2f} rad/s)")
    plt.grid(True)
    plt.legend()

    # Plot rod angular velocity
    plt.figure("Rod Angular Velocity")
    plt.plot(angles, rod_ang_velocities, 'g-', label='Rod Angular Velocity (rad/s)')
    plt.xlabel("t2 (degrees)")
    plt.ylabel("Rod Angular Velocity (rad/s)")
    plt.title(f"Rod Angular Velocity vs Crank Angle (ω={omega:.2f} rad/s)")
    plt.grid(True)
    plt.legend()

    # Plot mechanical advantage
    plt.figure("Mechanical Advantage")
    plt.plot(angles, mechanical_advantage, 'm-', label='Mechanical Advantage (ωL2 / v)')
    plt.xlabel("t2 (degrees)")
    plt.ylabel("Mechanical Advantage")
    plt.title("Output Force / Input Force")
    plt.grid(True)
    plt.legend()

    # Plot slider acceleration
    plt.figure("Slider Acceleration")
    plt.plot(angles, accelerations, 'r-', label='Slider Acceleration (m/s²)')
    plt.xlabel("t2 (degrees)")
    plt.ylabel("Slider Acceleration (m/s²)")
    plt.title("Slider Acceleration vs Crank Angle")
    plt.grid(True)
    plt.legend()

    plt.show()
