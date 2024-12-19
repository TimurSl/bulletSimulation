import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def calculate_cross_sectional_area(diameter_mm):
    diameter_m = diameter_mm / 1000  # Convert mm to meters
    return np.pi * (diameter_m / 2) ** 2

def air_density(h, T_kelvin, rho0=1.225):
    """
    Calculate air density based on altitude and temperature.

    Parameters:
        h (float): Altitude in meters.
        T_kelvin (float): Temperature in Kelvin.
        rho0 (float): Air density at sea level (default: 1.225 kg/m^3).

    Returns:
        float: Air density at the given altitude and temperature.
    """
    T0 = 288.15  # Standard temperature at sea level in Kelvin (15°C)
    return max(rho0 * (1 - 2.2e-5 * h) * (T0 / T_kelvin), 0.0)

def calculate_v0(base_v0, barrel_length, standard_length=20.0):
    """
    Calculate the initial velocity (v0) based on the barrel length.

    Parameters:
        base_v0 (float): Base velocity for the standard barrel length (m/s).
        barrel_length (float): Actual barrel length (cm).
        standard_length (float): Standard barrel length for the base velocity (cm).

    Returns:
        float: Adjusted initial velocity.
    """
    return base_v0 * (barrel_length / standard_length)

def simulate_bullet_trajectory(v0, theta, phi, wx, wy, wz, x0=0, y0=0, z0=350, dt=0.01, t_max=20.0, m=0.045, Cd=0.5,
                               A=np.pi * (12.7e-3 / 2) ** 2, rho0=1.225, g=9.81, T_kelvin=288.15):
    """
    Simulate the trajectory of a bullet considering wind, air resistance, and temperature.

    Parameters:
        v0 (float): Initial velocity of the bullet (m/s).
        theta (float): Elevation angle in radians.
        phi (float): Azimuthal angle in radians.
        wx (float): Wind velocity along the X-axis (m/s).
        wy (float): Wind velocity along the Y-axis (m/s).
        wz (float): Wind velocity along the Z-axis (m/s).
        x0 (float): Initial X-coordinate (default: 0).
        y0 (float): Initial Y-coordinate (default: 0).
        z0 (float): Initial Z-coordinate (default: 350).
        dt (float): Time step for the simulation (default: 0.01 s).
        t_max (float): Maximum simulation time (default: 20.0 s).
        m (float): Mass of the bullet (default: 0.045 kg).
        Cd (float): Drag coefficient (default: 0.5).
        A (float): Cross-sectional area of the bullet (default: 0.0001267 m^2).
        rho0 (float): Air density at sea level (default: 1.225 kg/m^3).
        g (float): Gravitational acceleration (default: 9.81 m/s^2).
        T_kelvin (float): Air temperature in Kelvin (default: 288.15 K, 15°C).

    Returns:
        tuple: Three lists containing the X, Y, and Z coordinates of the bullet's trajectory.
    """
    steps = int(t_max / dt)
    vx0 = v0 * np.cos(theta) * np.cos(phi)
    vy0 = v0 * np.cos(theta) * np.sin(phi)
    vz0 = v0 * np.sin(theta)

    x, y, z = [x0], [y0], [z0]
    vx, vy, vz = vx0, vy0, vz0

    for i in range(steps):
        h = z[-1]
        rho = air_density(h, T_kelvin, rho0)

        v_rel_x = vx - wx
        v_rel_y = vy - wy
        v_rel_z = vz - wz

        v_rel = np.sqrt(v_rel_x ** 2 + v_rel_y ** 2 + v_rel_z ** 2)

        Fx = -0.5 * Cd * A * rho * v_rel * v_rel_x
        Fy = -0.5 * Cd * A * rho * v_rel * v_rel_y
        Fz = -0.5 * Cd * A * rho * v_rel * v_rel_z - m * g

        vx += (Fx / m) * dt
        vy += (Fy / m) * dt
        vz += (Fz / m) * dt

        new_x = x[-1] + vx * dt
        new_y = y[-1] + vy * dt
        new_z = z[-1] + vz * dt

        if new_z <= 0:
            break

        x.append(new_x)
        y.append(new_y)
        z.append(new_z)

    return x, y, z