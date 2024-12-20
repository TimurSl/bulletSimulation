import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def calculate_cross_sectional_area(diameter_mm):
    diameter_m = diameter_mm / 1000  # Convert mm to meters
    return np.pi * (diameter_m / 2) ** 2

def air_density_with_humidity(h, T_kelvin, humidity, rho0=1.225):
    """
    Calculate air density based on altitude, temperature, and humidity.

    Parameters:
        h (float): Altitude in meters.
        T_kelvin (float): Temperature in Kelvin.
        humidity (float): Humidity in percentage (0-100).
        rho0 (float): Air density at sea level (default: 1.225 kg/m^3).

    Returns:
        float: Air density at the given altitude, temperature, and humidity.
    """
    T0 = 288.15  # Standard temperature at sea level in Kelvin (15°C)
    water_vapor_effect = 1 - (humidity / 100) * 0.015
    return max(rho0 * (1 - 2.2e-5 * h) * (T0 / T_kelvin) * water_vapor_effect, 0.0)

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

def calculate_deviation_with_angles(x, y, z, theta, phi, x0, y0, z0):
    """
    Calculate the deviation of the bullet impact from the intended trajectory line based on azimuth and elevation angles.

    Parameters:
        x (float): Final X-coordinate of the bullet impact.
        y (float): Final Y-coordinate of the bullet impact.
        z (float): Final Z-coordinate of the bullet impact.
        theta (float): Elevation angle in radians.
        phi (float): Azimuthal angle in radians.
        x0 (float): Initial X-coordinate of the shot.
        y0 (float): Initial Y-coordinate of the shot.
        z0 (float): Initial Z-coordinate of the shot.

    Returns:
        tuple: Deviation in the X, Y plane and along the trajectory line.
    """
    # Calculate the ideal impact point based on the trajectory angles
    direction_vector = np.array([np.cos(theta) * np.cos(phi), np.cos(theta) * np.sin(phi), np.sin(theta)])
    point_to_target_vector = np.array([x - x0, y - y0, z - z0])

    # Project point-to-target vector onto direction vector
    projection_length = np.dot(point_to_target_vector, direction_vector)
    projected_point = projection_length * direction_vector

    # Deviation is the vector difference
    deviation_vector = point_to_target_vector - projected_point
    return deviation_vector[0], deviation_vector[1], np.linalg.norm(deviation_vector)

def simulate_bullet_trajectory(v0, theta, phi, wx, wy, wz, humidity, x0=0, y0=0, z0=350, dt=0.01, t_max=20.0, m=0.045, Cd=0.5,
                               A=np.pi * (12.7e-3 / 2) ** 2, rho0=1.225, g=9.81, T_kelvin=288.15, latitude=45.0):
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
        latitude (float): Latitude of the location (default: 45.0 degrees).
        humidity (float): Humidity in percentage (0-100).

    Returns:
        tuple: Four lists containing the X, Y, and Z coordinates and deviation of the bullet's trajectory.
    """

    omega = 7.2921e-5
    lat_rad = np.radians(latitude)

    steps = int(t_max / dt)
    vx0 = v0 * np.cos(theta) * np.cos(phi)
    vy0 = v0 * np.cos(theta) * np.sin(phi)
    vz0 = v0 * np.sin(theta)

    x, y, z = [x0], [y0], [z0]
    vx, vy, vz = vx0, vy0, vz0

    for i in range(steps):
        h = z[-1]
        rho = air_density_with_humidity(h, T_kelvin, humidity, rho0)

        v_rel_x = vx - wx
        v_rel_y = vy - wy
        v_rel_z = vz - wz

        v_rel = np.sqrt(v_rel_x ** 2 + v_rel_y ** 2 + v_rel_z ** 2)

        Fx = -0.5 * Cd * A * rho * v_rel * v_rel_x
        Fy = -0.5 * Cd * A * rho * v_rel * v_rel_y
        Fz = -0.5 * Cd * A * rho * v_rel * v_rel_z - m * g

        # Coriolis force
        Fx_coriolis = -2 * m * omega * (vy * np.sin(lat_rad) - vz * np.cos(lat_rad))
        Fy_coriolis = 2 * m * omega * (vx * np.sin(lat_rad))
        Fz_coriolis = 2 * m * omega * (vx * np.cos(lat_rad))

        Fx_final = Fx + Fx_coriolis
        Fy_final = Fy + Fy_coriolis
        Fz_final = Fz - m * g + Fz_coriolis

        vx += (Fx_final / m) * dt
        vy += (Fy_final / m) * dt
        vz += (Fz_final / m) * dt

        new_x = x[-1] + vx * dt
        new_y = y[-1] + vy * dt
        new_z = z[-1] + vz * dt

        if new_z <= 0:
            break

        x.append(new_x)
        y.append(new_y)
        z.append(new_z)

    deviation = calculate_deviation_with_angles(x[-1], y[-1], z[-1], theta, phi, x0, y0, z0)
    return x, y, z, deviation
