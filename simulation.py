import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def calculate_cross_sectional_area(diameter_mm):
    """
    Calculate the cross-sectional area of the bullet.

    Parameters:
        diameter_mm (float): Diameter of the bullet in millimeters.

    Returns:
        float: Cross-sectional area in square meters.
    """
    diameter_m = diameter_mm / 1000  # Convert mm to meters
    return np.pi * (diameter_m / 2) ** 2

def air_density(h, rho0=1.225):
    """
    Calculate air density based on altitude.

    Parameters:
        h (float): Altitude in meters.
        rho0 (float): Air density at sea level (default: 1.225 kg/m^3).

    Returns:
        float: Air density at the given altitude.
    """
    return max(rho0 * (1 - 2.2e-5 * h), 0.0)

def simulate_bullet_trajectory(v0, theta, phi, wx, wy, wz, x0=0, y0=0, z0=350, dt=0.01, t_max=20.0, m=0.045, Cd=0.5,
                               A=np.pi * (12.7e-3 / 2) ** 2, rho0=1.225, g=9.81):
    """
    Simulate the trajectory of a bullet considering wind and air resistance.

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

    Returns:
        tuple: Three lists containing the X, Y, and Z coordinates of the bullet's trajectory.
    """
    steps = int(t_max / dt)

    # Initial velocities
    vx0 = v0 * np.cos(theta) * np.cos(phi)
    vy0 = v0 * np.cos(theta) * np.sin(phi)
    vz0 = v0 * np.sin(theta)

    # Trajectory lists
    x, y, z = [x0], [y0], [z0]
    vx, vy, vz = vx0, vy0, vz0

    for i in range(steps):
        h = z[-1]
        rho = air_density(h, rho0)

        # Relative velocity (including wind)
        v_rel_x = vx - wx
        v_rel_y = vy - wy
        v_rel_z = vz - wz

        v_rel = np.sqrt(v_rel_x ** 2 + v_rel_y ** 2 + v_rel_z ** 2)

        # Forces
        Fx = -0.5 * Cd * A * rho * v_rel * v_rel_x
        Fy = -0.5 * Cd * A * rho * v_rel * v_rel_y
        Fz = -0.5 * Cd * A * rho * v_rel * v_rel_z - m * g

        # Update velocities
        vx += (Fx / m) * dt
        vy += (Fy / m) * dt
        vz += (Fz / m) * dt

        # Update positions
        new_x = x[-1] + vx * dt
        new_y = y[-1] + vy * dt
        new_z = z[-1] + vz * dt

        if new_z <= 0:
            break

        x.append(new_x)
        y.append(new_y)
        z.append(new_z)

    return x, y, z

# # Example usage
# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
#
# # Parameters for the first simulation
# x, y, z = simulate_bullet_trajectory(
#     m=0.04277, # kg (0.04277 kg = 42.77 g)
#     v0=792.5, # m/s
#     theta=-0.0 * np.pi, # radians
#     phi=0.434 * np.pi, # radians
#     wx=-2.6482849519085416, # m/s currently in London
#     wy=-3.1561031056501894, # m/s currently in London
#     wz=0.0, # m/s (no wind in the Z-axis)
#     x0=0, # m
#     y0=0, # m
#     z0=66, # m height of the BIG BEN from sea level
#     t_max=10, # s (simulation time)
#     Cd=0.686, # drag coefficient for a bullet (0.686)
#     A=calculate_cross_sectional_area(12.954), # m^2 (cross-sectional area of a 12.954 mm bullet)
#     dt=0.001 # s (time step)
# )
# ax.plot(x, y, z, label='Trajectory 1')
#
# # Add a target cube (10x10x10 meters)
# target_center = [250, 200, 10]  # Center of the target cube
# cube_size = 10
#
# # Plot edges of the cube
# for i in range(2):
#     for j in range(2):
#         for k in range(2):
#             x_corner = target_center[0] + cube_size * (-0.5 if i == 0 else 0.5)
#             y_corner = target_center[1] + cube_size * (-0.5 if j == 0 else 0.5)
#             z_corner = target_center[2] + cube_size * (-0.5 if k == 0 else 0.5)
#             ax.scatter(x_corner, y_corner, z_corner, color="blue")
#
# # Add a red dot where the bullet lands
# if len(x) > 0:
#     ax.scatter(x[-1], y[-1], z[-1], color="red", label="Impact Point", s=50)
#
# # Axis labels and legend
# ax.set_xlabel('X (м)')
# ax.set_ylabel('Y (м)')
# ax.set_zlabel('Z (м)')
# ax.legend()
# plt.show()
