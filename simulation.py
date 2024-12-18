import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def calculate_cross_sectional_area(diameter_mm):
    diameter_m = diameter_mm / 1000  # Convert mm to meters
    return np.pi * (diameter_m / 2) ** 2

def air_density(h, rho0=1.225):
    return max(rho0 * (1 - 2.2e-5 * h), 0.0)

def simulate_bullet_trajectory(v0, theta, phi, wx, wy, wz, x0=0, y0=0, z0=350, dt=0.01, t_max=20.0, m=0.045, Cd=0.5,
                               A=np.pi * (12.7e-3 / 2) ** 2, rho0=1.225, g=9.81):
    steps = int(t_max / dt)
    vx0 = v0 * np.cos(theta) * np.cos(phi)
    vy0 = v0 * np.cos(theta) * np.sin(phi)
    vz0 = v0 * np.sin(theta)

    x, y, z = [x0], [y0], [z0]
    vx, vy, vz = vx0, vy0, vz0

    for i in range(steps):
        h = z[-1]
        rho = air_density(h, rho0)

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

# Interactive plot
def on_hover(event, x, y, z):
    if event.inaxes is not None:
        closest_index = np.argmin(
            np.sqrt((x - event.xdata)**2 + (y - event.ydata)**2)
        )
        print(f"Hovered position: X={x[closest_index]:.2f}, Y={y[closest_index]:.2f}, Z={z[closest_index]:.2f}")

# # Example usage
# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
#
# # Parameters for the first simulation
# x, y, z = simulate_bullet_trajectory(
#     m=0.04277,  # kg (0.04277 kg = 42.77 g)
#     v0=792.5,  # m/s
#     theta=0.05 * np.pi,  # radians
#     phi=0.25 * np.pi,  # radians
#     wx=-3.0, wy=-4.0, wz=0.0,
#     x0=0, y0=0, z0=66,  # Starting position
#     t_max=10, dt=0.001
# )
#
# # Plot trajectory
# ax.plot(x, y, z, label="Trajectory", color="blue")
#
# # Mark start and impact points
# ax.scatter(x[0], y[0], z[0], color="green", label="Shot Origin", s=50)
# ax.scatter(x[-1], y[-1], z[-1], color="red", label="Impact Point", s=50)
#
# # Add target cube
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
#             ax.scatter(x_corner, y_corner, z_corner, color="cyan")
#
# # Axis labels and legend
# ax.set_xlabel('X (м)')
# ax.set_ylabel('Y (м)')
# ax.set_zlabel('Z (м)')
# ax.legend()
#
# # Connect the hover event
# fig.canvas.mpl_connect("motion_notify_event", lambda event: on_hover(event, np.array(x), np.array(y), np.array(z)))
#
# plt.show()
