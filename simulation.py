import time

import numpy as np


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


def compute_derivatives(state, params):
    x, y, z, vx, vy, vz = state
    m = params['m']
    Cd = params['Cd']
    A = params['A']
    rho0 = params['rho0']
    g = params['g']
    T_kelvin = params['T_kelvin']
    humidity = params['humidity']
    wx = params['wx']
    wy = params['wy']
    wz = params['wz']
    phi = params['phi']
    theta = params['theta']
    latitude = params['latitude']
    R = params['R']
    C_m = params['C_m']
    C_spin = params['C_spin']
    rifling_sign = params['rifling_sign']

    # Вычисляем плотность воздуха
    h = z
    rho = air_density_with_humidity(h, T_kelvin, humidity, rho0)

    # Относительная скорость
    v_rel_x = vx - wx
    v_rel_y = vy - wy
    v_rel_z = vz - wz

    v_rel_vec = np.array([v_rel_x, v_rel_y, v_rel_z])
    v_rel = np.linalg.norm(v_rel_vec)

    v = np.sqrt(vx**2 + vy**2 + vz**2)

    # Сила сопротивления
    Fx = -0.5 * Cd * A * rho * v_rel * v_rel_x
    Fy = -0.5 * Cd * A * rho * v_rel * v_rel_y
    Fz = -0.5 * Cd * A * rho * v_rel * v_rel_z - m * g

    # Магнус эффект
    if v > 0:
        omega_mag = (2 * np.pi * v) / R
        v_norm = np.linalg.norm([vx, vy, vz])
        if v_norm > 1e-8:
            omega_vector = omega_mag * np.array([vx, vy, vz]) / v_norm
        else:
            omega_vector = np.array([0.0, 0.0, 0.0])
    else:
        omega_vector = np.array([0.0, 0.0, 0.0])

    cos_theta = np.dot(omega_vector, v_rel_vec) / (np.linalg.norm(omega_vector) * v_rel + 1e-8)
    if abs(cos_theta) > 0.99:
        omega_vector += np.array([0.01, -0.01, 0.01])

    if v_rel > 0:
        F_magnus = C_m * 0.5 * rho * A * v_rel**2 * np.cross(omega_vector, v_rel_vec / v_rel)
    else:
        F_magnus = np.array([0, 0, 0])

    Fx_magnus, Fy_magnus, Fz_magnus = F_magnus

    # Сила Кориолиса
    omega_coriolis = 7.2921e-5
    lat_rad = np.radians(latitude)
    Fx_coriolis = -2 * m * omega_coriolis * (vy * np.sin(lat_rad) - vz * np.cos(lat_rad))
    Fy_coriolis =  2 * m * omega_coriolis * (vx * np.sin(lat_rad))
    Fz_coriolis =  2 * m * omega_coriolis * (vx * np.cos(lat_rad))


    # Сила вращения
    Fx_spin = Fy_spin = Fz_spin = 0.0

    if v_rel > 1e-8:  # Ensure nonzero velocity
        # unit vector of velocity
        v_hat = v_rel_vec / v_rel
        z_hat = np.array([0, 0, 1])
        # direction of spin drift
        drift_dir = np.cross(v_hat, z_hat)  # v̂ × ẑ
        drift_dir_norm = np.linalg.norm(drift_dir)
        if drift_dir_norm > 1e-8:
            drift_dir /= drift_dir_norm
        else:
            drift_dir = np.array([0, 0, 0])
        # Force magnitude ~ C_spin * v^2 (tunable)
        F_spin = C_spin * (v_rel ** 2) * rifling_sign
        Fx_spin, Fy_spin, Fz_spin = F_spin * drift_dir


    # Sum all forces
    Fx_final = Fx + Fx_magnus + Fx_coriolis + Fx_spin
    Fy_final = Fy + Fy_magnus + Fy_coriolis + Fy_spin
    Fz_final = Fz + Fz_magnus + Fz_coriolis + Fz_spin

    # Производные
    dxdt = vx
    dydt = vy
    dzdt = vz
    dvxdt = Fx_final / m
    dvydt = Fy_final / m
    dvzdt = Fz_final / m

    return np.array([dxdt, dydt, dzdt, dvxdt, dvydt, dvzdt])

def simulate_bullet_trajectory(v0, theta, phi, wx, wy, wz, humidity, R, C_m,
                                               x0=0, y0=0, z0=350, dt=0.01, t_max=20.0, m=0.045, Cd=0.5,
                                               A=np.pi * (12.7e-3 / 2) ** 2, rho0=1.225, g=9.81, T_kelvin=288.15, latitude=45.0,
                                               C_spin=0.0001, rifling_sign=1):
    vx0 = v0 * np.cos(theta) * np.cos(phi)
    vy0 = v0 * np.cos(theta) * np.sin(phi)
    vz0 = v0 * np.sin(theta)

    state = np.array([x0, y0, z0, vx0, vy0, vz0])

    params = {
        'm': m, 'Cd': Cd, 'A': A, 'rho0': rho0, 'g': g, 'T_kelvin': T_kelvin,
        'humidity': humidity, 'wx': wx, 'wy': wy, 'wz': wz, 'phi': phi, 'theta': theta,
        'latitude': latitude, 'R': R, 'C_m': C_m, 'C_spin': C_spin, 'rifling_sign': rifling_sign
    }

    steps = int(t_max / dt)
    X, Y, Z = [x0], [y0], [z0]

    for i in range(steps):
        deriv = compute_derivatives(state, params)
        # Простой метод Эйлера для интегрирования
        state = state + deriv * dt
        x, y, z, vx, vy, vz = state

        if z <= 0:
            break

        X.append(x)
        Y.append(y)
        Z.append(z)

    deviation = calculate_deviation_with_angles(X[-1], Y[-1], Z[-1], theta, phi, x0, y0, z0)
    return X, Y, Z, deviation
