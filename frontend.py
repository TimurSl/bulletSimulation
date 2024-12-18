import tkinter as tk
from tkinter import ttk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.pyplot as plt
import numpy as np
from simulation import simulate_bullet_trajectory, calculate_cross_sectional_area


def run_simulation():
    """Запускает симуляцию и обновляет график."""
    try:
        # Получаем параметры из полей ввода
        v0 = float(v0_entry.get())
        theta = float(theta_entry.get()) * (np.pi / 180)  # Перевод градусов в радианы
        phi = float(phi_entry.get()) * (np.pi / 180)
        wx = float(wx_entry.get())
        wy = float(wy_entry.get())
        wz = float(wz_entry.get())
        x0 = float(x0_entry.get())
        y0 = float(y0_entry.get())
        z0 = float(z0_entry.get())
        t_max = float(t_max_entry.get())
        dt = float(dt_entry.get())
        m = float(m_entry.get())
        Cd = float(Cd_entry.get())
        diameter_mm = float(diameter_entry.get())
        A = calculate_cross_sectional_area(diameter_mm)

        # Запуск симуляции
        global simulation_data
        simulation_data = simulate_bullet_trajectory(
            v0=v0, theta=theta, phi=phi,
            wx=wx, wy=wy, wz=wz,
            x0=x0, y0=y0, z0=z0,
            t_max=t_max, dt=dt,
            m=m, Cd=Cd, A=A
        )

        # Обновление графика
        update_plot()
    except Exception as e:
        error_label.config(text=f"Ошибка: {e}")


def update_plot():
    """Обновляет график траектории пули."""
    if simulation_data is None:
        return

    x, y, z = simulation_data

    # Очистка графика
    ax.clear()

    # Построение траектории
    ax.plot(x, y, z, label="Trajectory")
    ax.scatter(x[-1], y[-1], z[-1], color="red", label="Impact Point", s=50)
    ax.set_xlabel("X (м)")
    ax.set_ylabel("Y (м)")
    ax.set_zlabel("Z (м)")
    ax.legend()
    canvas.draw()


def toggle_auto_simulation():
    """Включает или отключает автоматическое обновление при изменении параметров."""
    if auto_simulation_var.get():
        start_auto_simulation()
    else:
        root.after_cancel(auto_simulation_id)


def start_auto_simulation():
    """Автоматически запускает симуляцию каждые 0.01 секунд."""
    run_simulation()
    global auto_simulation_id
    auto_simulation_id = root.after(10, start_auto_simulation)


def on_parameter_change(*args):
    """Запускает симуляцию при изменении параметров, если автообновление активно."""
    if auto_simulation_var.get():
        run_simulation()


# Настройка интерфейса
root = tk.Tk()
root.title("Bullet Trajectory Simulator (Realtime Updates)")

# Ввод параметров
param_frame = tk.Frame(root)
param_frame.pack(side=tk.LEFT, padx=10, pady=10)

ttk.Label(param_frame, text="Параметры симуляции").grid(row=0, column=0, columnspan=2)

params = [
    ("Начальная скорость (v0, м/с):", "v0_entry"),
    ("Угол возвышения (theta, °):", "theta_entry"),
    ("Азимутальный угол (phi, °):", "phi_entry"),
    ("Скорость ветра по X (wx, м/с):", "wx_entry"),
    ("Скорость ветра по Y (wy, м/с):", "wy_entry"),
    ("Скорость ветра по Z (wz, м/с):", "wz_entry"),
    ("Начальная позиция X (x0, м):", "x0_entry"),
    ("Начальная позиция Y (y0, м):", "y0_entry"),
    ("Начальная позиция Z (z0, м):", "z0_entry"),
    ("Время симуляции (t_max, с):", "t_max_entry"),
    ("Шаг времени (dt, с):", "dt_entry"),
    ("Масса пули (m, кг):", "m_entry"),
    ("Коэф. сопротивления (Cd):", "Cd_entry"),
    ("Диаметр пули (мм):", "diameter_entry")
]

entries = {}
for i, (label, var_name) in enumerate(params, start=1):
    ttk.Label(param_frame, text=label).grid(row=i, column=0, sticky=tk.W, pady=2)
    entry = ttk.Entry(param_frame)
    entry.grid(row=i, column=1, pady=2)
    entries[var_name] = entry
    entry.bind("<KeyRelease>", on_parameter_change)  # Вызов симуляции при изменении

# Сохранение ссылок на поля ввода
v0_entry = entries["v0_entry"]
theta_entry = entries["theta_entry"]
phi_entry = entries["phi_entry"]
wx_entry = entries["wx_entry"]
wy_entry = entries["wy_entry"]
wz_entry = entries["wz_entry"]
x0_entry = entries["x0_entry"]
y0_entry = entries["y0_entry"]
z0_entry = entries["z0_entry"]
t_max_entry = entries["t_max_entry"]
dt_entry = entries["dt_entry"]
m_entry = entries["m_entry"]
Cd_entry = entries["Cd_entry"]
diameter_entry = entries["diameter_entry"]

# Установка значений по умолчанию
v0_entry.insert(0, "792.5")
theta_entry.insert(0, "0")
phi_entry.insert(0, "45")
wx_entry.insert(0, "-2.6")
wy_entry.insert(0, "-3.1")
wz_entry.insert(0, "0")
x0_entry.insert(0, "0")
y0_entry.insert(0, "0")
z0_entry.insert(0, "66")
t_max_entry.insert(0, "10")
dt_entry.insert(0, "0.01")
m_entry.insert(0, "0.04277")
Cd_entry.insert(0, "0.686")
diameter_entry.insert(0, "12.954")

# Галочка для автообновления
auto_simulation_var = tk.BooleanVar(value=False)
auto_simulation_id = None
auto_simulation_checkbox = ttk.Checkbutton(
    param_frame, text="Реагировать на изменения и запускать каждые 0.01 секунд",
    variable=auto_simulation_var, command=toggle_auto_simulation
)
auto_simulation_checkbox.grid(row=len(params) + 1, column=0, columnspan=2, pady=10)

# Кнопка запуска симуляции
run_button = ttk.Button(param_frame, text="Запустить симуляцию", command=run_simulation)
run_button.grid(row=len(params) + 2, column=0, columnspan=2, pady=10)

# Метка для ошибок
error_label = ttk.Label(param_frame, text="", foreground="red")
error_label.grid(row=len(params) + 3, column=0, columnspan=2)

# График
fig = plt.figure(figsize=(7, 5))
ax = fig.add_subplot(111, projection="3d")
canvas = FigureCanvasTkAgg(fig, master=root)
canvas.get_tk_widget().pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)

# Глобальные переменные для симуляции
simulation_data = None

root.mainloop()
