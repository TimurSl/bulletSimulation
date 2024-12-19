from PySide6.QtWidgets import (QApplication, QMainWindow, QVBoxLayout, QHBoxLayout, QWidget, QLabel, QLineEdit, QPushButton, QComboBox, QSpinBox, QSplitter)
from PySide6.QtCore import QTimer
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
import matplotlib.pyplot as plt
import numpy as np
from simulation import simulate_bullet_trajectory, calculate_cross_sectional_area

from PySide6.QtWidgets import QSlider
from PySide6.QtCore import Qt

# Пресеты боеприпасов
AMMO_PRESETS = {
    # .50 BMG
    ".50 BMG | 647 gr Speer (G1)": {
        "Cd": 0.686,
        "diameter_mm": 12.98,
        "m": 0.042,
        "v0": 928
    },
    ".50 BMG | 647 gr Speer (G7)": {
        "Cd": 0.351,
        "diameter_mm": 12.98,
        "m": 0.042,
        "v0": 928
    },
    ".50 BMG | 655 gr ADI (G1)": {
        "Cd": 0.686,
        "diameter_mm": 12.98,
        "m": 0.042,
        "v0": 923
    },
    ".50 BMG | 655 gr ADI (G7)": {
        "Cd": 0.351,
        "diameter_mm": 12.98,
        "m": 0.042,
        "v0": 923
    },
    ".50 BMG | 700 gr Barnes (G1)": {
        "Cd": 0.686,
        "diameter_mm": 12.98,
        "m": 0.045,
        "v0": 908
    },
    ".50 BMG | 700 gr Barnes (G7)": {
        "Cd": 0.351,
        "diameter_mm": 12.98,
        "m": 0.045,
        "v0": 908
    },
    ".50 BMG | 750 gr Hornady (G1)": {
        "Cd": 0.686,
        "diameter_mm": 12.98,
        "m": 0.049,
        "v0": 860
    },
    ".50 BMG | 750 gr Hornady (G7)": {
        "Cd": 0.351,
        "diameter_mm": 12.98,
        "m": 0.049,
        "v0": 860
    },
    ".50 BMG | 800 gr Barnes (G1)": {
        "Cd": 0.686,
        "diameter_mm": 12.98,
        "m": 0.052,
        "v0": 882
    },
    ".50 BMG | 800 gr Barnes (G7)": {
        "Cd": 0.351,
        "diameter_mm": 12.98,
        "m": 0.052,
        "v0": 882
    },

    # 7.92x57 Mauser
    "7.92x57 Mauser | 11.3 g HPBT (G1)": {
        "Cd": 0.290,
        "diameter_mm": 7.92,
        "m": 0.0113,
        "v0": 797
    },
    "7.92x57 Mauser | 11.3 g HPBT (G7)": {
        "Cd": 0.148,
        "diameter_mm": 7.92,
        "m": 0.0113,
        "v0": 797
    },
    "7.92x57 Mauser | 11.7 g FMJ (G1)": {
        "Cd": 0.290,
        "diameter_mm": 7.92,
        "m": 0.0117,
        "v0": 786
    },
    "7.92x57 Mauser | 11.7 g FMJ (G7)": {
        "Cd": 0.148,
        "diameter_mm": 7.92,
        "m": 0.0117,
        "v0": 786
    },
    "7.92x57 Mauser | 11.7 g SP (800 m/s, G1)": {
        "Cd": 0.290,
        "diameter_mm": 7.92,
        "m": 0.0117,
        "v0": 800
    },
    "7.92x57 Mauser | 11.7 g SP (800 m/s, G7)": {
        "Cd": 0.148,
        "diameter_mm": 7.92,
        "m": 0.0117,
        "v0": 800
    },
    "7.92x57 Mauser | 9.7 g FMJ (G1)": {
        "Cd": 0.290,
        "diameter_mm": 7.92,
        "m": 0.0097,
        "v0": 865
    },
    "7.92x57 Mauser | 9.7 g FMJ (G7)": {
        "Cd": 0.148,
        "diameter_mm": 7.92,
        "m": 0.0097,
        "v0": 865
    },
    "7.92x57 Mauser | 11.7 g SP (805 m/s, G1)": {
        "Cd": 0.290,
        "diameter_mm": 7.92,
        "m": 0.0117,
        "v0": 805
    },
    "7.92x57 Mauser | 11.7 g SP (805 m/s, G7)": {
        "Cd": 0.148,
        "diameter_mm": 7.92,
        "m": 0.0117,
        "v0": 805
    }
}


# Пресеты окружения
ENVIRONMENT_PRESETS = {
    "Calm": {"wx": 0, "wy": 0, "wz": 0},
    "Windy": {"wx": 10, "wy": -5, "wz": 0},
    "Rain": {"wx": 2, "wy": 2, "wz": 0},
    "Storm": {"wx": 20, "wy": 15, "wz": 5}
}

class BulletTrajectorySimulator(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Bullet Trajectory Simulator (Presets)")

        # Main layout
        self.main_widget = QWidget()
        self.setCentralWidget(self.main_widget)
        self.layout = QHBoxLayout()
        self.main_widget.setLayout(self.layout)

        # Добавляем QSplitter
        splitter = QSplitter()
        self.layout.addWidget(splitter)

        # Левая часть: параметры
        param_widget = QWidget()
        param_widget.setLayout(QVBoxLayout())
        splitter.addWidget(param_widget)
        self.param_layout = param_widget.layout()  # Сохраняем ссылку на макет параметров

        # Правая часть: график
        self.figure, self.ax = plt.subplots(subplot_kw={"projection": "3d"})
        self.canvas = FigureCanvas(self.figure)
        splitter.addWidget(self.canvas)

        # Устанавливаем начальные размеры
        splitter.setSizes([300, 700])  # Левая часть (параметры): 300px, правая часть (график): 700px

        self.default_limits = {"x": None, "y": None, "z": None}

        # Ammo presets
        self.ammo_label = QLabel("Выберите боеприпас:")
        self.param_layout.addWidget(self.ammo_label)
        self.ammo_dropdown = QComboBox()
        self.ammo_dropdown.addItems(AMMO_PRESETS.keys())
        self.param_layout.addWidget(self.ammo_dropdown)
        self.ammo_button = QPushButton("Применить боеприпас")
        self.ammo_button.clicked.connect(self.apply_ammo_preset)
        self.param_layout.addWidget(self.ammo_button)

        # Environment presets
        self.env_label = QLabel("Выберите окружение:")
        self.param_layout.addWidget(self.env_label)
        self.env_dropdown = QComboBox()
        self.env_dropdown.addItems(ENVIRONMENT_PRESETS.keys())
        self.param_layout.addWidget(self.env_dropdown)
        self.env_button = QPushButton("Применить окружение")
        self.env_button.clicked.connect(self.apply_environment_preset)
        self.param_layout.addWidget(self.env_button)

        # Simulation parameters
        self.entries = {}
        self.params = [
            ("Начальная скорость (v0, м/с):", "v0"),
            ("Угол возвышения (theta, °):", "theta"),
            ("Азимутальный угол (phi, °):", "phi"),
            ("Скорость ветра по X (wx, м/с):", "wx"),
            ("Скорость ветра по Y (wy, м/с):", "wy"),
            ("Скорость ветра по Z (wz, м/с):", "wz"),
            ("Начальная позиция X (x0, м):", "x0"),
            ("Начальная позиция Y (y0, м):", "y0"),
            ("Начальная позиция Z (z0, м):", "z0"),
            ("Время симуляции (t_max, с):", "t_max"),
            ("Шаг времени (dt, с):", "dt"),
            ("Масса пули (m, кг):", "m"),
            ("Коэф. сопротивления (Cd):", "Cd"),
            ("Диаметр пули (мм):", "diameter_mm"),
            ("Температура воздуха (T, K)", "T_kelvin")
        ]

        for label, key in self.params:
            param_label = QLabel(label)
            self.param_layout.addWidget(param_label)
            param_entry = QLineEdit()
            self.param_layout.addWidget(param_entry)
            self.entries[key] = param_entry

        self.scale_label = QLabel("Масштаб графика:")
        self.param_layout.addWidget(self.scale_label)

        self.scale_sliders = {}
        for axis in ["x", "y", "z"]:
            label = QLabel(f"Масштаб {axis}:")
            self.param_layout.addWidget(label)

            slider = QSlider(Qt.Horizontal)
            slider.setMinimum(5)  # 1/4 от среднего значения
            slider.setMaximum(1000)  # 4x среднего значения
            slider.setValue(100)  # Исходное значение (1x)
            slider.valueChanged.connect(self.update_axis_scale)
            self.param_layout.addWidget(slider)
            self.scale_sliders[axis] = slider


        # Общий масштаб
        self.global_scale_label = QLabel("Общий масштаб:")
        self.param_layout.addWidget(self.global_scale_label)
        self.global_scale_slider = QSlider(Qt.Horizontal)
        self.global_scale_slider.setMinimum(5)
        self.global_scale_slider.setMaximum(1000)
        self.global_scale_slider.setValue(100)
        self.global_scale_slider.valueChanged.connect(self.update_global_scale)
        self.param_layout.addWidget(self.global_scale_slider)

        # Кнопка сброса масштаба
        self.reset_button = QPushButton("Сбросить масштаб")
        self.reset_button.clicked.connect(self.reset_scale)
        self.param_layout.addWidget(self.reset_button)

        # Interval and auto-update
        self.interval_label = QLabel("Интервал автообновления (мс):")
        self.param_layout.addWidget(self.interval_label)
        self.interval_spinbox = QSpinBox()
        self.interval_spinbox.setRange(100, 10000)
        self.interval_spinbox.setValue(1000)
        self.param_layout.addWidget(self.interval_spinbox)

        self.auto_update_button = QPushButton("Запустить автообновление")
        self.auto_update_button.setCheckable(True)
        self.auto_update_button.clicked.connect(self.toggle_auto_update)
        self.param_layout.addWidget(self.auto_update_button)



        # Run button
        self.run_button = QPushButton("Запустить симуляцию")
        self.run_button.clicked.connect(self.run_simulation)
        self.param_layout.addWidget(self.run_button)

        # Timer for auto-update
        self.timer = QTimer()
        self.timer.timeout.connect(self.run_simulation)

    def apply_ammo_preset(self):
        ammo_name = self.ammo_dropdown.currentText()
        preset = AMMO_PRESETS.get(ammo_name, {})
        if preset:
            self.entries["Cd"].setText(str(preset["Cd"]))
            self.entries["diameter_mm"].setText(str(preset["diameter_mm"]))
            self.entries["m"].setText(str(preset["m"]))
            self.entries["v0"].setText(str(preset["v0"]))

    def apply_environment_preset(self):
        env_name = self.env_dropdown.currentText()
        preset = ENVIRONMENT_PRESETS.get(env_name, {})
        if preset:
            self.entries["wx"].setText(str(preset["wx"]))
            self.entries["wy"].setText(str(preset["wy"]))
            self.entries["wz"].setText(str(preset["wz"]))

    def toggle_auto_update(self):
        if self.auto_update_button.isChecked():
            interval = self.interval_spinbox.value()
            self.timer.start(interval)
            self.auto_update_button.setText("Остановить автообновление")
        else:
            self.timer.stop()
            self.auto_update_button.setText("Запустить автообновление")

    def run_simulation(self):
        try:
            # Get parameters from inputs
            params = {key: float(self.entries[key].text()) for key in self.entries}
            params["theta"] *= np.pi / 180  # Convert degrees to radians
            params["phi"] *= np.pi / 180
            params["A"] = calculate_cross_sectional_area(params["diameter_mm"])

            # Run simulation
            x, y, z = simulate_bullet_trajectory(
                v0=params["v0"], theta=params["theta"], phi=params["phi"],
                wx=params["wx"], wy=params["wy"], wz=params["wz"],
                x0=params["x0"], y0=params["y0"], z0=params["z0"],
                t_max=params["t_max"], dt=params["dt"],
                m=params["m"], Cd=params["Cd"], A=params["A"], T_kelvin=params["T_kelvin"]
            )

            # Update plot
            self.ax.clear()
            self.ax.plot(x, y, z, label="Trajectory", color="blue")
            self.ax.scatter(x[0], y[0], z[0], color="green", label="Shot Origin", s=50)
            self.ax.scatter(x[-1], y[-1], z[-1], color="red", label="Impact Point", s=50)

            x, y, z = simulate_bullet_trajectory(
                v0=params["v0"], theta=params["theta"], phi=params["phi"],
                wx=0, wy=0, wz=0,
                x0=params["x0"], y0=params["y0"], z0=params["z0"],
                t_max=params["t_max"], dt=params["dt"],
                m=params["m"], Cd=params["Cd"], A=params["A"], T_kelvin=params["T_kelvin"]
            )
            self.ax.plot(x, y, z, label="Trajectory (No Wind)", color="green")

            # Save axis limits
            self.default_limits = {
                "x": self.ax.get_xlim(),
                "y": self.ax.get_ylim(),
                "z": self.ax.get_zlim()
            }

            # Axis labels
            self.ax.set_xlabel("X (м)")
            self.ax.set_ylabel("Y (м)")
            self.ax.set_zlabel("Z (м)")
            self.ax.legend()

            self.canvas.draw()
        except Exception as e:
            print(f"Ошибка: {e}")


    def update_axis_scale(self):
        """
        Update the scale of individual axes based on their sliders.
        """
        for axis, slider in self.scale_sliders.items():
            if self.default_limits[axis] is None:
                continue

            scale_factor = slider.value() / 100
            center = (self.default_limits[axis][0] + self.default_limits[axis][1]) / 2
            range_val = (self.default_limits[axis][1] - self.default_limits[axis][0]) * scale_factor
            getattr(self.ax, f"set_{axis}lim")([center - range_val / 2, center + range_val / 2])
        self.canvas.draw()

    def update_global_scale(self):
        """
        Update all axis scales based on the global scale slider.
        """
        global_scale = self.global_scale_slider.value()
        for axis, slider in self.scale_sliders.items():
            slider.setValue(global_scale)

    def reset_scale(self):
        """
        Reset all scales to their default values.
        """
        self.global_scale_slider.setValue(100)
        for slider in self.scale_sliders.values():
            slider.setValue(100)


if __name__ == "__main__":
    app = QApplication([])
    window = BulletTrajectorySimulator()
    window.show()
    app.exec()
