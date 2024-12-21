import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PySide6.QtCore import QTimer
from PySide6.QtCore import Qt
from PySide6.QtWidgets import (QApplication, QMainWindow, QVBoxLayout, QHBoxLayout, QWidget, QLabel, QLineEdit,
                               QPushButton, QComboBox, QSpinBox, QSplitter, QScrollArea, QGroupBox)
from PySide6.QtWidgets import QSlider
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas

from simulation import simulate_bullet_trajectory, calculate_cross_sectional_area
import chardet

def load_ammo_presets_from_csv(file_path):
    """
    Загрузить данные о боеприпасах из CSV в формат AMMO_PRESETS с использованием pandas.

    Parameters:
        file_path (str): Путь к CSV-файлу.

    Returns:
        dict: Структура AMMO_PRESETS.
    """
    try:
        # Определяем кодировку файла

        with open(file_path, 'rb') as f:
            result = chardet.detect(f.read())
            encoding = result['encoding']

        # Чтение CSV с помощью pandas с найденной кодировкой
        df = pd.read_csv(file_path, encoding=encoding, dtype={'R': 'float'})
        ammo_presets = {
            row['Type']: {
                "Cd": float(row["Cd"]),
                "diameter_mm": float(row["diameter_mm"]),
                "m": float(row["m"]),
                "v0": int(row["v0"]),
                "R": float(row["R"]) if pd.notna(row["R"]) else None
            } for _, row in df.iterrows()
        }
        return ammo_presets
    except Exception as e:
        print(f"Ошибка при загрузке CSV: {e}")
        return {}




if os.path.exists("ammo_presets.csv"):
    AMMO_PRESETS = load_ammo_presets_from_csv("ammo_presets.csv")
else:
    AMMO_PRESETS = {
    # .50 BMG
    ".50 BMG | 647 gr Speer (G1)": {
        "Cd": 0.686,
        "diameter_mm": 12.98,
        "m": 0.042,
        "v0": 928,
        "R": 0.381,
    },
    ".50 BMG | 647 gr Speer (G7)": {
        "Cd": 0.351,
        "diameter_mm": 12.98,
        "m": 0.042,
        "v0": 928,
        "R": 0.381,

    },
    ".50 BMG | 655 gr ADI (G1)": {
        "Cd": 0.686,
        "diameter_mm": 12.98,
        "m": 0.042,
        "v0": 923,
        "R": 0.381,
    },
    ".50 BMG | 655 gr ADI (G7)": {
        "Cd": 0.351,
        "diameter_mm": 12.98,
        "m": 0.042,
        "v0": 923,
        "R": 0.381,
    },
    ".50 BMG | 700 gr Barnes (G1)": {
        "Cd": 0.686,
        "diameter_mm": 12.98,
        "m": 0.045,
        "v0": 908,
        "R": 0.381,
    },
    ".50 BMG | 700 gr Barnes (G7)": {
        "Cd": 0.351,
        "diameter_mm": 12.98,
        "m": 0.045,
        "R": 0.381,
        "v0": 908,
    },
    ".50 BMG | 750 gr Hornady (G1)": {
        "Cd": 0.686,
        "diameter_mm": 12.98,
        "m": 0.049,
        "R": 0.381,
        "v0": 860,
    },
    ".50 BMG | 750 gr Hornady (G7)": {
        "Cd": 0.351,
        "diameter_mm": 12.98,
        "m": 0.049,
        "R": 0.381,
        "v0": 860,
    },
    ".50 BMG | 800 gr Barnes (G1)": {
        "Cd": 0.686,
        "diameter_mm": 12.98,
        "m": 0.052,
        "R": 0.381,
        "v0": 882
    },
    ".50 BMG | 800 gr Barnes (G7)": {
        "Cd": 0.351,
        "diameter_mm": 12.98,
        "R": 0.381,
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

scrollbar_style = """
        QScrollBar:vertical {
            width: 10px;
        }
        QScrollBar:horizontal {
            height: 10px;
        }
        """


class BulletTrajectorySimulator(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Bullet Trajectory Simulator (Presets)")

        self.main_widget = QWidget()
        self.setCentralWidget(self.main_widget)
        self.layout = QHBoxLayout()
        self.main_widget.setLayout(self.layout)

        splitter = QSplitter()
        self.layout.addWidget(splitter)

        param_widget = QWidget()
        param_widget.setLayout(QVBoxLayout())
        splitter.addWidget(param_widget)
        self.param_layout = param_widget.layout()

        # Scroll area for parameters
        self.scroll_area = QScrollArea()
        self.scroll_area.setWidgetResizable(True)
        self.scroll_content = QWidget()
        self.scroll_layout = QVBoxLayout()
        self.scroll_content.setLayout(self.scroll_layout)
        self.scroll_area.setWidget(self.scroll_content)

        # Apply scrollbar style
        self.scroll_area.setStyleSheet(scrollbar_style)

        self.param_layout.addWidget(self.scroll_area)



        # Group box for presets
        self.preset_group = QGroupBox("Пресеты")
        self.preset_layout = QVBoxLayout()
        self.preset_group.setLayout(self.preset_layout)
        self.scroll_layout.addWidget(self.preset_group)

        # Group box for simulation settings
        self.simulation_group = QGroupBox("Настройки симуляции")
        self.simulation_layout = QVBoxLayout()
        self.simulation_group.setLayout(self.simulation_layout)
        self.scroll_layout.addWidget(self.simulation_group)

        # Group box for scales
        self.scale_group = QGroupBox("Масштабы")
        self.scale_layout = QVBoxLayout()
        self.scale_group.setLayout(self.scale_layout)
        self.scroll_layout.addWidget(self.scale_group)

        # Правая часть: график
        self.figure, self.ax = plt.subplots(subplot_kw={"projection": "3d"})
        self.canvas = FigureCanvas(self.figure)
        splitter.addWidget(self.canvas)

        # Устанавливаем начальные размеры
        splitter.setSizes([300, 700])  # Левая часть (параметры): 300px, правая часть (график): 700px

        self.default_limits = {"x": None, "y": None, "z": None}

        # Ammo presets
        self.ammo_label = QLabel("Выберите боеприпас:")
        self.preset_layout.addWidget(self.ammo_label)
        self.ammo_dropdown = QComboBox()
        self.ammo_dropdown.addItems(AMMO_PRESETS.keys())
        self.preset_layout.addWidget(self.ammo_dropdown)
        self.ammo_button = QPushButton("Применить боеприпас")
        self.ammo_button.clicked.connect(self.apply_ammo_preset)
        self.preset_layout.addWidget(self.ammo_button)

        # Environment presets
        self.env_label = QLabel("Выберите окружение:")
        self.preset_layout.addWidget(self.env_label)
        self.env_dropdown = QComboBox()
        self.env_dropdown.addItems(ENVIRONMENT_PRESETS.keys())
        self.preset_layout.addWidget(self.env_dropdown)
        self.env_button = QPushButton("Применить окружение")
        self.env_button.clicked.connect(self.apply_environment_preset)
        self.preset_layout.addWidget(self.env_button)

        # Simulation parameters
        self.entries = {}
        self.params = [
            ("Начальная скорость (v0, м/с):", "v0", 900),
            ("Угол возвышения (theta, °):", "theta", 0),
            ("Азимутальный угол (phi, °):", "phi", 0),
            ("Скорость ветра по X (wx, м/с):", "wx", 0),
            ("Скорость ветра по Y (wy, м/с):", "wy", 0),
            ("Скорость ветра по Z (wz, м/с):", "wz", 0),
            ("Начальная позиция X (x0, м):", "x0", 0),
            ("Начальная позиция Y (y0, м):", "y0", 0),
            ("Начальная позиция Z (z0, м):", "z0", 1.5),
            ("Время симуляции (t_max, с):", "t_max", 10),
            ("Шаг времени (dt, с):", "dt", 0.01),
            ("Масса пули (m, кг):", "m", 0.045),
            ("Шаг нарезов снаряда (R, м):", "R", 0.3),
            ("Коэф. спин-дрейфа (C_spin):", "C_spin", 1e-6),
            ("Направление вращения нарезов (rifling_sign):", "rifling_sign", 1),
            ("Коэф. сопротивления (Cd):", "Cd", 0.5),
            ("Коэф. Магнуса (Cl, кг·м^2):", "Cl", 1e-6),
            ("Кросс-секция снаряда (A, м^2):", "A", 0),
            ("Диаметр снаряда (мм):", "diameter_mm", 12.7),
            ("Температура воздуха (T, K)", "T_kelvin", 288.15),
            ("Влажность (%)", "humidity", 10),
            ("Широта (latitude, °)", "latitude", 45),
            ("Ускорение свободного падения (g, м/с^2)", "g", 9.81)
        ]

        for label, key, default_value in self.params:
            param_label = QLabel(label)
            self.simulation_layout.addWidget(param_label)
            param_entry = QLineEdit()
            param_entry.setText(str(default_value))
            self.simulation_layout.addWidget(param_entry)
            self.entries[key] = param_entry


        self.scale_sliders = {}
        for axis in ["x", "y", "z"]:
            label = QLabel(f"Масштаб {axis}:")
            self.scale_layout.addWidget(label)

            slider = QSlider(Qt.Horizontal)
            slider.setMinimum(5)  # 1/4 от среднего значения
            slider.setMaximum(1000)  # 4x среднего значения
            slider.setValue(100)  # Исходное значение (1x)
            slider.valueChanged.connect(self.update_axis_scale)
            self.scale_layout.addWidget(slider)
            self.scale_sliders[axis] = slider


        # Общий масштаб
        self.global_scale_label = QLabel("Общий масштаб:")
        self.scale_layout.addWidget(self.global_scale_label)
        self.global_scale_slider = QSlider(Qt.Horizontal)
        self.global_scale_slider.setMinimum(5)
        self.global_scale_slider.setMaximum(1000)
        self.global_scale_slider.setValue(100)
        self.global_scale_slider.valueChanged.connect(self.update_global_scale)
        self.scale_layout.addWidget(self.global_scale_slider)

        # Кнопка сброса масштаба
        self.reset_button = QPushButton("Сбросить масштаб")
        self.reset_button.clicked.connect(self.reset_scale)
        self.scale_layout.addWidget(self.reset_button)

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
            if preset["R"] is not None:
                self.entries["R"].setText(str(preset["R"]))

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

            if params["A"] == "" or params["A"] == 0:
                params["A"] = calculate_cross_sectional_area(params["diameter_mm"])

            # Run simulation
            x, y, z, deviation = simulate_bullet_trajectory(
                v0=params["v0"], theta=params["theta"], phi=params["phi"],
                wx=params["wx"], wy=params["wy"], wz=params["wz"],
                x0=params["x0"], y0=params["y0"], z0=params["z0"],
                t_max=params["t_max"], dt=params["dt"],
                m=params["m"], Cd=params["Cd"], A=params["A"], T_kelvin=params["T_kelvin"],
                humidity=params["humidity"], latitude=params["latitude"], g=params["g"], R=params["R"],
                C_m=params["Cl"],
                C_spin=params["C_spin"], rifling_sign=params["rifling_sign"]
            )

            # Update plot
            self.ax.clear()
            self.ax.plot(x, y, z, label="Trajectory", color="blue")
            self.ax.scatter(x[0], y[0], z[0], color="green", label="Shot Origin", s=50)
            self.ax.scatter(x[-1], y[-1], z[-1], color="red", label="Impact Point", s=50)

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

            # Remove previous deviation widgets if they exist
            if hasattr(self, "deviation_canvas"):
                self.layout.removeWidget(self.deviation_canvas)
                self.deviation_canvas.deleteLater()
                self.deviation_canvas = None

            if hasattr(self, "deviation_label"):
                self.layout.removeWidget(self.deviation_label)
                self.deviation_label.deleteLater()
                self.deviation_label = None

            # Add 2D deviation plot
            self.deviation_fig, self.deviation_ax = plt.subplots()

            # Deviation is based on 3D trajectory
            self.deviation_ax.scatter(0, 0, color="green", label="Target")
            self.deviation_ax.scatter(deviation[1], 0, color="red", label="Impact")
            self.deviation_ax.annotate(f"Deviation: {deviation[1]:.4f} m", (deviation[1], 0), textcoords="offset points", xytext=(0, 10), ha='center')

            # Adjust plot limits to ensure square shape
            max_deviation = max(abs(deviation[1]), 1)  # Ensure at least 1 meter square
            self.deviation_ax.set_xlim([-max_deviation, max_deviation])
            self.deviation_ax.set_ylim([-max_deviation, max_deviation])
            self.deviation_ax.set_aspect('equal', 'box')

            # Add labels and legend
            self.deviation_ax.set_xlabel("Deviation in Y (m)")
            self.deviation_ax.set_ylabel("Deviation in X (m)")
            self.deviation_ax.legend()

            # Add numerical deviation display
            #self.deviation_label = QLabel(f"Deviation: Y = {deviation[1]:.2f} m")
            #self.layout.addWidget(self.deviation_label)

            # Create a new canvas for the deviation plot
            self.deviation_canvas = FigureCanvas(self.deviation_fig)
            self.layout.addWidget(self.deviation_canvas)

            self.update_axis_scale()

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
