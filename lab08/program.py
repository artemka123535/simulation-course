import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import poisson
from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, 
    QHBoxLayout, QLabel, QPushButton, QDoubleSpinBox, 
    QSpinBox, QFrame, QTextEdit
)
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas

STYLE_SHEET = """
QMainWindow { background-color: #121212; }
QLabel { color: #00ff41; font-family: 'Courier New'; font-weight: bold; }
QGroupBox { border: 1px solid #00ff41; color: #00ff41; margin-top: 10px; }
QPushButton { 
    background-color: #003300; color: #00ff41; 
    border: 1px solid #00ff41; padding: 10px; border-radius: 4px;
}
QPushButton:hover { background-color: #004400; }
QDoubleSpinBox, QSpinBox { 
    background-color: #1a1a1a; color: #00ff41; 
    border: 1px solid #004400; 
}
QTextEdit { background-color: #0a0a0a; color: #00ff41; border: none; font-size: 11px; }
"""

class PoissonEngine:
    @staticmethod
    def simulate(intensity, t_interval, n_trials):
        results = []
        for _ in range(n_trials):
            time = 0
            count = 0
            while True:
                dt = -np.log(1.0 - np.random.uniform()) / intensity
                time += dt
                if time > t_interval:
                    break
                count += 1
            results.append(count)
        return np.array(results)

class LabWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("События на сервере")
        self.resize(1100, 700)
        self.setStyleSheet(STYLE_SHEET)

        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        layout = QHBoxLayout(central_widget)

        control_panel = QVBoxLayout()
        
        control_panel.addWidget(QLabel("ИНТЕНСИВНОСТЬ (λ):"))
        self.sb_lambda = QDoubleSpinBox()
        self.sb_lambda.setRange(0.1, 50.0)
        self.sb_lambda.setValue(5.0)
        control_panel.addWidget(self.sb_lambda)

        control_panel.addWidget(QLabel("ИНТЕРВАЛ (T):"))
        self.sb_t = QDoubleSpinBox()
        self.sb_t.setRange(0.1, 20.0)
        self.sb_t.setValue(2.0)
        control_panel.addWidget(self.sb_t)

        control_panel.addWidget(QLabel("ЭКСПЕРИМЕНТОВ (N):"))
        self.sb_n = QSpinBox()
        self.sb_n.setRange(100, 100000)
        self.sb_n.setValue(5000)
        control_panel.addWidget(self.sb_n)

        self.btn_run = QPushButton("ЗАПУСТИТЬ АНАЛИЗ")
        self.btn_run.clicked.connect(self.process)
        control_panel.addWidget(self.btn_run)

        self.results_log = QTextEdit()
        self.results_log.setReadOnly(True)
        control_panel.addWidget(self.results_log)
        
        layout.addLayout(control_panel, 1)

        self.canvas = FigureCanvas(plt.Figure(figsize=(7, 5), facecolor='#121212'))
        layout.addWidget(self.canvas, 3)

    def process(self):
        lam = self.sb_lambda.value()
        t = self.sb_t.value()
        n = self.sb_n.value()
        expected_mu = lam * t

        data = PoissonEngine.simulate(lam, t, n)
        
        mean_emp = np.mean(data)
        var_emp = np.var(data)

        self.canvas.figure.clear()
        ax = self.canvas.figure.add_subplot(111)
        ax.set_facecolor('#1a1a1a')
        
        bins = np.arange(min(data), max(data) + 2) - 0.5
        ax.hist(data, bins=bins, density=True, color='#004400', edgecolor='#00ff41', alpha=0.7, label='Эмпирические данные')
        
        x = np.arange(0, max(data) + 1)
        y = poisson.pmf(x, expected_mu)
        ax.step(x, y, where='mid', color='#ff00ff', label=f'Теория (λT={expected_mu:.1f})', linewidth=2)

        ax.set_title("Распределение числа заявок", color='#00ff41')
        ax.tick_params(colors='#00ff41')
        ax.legend()
        self.canvas.draw()

        report = (
            f"--- РЕЗУЛЬТАТЫ ---\n"
            f"Теоретическое λT: {expected_mu:.3f}\n"
            f"Среднее (M): {mean_emp:.3f}\n"
            f"Дисперсия (D): {var_emp:.3f}\n"
            f"Разность |M-D|: {abs(mean_emp - var_emp):.4f}\n"
        )
        self.results_log.setText(report)

if __name__ == "__main__":
    app = QApplication(sys.argv)
    w = LabWindow()
    w.show()
    sys.exit(app.exec())