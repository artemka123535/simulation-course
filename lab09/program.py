import sys
import numpy as np
import matplotlib.pyplot as plt
from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, 
    QHBoxLayout, QLabel, QPushButton, QDoubleSpinBox, 
    QSpinBox, QTextEdit
)
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas

STYLE_SHEET = """
QMainWindow { background-color: #121212; }
QLabel { color: #00ff41; font-family: 'Courier New'; font-weight: bold; }
QGroupBox { border: 1px solid #00ff41; color: #00ff41; margin-top: 10px; }
QPushButton { 
    background-color: #003300; color: #00ff41; 
    border: 1px solid #00ff41; padding: 10px; border-radius: 4px; font-weight: bold;
}
QPushButton:hover { background-color: #004400; }
QDoubleSpinBox, QSpinBox { 
    background-color: #1a1a1a; color: #00ff41; 
    border: 1px solid #004400; padding: 2px;
}
QTextEdit { background-color: #0a0a0a; color: #00ff41; border: 1px solid #003300; font-size: 12px; font-family: 'Courier New'; }
"""

class MM1LossEngine:
    @staticmethod
    def simulate(lam, mu, n_requests):
        accepted = 0
        rejected = 0
        
        inter_arrivals = np.random.exponential(1.0 / lam, n_requests)
        service_times = np.random.exponential(1.0 / mu, n_requests)
        
        arrival_times = np.cumsum(inter_arrivals)
        
        server_free_time = 0.0
        
        for i in range(n_requests):
            arrival = arrival_times[i]
            if arrival >= server_free_time:
                accepted += 1
                server_free_time = arrival + service_times[i]
            else:
                rejected += 1
                
        p0_emp = accepted / n_requests
        p1_emp = rejected / n_requests
        
        return p0_emp, p1_emp

class QueuingLabWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("СМО M/M/1/")
        self.resize(1100, 700)
        self.setStyleSheet(STYLE_SHEET)

        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        layout = QHBoxLayout(central_widget)

        control_panel = QVBoxLayout()
        
        control_panel.addWidget(QLabel("ПОСТУПЛЕНИЕ λ (заявок/ед.вр):"))
        self.sb_lambda = QDoubleSpinBox()
        self.sb_lambda.setRange(0.1, 100.0)
        self.sb_lambda.setValue(5.0)
        self.sb_lambda.setSingleStep(0.5)
        control_panel.addWidget(self.sb_lambda)

        control_panel.addWidget(QLabel("ОБСЛУЖИВАНИЕ μ (заявок/ед.вр):"))
        self.sb_mu = QDoubleSpinBox()
        self.sb_mu.setRange(0.1, 100.0)
        self.sb_mu.setValue(6.0)
        self.sb_mu.setSingleStep(0.5)
        control_panel.addWidget(self.sb_mu)

        control_panel.addWidget(QLabel("КОЛИЧЕСТВО ЗАЯВОК (N):"))
        self.sb_n = QSpinBox()
        self.sb_n.setRange(1, 1000000)
        self.sb_n.setValue(10000)
        self.sb_n.setSingleStep(1000)
        control_panel.addWidget(self.sb_n)

        self.btn_run = QPushButton("ЗАПУСТИТЬ СИМУЛЯЦИЮ")
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
        mu = self.sb_mu.value()
        n = self.sb_n.value()

        p0_emp, p1_emp = MM1LossEngine.simulate(lam, mu, n)
        
        rho = lam / mu
        p0_theor = 1.0 / (1.0 + rho)
        p1_theor = rho / (1.0 + rho)

        self.canvas.figure.clear()
        ax = self.canvas.figure.add_subplot(111)
        ax.set_facecolor('#1a1a1a')
        
        labels = ['Свободен (P0)', 'Занят (P1)']
        emp_vals = [p0_emp, p1_emp]
        theor_vals = [p0_theor, p1_theor]
        
        x = np.arange(len(labels))
        width = 0.35
        
        ax.bar(x - width/2, emp_vals, width, label='Эмпирика', color='#00ff41', edgecolor='#004400')
        ax.bar(x + width/2, theor_vals, width, label='Теория', color='#ff00ff', edgecolor='#550055', alpha=0.7)

        ax.set_title(f"Сравнение вероятностей состояний СМО (ρ = {rho:.2f})", color='#00ff41')
        ax.set_xticks(x)
        ax.set_xticklabels(labels, color='#00ff41')
        ax.tick_params(axis='y', colors='#00ff41')
        ax.legend()
        ax.grid(color='#003300', linestyle='--', alpha=0.5)
        
        self.canvas.draw()
        report = (
            f"--- СТАТИСТИКА СМО ---\n\n"
            f"Нагрузка на систему (ρ = λ/μ): {rho:.3f}\n\n"
            
            f"--- ВЕРОЯТНОСТЬ P0 ---\n"
            f"Теория:   {p0_theor:.4f}\n"
            f"Эмпирика: {p0_emp:.4f}\n"
            f"Разница:  {abs(p0_theor - p0_emp):.4f}\n\n"
            
            f"--- ВЕРОЯТНОСТЬ P1 ---\n"
            f"Теория:   {p1_theor:.4f}\n"
            f"Эмпирика: {p1_emp:.4f}\n"
            f"Разница:  {abs(p1_theor - p1_emp):.4f}\n\n"
            
            f"Абсолютная пропускная способность:\n"
            f"A = λ * P0 = {(lam * p0_emp):.3f} заявок/ед.вр."
        )
        self.results_log.setText(report)

if __name__ == "__main__":
    app = QApplication(sys.argv)
    w = QueuingLabWindow()
    w.show()
    sys.exit(app.exec())