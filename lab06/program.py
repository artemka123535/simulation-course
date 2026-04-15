import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from scipy.stats import chi2, norm
from PyQt6.QtWidgets import (QApplication, QMainWindow, QWidget, QVBoxLayout, 
                             QHBoxLayout, QPushButton, QTableWidget, QTableWidgetItem, 
                             QComboBox, QLabel, QSpinBox, QFrame, QHeaderView, QTabWidget)

LIGHT_STYLE = """
    QMainWindow, QWidget { 
        background-color: #f5f5f7; 
        color: #1d1d1f; 
        font-family: 'Segoe UI', 'Helvetica Neue', sans-serif; 
    }
    QTabWidget::pane { border: 1px solid #d2d2d7; background: white; border-radius: 5px; }
    QTabBar::tab {
        background: #e5e5e7; padding: 10px 20px; border-top-left-radius: 5px; border-top-right-radius: 5px;
        margin-right: 2px;
    }
    QTabBar::tab:selected { background: white; border-bottom: 2px solid #0071e3; font-weight: bold; }
    
    QTableWidget { 
        background-color: white; border: 1px solid #d2d2d7; 
        gridline-color: #f2f2f2; selection-background-color: #e8f0fe; 
    }
    QHeaderView::section { 
        background-color: #fbfbfd; color: #86868b; 
        padding: 8px; border: 1px solid #d2d2d7; font-weight: bold;
    }
    
    QPushButton { 
        background-color: #0071e3; color: white; border: none; 
        border-radius: 8px; padding: 12px; font-weight: 600; 
    }
    QPushButton:hover { background-color: #0077ed; }
    QPushButton:pressed { background-color: #005bb7; }

    QSpinBox, QComboBox { 
        background-color: white; border: 1px solid #d2d2d7; 
        border-radius: 6px; padding: 6px; selection-background-color: #0071e3; 
    }
    QLabel { color: #1d1d1f; }
"""

class StatCard(QFrame):
    def __init__(self, title):
        super().__init__()
        self.setFrameShape(QFrame.Shape.StyledPanel)
        self.setStyleSheet("""
            QFrame { 
                background-color: white; border-radius: 12px; 
                border: 1px solid #d2d2d7; border-left: 5px solid #0071e3;
            }
        """)
        layout = QVBoxLayout(self)
        self.title_lbl = QLabel(title.upper())
        self.title_lbl.setStyleSheet("color: #86868b; font-size: 11px; border: none; font-weight: bold;")
        self.val_lbl = QLabel("0.0000")
        self.val_lbl.setStyleSheet("color: #1d1d1f; font-size: 20px; font-weight: 700; border: none; margin: 4px 0;")
        self.status_lbl = QLabel("—")
        self.status_lbl.setStyleSheet("font-size: 11px; border: none;")
        
        layout.addWidget(self.title_lbl)
        layout.addWidget(self.val_lbl)
        layout.addWidget(self.status_lbl)

    def update_val(self, value, error_pct=None, is_pass=None):
        self.val_lbl.setText(f"{value:.4f}")
        if error_pct is not None:
            self.status_lbl.setText(f"Погрешность: {error_pct:.2f}%")
            self.status_lbl.setStyleSheet(f"color: {'#248a3d' if error_pct < 5 else '#b38f00'}; border: none; font-weight: 500;")
        if is_pass is not None:
            status_text = "ПРОЙДЕН ✅" if is_pass else "ОТКЛОНЕН ❌"
            self.status_lbl.setText(status_text)
            self.status_lbl.setStyleSheet(f"color: {'#248a3d' if is_pass else '#d70015'}; border: none; font-weight: bold;")

class ModernLabApp(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Моделирование случайных величин")
        self.resize(1200, 800)
        self.setStyleSheet(LIGHT_STYLE)
        self.init_ui()

    def init_ui(self):
        central = QWidget()
        self.setCentralWidget(central)
        main_layout = QVBoxLayout(central)

        self.tabs = QTabWidget()
        main_layout.addWidget(self.tabs)

        self.tab_d = QWidget()
        self.init_discrete_tab()
        self.tabs.addTab(self.tab_d, "Дискретные величины")

        self.tab_n = QWidget()
        self.init_normal_tab()
        self.tabs.addTab(self.tab_n, "Нормальное распределение")

    def init_discrete_tab(self):
        layout = QHBoxLayout(self.tab_d)
        
        side_panel = QVBoxLayout()
        side_panel.setContentsMargins(15, 15, 15, 15)
        side_panel.setSpacing(10)
        
        side_panel.addWidget(QLabel("<b>Настройки выборки</b>"))
        
        self.n_input = QSpinBox()
        self.n_input.setRange(10, 1000000)
        self.n_input.setValue(1000)
        side_panel.addWidget(QLabel("Размер (N):"))
        side_panel.addWidget(self.n_input)

        self.presets = QComboBox()
        self.presets.addItems(["Пустая таблица", "Равномерное", "Бернулли"])
        self.presets.currentIndexChanged.connect(self.load_preset)
        side_panel.addWidget(QLabel("Выбрать пресет:"))
        side_panel.addWidget(self.presets)

        self.table = QTableWidget(5, 3)
        self.table.setHorizontalHeaderLabels(["X", "P теор.", "P эмп."])
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        side_panel.addWidget(self.table)

        self.btn_calc = QPushButton("Сгенерировать")
        self.btn_calc.clicked.connect(self.calculate_discrete)
        side_panel.addWidget(self.btn_calc)

        self.card_m = StatCard("Математическое ожидание")
        self.card_v = StatCard("Дисперсия")
        self.card_chi = StatCard("Критерий χ²")
        side_panel.addWidget(self.card_m)
        side_panel.addWidget(self.card_v)
        side_panel.addWidget(self.card_chi)
        side_panel.addStretch()

        self.fig_d, self.ax_d = plt.subplots(facecolor='white', tight_layout=True)
        self.canvas_d = FigureCanvas(self.fig_d)
        
        layout.addLayout(side_panel, 1)
        layout.addWidget(self.canvas_d, 2)

    def init_normal_tab(self):
        layout = QHBoxLayout(self.tab_n)
        
        side_panel = QVBoxLayout()
        side_panel.setContentsMargins(15, 15, 15, 15)
        
        side_panel.addWidget(QLabel("<b>Настройки выборки</b>"))
        self.n_norm = QSpinBox()
        self.n_norm.setRange(100, 1000000)
        self.n_norm.setValue(10000)
        side_panel.addWidget(QLabel("Объем выборки N:"))
        side_panel.addWidget(self.n_norm)
        
        btn = QPushButton("Генерировать")
        btn.clicked.connect(self.calculate_normal)
        side_panel.addWidget(btn)

        self.card_m_norm = StatCard("Математическое ожидание")
        self.card_v_norm = StatCard("Дисперсия")
        side_panel.addWidget(self.card_m_norm)
        side_panel.addWidget(self.card_v_norm)
        side_panel.addStretch()
        
        layout.addLayout(side_panel, 1)

        self.fig_n, self.ax_n = plt.subplots(facecolor='white', tight_layout=True)
        self.canvas_n = FigureCanvas(self.fig_n)
        layout.addWidget(self.canvas_n, 2)

    def load_preset(self):
        idx = self.presets.currentIndex()
        self.table.setRowCount(0)
        data = []
        if idx == 1: data = [(i, 0.2) for i in range(1, 6)]
        elif idx == 2: data = [(0, 0.45), (1, 0.55)]
        
        for x, p in data:
            row = self.table.rowCount()
            self.table.insertRow(row)
            self.table.setItem(row, 0, QTableWidgetItem(str(x)))
            self.table.setItem(row, 1, QTableWidgetItem(str(p)))

    def calculate_discrete(self):
        try:
            n = self.n_input.value()
            rows = self.table.rowCount()
            x = np.array([float(self.table.item(i, 0).text()) for i in range(rows)])
            p = np.array([float(self.table.item(i, 1).text()) for i in range(rows)])
            p /= p.sum()

            samples = np.random.choice(x, size=n, p=p)
            
            unique, counts = np.unique(samples, return_counts=True)
            obs_map = dict(zip(unique, counts))
            
            m_th, d_th = np.sum(x * p), np.sum((x**2) * p) - np.sum(x * p)**2
            m_emp, d_emp = np.mean(samples), np.var(samples)
            
            err_m = abs(m_th - m_emp) / (abs(m_th) if m_th != 0 else 1) * 100
            err_d = abs(d_th - d_emp) / (d_th if d_th != 0 else 1) * 100

            chi_stat = 0
            for i in range(rows):
                obs = obs_map.get(x[i], 0)
                exp = n * p[i]
                chi_stat += ((obs - exp)**2) / exp
                self.table.setItem(i, 2, QTableWidgetItem(f"{obs/n:.4f}"))

            self.card_m.update_val(m_emp, error_pct=err_m)
            self.card_v.update_val(d_emp, error_pct=err_d)
            self.card_chi.update_val(chi_stat, is_pass=(chi_stat < chi2.ppf(0.95, rows-1)))

            self.ax_d.clear()
            self.ax_d.bar(x, [obs_map.get(val, 0)/n for val in x], color='#0071e3', alpha=0.15, label='Эмпирические')
            self.ax_d.stem(x, p, linefmt='#0071e3', markerfmt='o', basefmt=" ", label='Теоретические')
            self.ax_d.set_title(f"Сравнение распределений (N={n})", pad=15, fontsize=12, fontweight='bold')
            self.ax_d.grid(axis='y', linestyle='--', alpha=0.3)
            self.ax_d.legend()
            self.canvas_d.draw()
        except: pass

    def calculate_normal(self):
        n = self.n_norm.value()
        u1, u2 = np.random.rand(n), np.random.rand(n)
        samples = np.sqrt(-2 * np.log(u1)) * np.cos(2 * np.pi * u2)

        m_emp = np.mean(samples)
        v_emp = np.var(samples)
        
        err_m = abs(m_emp - 0) * 100 
        err_v = abs(v_emp - 1) / 1 * 100

        self.card_m_norm.update_val(m_emp, error_pct=err_m)
        self.card_v_norm.update_val(v_emp, error_pct=err_v)

        self.ax_n.clear()
        self.ax_n.hist(samples, bins=60, density=True, color='#0071e3', alpha=0.15, edgecolor='#0071e3')
        
        x_axis = np.linspace(-4, 4, 200)
        self.ax_n.plot(x_axis, norm.pdf(x_axis), color='#d70015', lw=2, label='Теоретическая плотность')
        self.ax_n.set_title(f"Анализ нормального распределения (N={n})", fontsize=12, fontweight='bold')
        self.ax_n.grid(axis='both', linestyle='--', alpha=0.2)
        self.ax_n.legend()
        self.canvas_n.draw()

if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = ModernLabApp()
    window.show()
    sys.exit(app.exec())