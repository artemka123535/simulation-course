import tkinter as tk
import customtkinter as ctk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import csv
import threading
import time

ctk.set_appearance_mode("Dark")
ctk.set_default_color_theme("blue")

class WeatherMarkovApp(ctk.CTk):
    def __init__(self):
        super().__init__()
        self.title("Марковская модель погоды")
        self.geometry("1400x900")

        self.running = False
        self.day = 0
        self.history = []
        self.states = {1: "Ясно ☀️", 2: "Облачно ☁️", 3: "Пасмурно 🌧️"}
        self.current_state = 1
        self.csv_file = "lab07/weather_analysis.csv"
        
        self.setup_ui()

    def setup_ui(self):
        self.grid_columnconfigure(1, weight=1)
        self.grid_rowconfigure(0, weight=1)

        self.sidebar = ctk.CTkFrame(self, width=250)
        self.sidebar.grid(row=0, column=0, sticky="nsew", padx=10, pady=10)
        
        ctk.CTkLabel(self.sidebar, text="Матрица переходов", font=("Arial", 16, "bold")).pack(pady=10)

        self.matrix_inputs = []
        for i in range(3):
            f = ctk.CTkFrame(self.sidebar, fg_color="transparent")
            f.pack(pady=2)
            row_entries = []
            for j in range(3):
                e = ctk.CTkEntry(f, width=50)
                e.insert(0, str([0.7, 0.2, 0.1][j] if i==0 else [0.3, 0.4, 0.3][j] if i==1 else [0.2, 0.3, 0.5][j]))
                e.pack(side="left", padx=2)
                row_entries.append(e)
            self.matrix_inputs.append(row_entries)

        ctk.CTkLabel(self.sidebar, text="Задержка:").pack(pady=(20, 0))
        self.speed_slider = ctk.CTkSlider(self.sidebar, from_=0.0, to=0.5)
        self.speed_slider.set(0.05)
        self.speed_slider.pack(pady=10, padx=20)

        self.start_btn = ctk.CTkButton(self.sidebar, text="Запустить", command=self.toggle_simulation, fg_color="#2ecc71")
        self.start_btn.pack(pady=20, padx=20)
        
        self.reset_btn = ctk.CTkButton(self.sidebar, text="Сброс данных", command=self.reset_sim, fg_color="#e74c3c")
        self.reset_btn.pack(pady=5, padx=20)

        self.main_content = ctk.CTkFrame(self, fg_color="transparent")
        self.main_content.grid(row=0, column=1, sticky="nsew", padx=15, pady=15)

        self.info_card = ctk.CTkFrame(self.main_content)
        self.info_card.pack(fill="x", pady=(0, 10))
        self.day_text = ctk.CTkLabel(self.info_card, text="День: 0", font=("Arial", 20))
        self.day_text.pack(side="left", padx=30, pady=20)
        self.state_text = ctk.CTkLabel(self.info_card, text="ОЖИДАНИЕ", font=("Arial", 32, "bold"))
        self.state_text.pack(side="right", padx=30, pady=20)

        self.fig, (self.ax_line, self.ax_bar) = plt.subplots(2, 1, figsize=(8, 8), facecolor='#1a1a1a')
        self.fig.tight_layout(pad=5.0)
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.main_content)
        self.canvas.get_tk_widget().pack(fill="both", expand=True)

    def calculate_stationary(self, P):
        try:
            n = P.shape[0]
            A = np.transpose(P) - np.eye(n)
            A[-1] = np.ones(n)
            B = np.zeros(n)
            B[-1] = 1
            return np.linalg.solve(A, B)
        except: return np.array([0.33, 0.33, 0.34])

    def get_matrix(self):
        try: return np.array([[float(e.get()) for e in row] for row in self.matrix_inputs])
        except: return None

    def toggle_simulation(self):
        if not self.running:
            m = self.get_matrix()
            if m is None or not np.allclose(m.sum(axis=1), 1):
                tk.messagebox.showerror("Ошибка", "Сумма строки в матрице должна быть 1.0")
                return
            self.running = True
            self.start_btn.configure(text="Стоп", fg_color="#f39c12")
            threading.Thread(target=self.simulation_loop, daemon=True).start()
        else:
            self.running = False
            self.start_btn.configure(text="Запустить", fg_color="#2ecc71")

    def simulation_loop(self):
        matrix = self.get_matrix()
        theo_probs = self.calculate_stationary(matrix)
        
        while self.running:
            self.day += 1
            self.history.append(self.current_state)
            
            self.day_text.configure(text=f"День: {self.day}")
            self.state_text.configure(text=self.states[self.current_state])
            
            self.update_charts(theo_probs)

            with open(self.csv_file, 'a', newline='', encoding='utf-16') as f:
                writer = csv.writer(f)
                writer.writerow([self.day, self.current_state, self.states[self.current_state]])

            self.current_state = np.random.choice([1, 2, 3], p=matrix[self.current_state - 1])
            time.sleep(self.speed_slider.get())

    def update_charts(self, theo_probs):
        self.ax_line.clear()
        self.ax_line.set_facecolor('#1a1a1a')
        self.ax_line.plot(self.history[-100:], color='#3498db', marker='o', markersize=3, label="История (посл. 100 дней)")
        self.ax_line.set_yticks([1, 2, 3])
        self.ax_line.set_yticklabels(["Ясно", "Облачно", "Пасмурно"], color="white")
        self.ax_line.tick_params(colors='white')
        self.ax_line.grid(True, alpha=0.1)

        self.ax_bar.clear()
        self.ax_bar.set_facecolor('#1a1a1a')
        
        labels = ['Ясно', 'Облачно', 'Пасмурно']
        x = np.arange(len(labels))
        width = 0.35
        
        counts = [self.history.count(i)/len(self.history) for i in [1, 2, 3]]
        
        rects1 = self.ax_bar.bar(x - width/2, counts, width, label='Эмпирическое', color='#3498db')
        rects2 = self.ax_bar.bar(x + width/2, theo_probs, width, label='Теоретическое', color='#f1c40f')

        self.ax_bar.set_ylabel('Вероятность', color='white')
        self.ax_bar.set_title(f'Сравнение распределений (n={self.day})', color='white', pad=10)
        self.ax_bar.set_xticks(x)
        self.ax_bar.set_xticklabels(labels, color='white')
        self.ax_bar.legend()
        self.ax_bar.tick_params(colors='white')
        self.ax_bar.set_ylim(0, 1.1)

        self.ax_bar.bar_label(rects1, padding=3, color='#3498db', fmt='%.3f')
        self.ax_bar.bar_label(rects2, padding=3, color='#f1c40f', fmt='%.3f')

        self.canvas.draw()

    def reset_sim(self):
        self.running = False
        self.day = 0
        self.history = []
        self.current_state = 1
        self.day_text.configure(text="День: 0")
        self.state_text.configure(text="СБРОШЕНО")
        
        with open(self.csv_file, 'w', newline='', encoding='utf-16') as f:
            writer = csv.writer(f)
            writer.writerow(["Day", "State", "Name"])
            
        self.ax_line.clear()
        self.ax_bar.clear()
        self.canvas.draw()

if __name__ == "__main__":
    app = WeatherMarkovApp()
    app.mainloop()