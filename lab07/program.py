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

class WeatherContinuousMarkov(ctk.CTk):
    def __init__(self):
        super().__init__()
        self.title("Марковская модель погоды")
        self.geometry("1400x900")

        self.running = False
        self.total_time = 0.0
        self.history_time = [0.0]
        self.history_states = [1]
        self.time_spent = {1: 0.0, 2: 0.0, 3: 0.0}
        
        self.states = {1: "Ясно ☀️", 2: "Облачно ☁️", 3: "Пасмурно 🌧️"}
        self.current_state = 1
        self.csv_file = "lab07/weather_analysis.csv"
        
        self.reset_csv()
        self.setup_ui()

    def reset_csv(self):
        with open(self.csv_file, 'w', newline='', encoding='utf-16') as f:
            writer = csv.writer(f)
            writer.writerow(["Total_Time", "State_ID", "Duration"])

    def setup_ui(self):
        self.grid_columnconfigure(1, weight=1)
        self.grid_rowconfigure(0, weight=1)

        self.sidebar = ctk.CTkFrame(self, width=300)
        self.sidebar.grid(row=0, column=0, sticky="nsew", padx=15, pady=15)
        
        ctk.CTkLabel(self.sidebar, text="Матрица интенсивностей", font=("Arial", 18, "bold")).pack(pady=15)

        self.matrix_inputs = []
        labels = ["Ясно", "Облачно", "Пасмурно"]
        for i in range(3):
            ctk.CTkLabel(self.sidebar, text=f"Интенсивности из {labels[i]}:", font=("Arial", 13)).pack(anchor="w", padx=25, pady=(15, 0))
            f = ctk.CTkFrame(self.sidebar, fg_color="transparent")
            f.pack(pady=5)
            row_entries = []
            for j in range(3):
                e = ctk.CTkEntry(f, width=70)
                default_val = "0.5" if i != j else "0.0"
                if i == 0 and j == 1: default_val = "0.6"
                if i == 2 and j == 1: default_val = "0.8"
                
                e.insert(0, default_val)
                if i == j: 
                    e.configure(state="disabled", fg_color="#2c3e50") 
                e.pack(side="left", padx=3)
                row_entries.append(e)
            self.matrix_inputs.append(row_entries)

        ctk.CTkLabel(self.sidebar, text="Визуальный масштаб времени:").pack(pady=(30, 0))
        self.speed_slider = ctk.CTkSlider(self.sidebar, from_=0.1, to=5.0)
        self.speed_slider.set(1.5)
        self.speed_slider.pack(pady=10, padx=25)

        self.start_btn = ctk.CTkButton(self.sidebar, text="Запустить симуляцию", command=self.toggle_simulation, 
                                       fg_color="#2ecc71", font=("Arial", 14, "bold"), height=45)
        self.start_btn.pack(pady=30, padx=25)
        
        self.reset_btn = ctk.CTkButton(self.sidebar, text="Сбросить всё", command=self.reset_sim, 
                                       fg_color="#e74c3c", font=("Arial", 14, "bold"), height=45)
        self.reset_btn.pack(pady=5, padx=25)

        self.main_content = ctk.CTkFrame(self, fg_color="transparent")
        self.main_content.grid(row=0, column=1, sticky="nsew", padx=15, pady=15)

        self.info_card = ctk.CTkFrame(self.main_content, height=120)
        self.info_card.pack(fill="x", pady=(0, 15))
        
        self.time_label = ctk.CTkLabel(self.info_card, text="Модельное время: 0.00", font=("Arial", 22))
        self.time_label.pack(side="left", padx=40, pady=25)
        
        self.state_label = ctk.CTkLabel(self.info_card, text="ОЖИДАНИЕ", font=("Arial", 36, "bold"), width=400)
        self.state_label.pack(side="right", padx=40, pady=25)

        self.fig, (self.ax_line, self.ax_bar) = plt.subplots(2, 1, figsize=(10, 8), facecolor='#1a1a1a')
        self.fig.subplots_adjust(left=0.15, bottom=0.1, right=0.95, top=0.9, hspace=0.4)
        
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.main_content)
        self.canvas.get_tk_widget().pack(fill="both", expand=True)

    def get_q_matrix(self):
        try:
            Q = np.zeros((3, 3))
            for i in range(3):
                row_sum = 0
                for j in range(3):
                    if i != j:
                        val = float(self.matrix_inputs[i][j].get())
                        Q[i, j] = val
                        row_sum += val
                Q[i, i] = -row_sum 
            return Q
        except ValueError:
            return None

    def calculate_stationary(self, Q):
        try:
            A = Q.T
            print(A)
            A[-1] = np.ones(3)
            b = np.zeros(3)
            b[-1] = 1
            return np.linalg.solve(A, b)
        except Exception:
            return np.array([0.33, 0.33, 0.34])

    def toggle_simulation(self):
        if not self.running:
            if self.get_q_matrix() is None:
                tk.messagebox.showerror("Ошибка", "Введите корректные числа в матрицу")
                return
            self.running = True
            self.start_btn.configure(text="Остановить", fg_color="#f39c12")
            threading.Thread(target=self.simulation_loop, daemon=True).start()
        else:
            self.running = False
            self.start_btn.configure(text="Запустить симуляцию", fg_color="#2ecc71")

    def simulation_loop(self):
        while self.running:
            Q = self.get_q_matrix()
            theo_probs = self.calculate_stationary(Q)
            
            curr_i = self.current_state - 1
            lambda_out = -Q[curr_i, curr_i]
            
            if lambda_out > 0:
                dt = np.random.exponential(1.0 / lambda_out)
            else:
                dt = 1.0
            
            visual_speed = self.speed_slider.get()
            time.sleep(dt / visual_speed) 
            
            if not self.running: break

            self.total_time += dt
            self.time_spent[self.current_state] += dt
            self.history_states.append(self.current_state)
            self.history_time.append(self.total_time)

            probabilities = []
            other_states = []
            for j in range(3):
                if curr_i != j:
                    probabilities.append(Q[curr_i, j] / lambda_out)
                    other_states.append(j + 1)
            
            probabilities = np.array(probabilities)
            probabilities /= probabilities.sum()

            self.update_ui_state(theo_probs)
            
            with open(self.csv_file, 'a', newline='', encoding='utf-16') as f:
                writer = csv.writer(f)
                writer.writerow([round(self.total_time, 3), self.current_state, round(dt, 3)])

            self.current_state = np.random.choice(other_states, p=probabilities)

    def update_ui_state(self, theo_probs):
        self.time_label.configure(text=f"Модельное время: {self.total_time:.2f}")
        self.state_label.configure(text=self.states[self.current_state])
        
        self.ax_line.clear()
        self.ax_line.set_facecolor('#1a1a1a')

        self.ax_line.step(self.history_time[-40:], self.history_states[-40:], where='post', color='#3498db', linewidth=2)
        self.ax_line.set_yticks([1, 2, 3])
        self.ax_line.set_yticklabels(["Ясно", "Облачно", "Пасмурно"], color="white", fontsize=11)
        self.ax_line.tick_params(colors='white')
        self.ax_line.set_title("Траектория процесса", color="white")

        self.ax_bar.clear()
        self.ax_bar.set_facecolor('#1a1a1a')
        labels = ['Ясно', 'Облачно', 'Пасмурно']
        emp_probs = [self.time_spent[i]/self.total_time if self.total_time > 0 else 0 for i in [1, 2, 3]]
        
        x = np.arange(len(labels))
        width = 0.35
        rects1 = self.ax_bar.bar(x - width/2, emp_probs, width, label='Эмп.', color='#3498db')
        rects2 = self.ax_bar.bar(x + width/2, theo_probs, width, label='Теор.', color='#f1c40f')

        self.ax_bar.set_xticks(x)
        self.ax_bar.set_xticklabels(labels, color='white')
        self.ax_bar.tick_params(colors='white')
        self.ax_bar.legend(facecolor='#1a1a1a', labelcolor='white')
        self.ax_bar.set_ylim(0, 1.1)
        self.ax_bar.bar_label(rects1, padding=3, color='#3498db', fmt='%.3f')
        self.ax_bar.bar_label(rects2, padding=3, color='#f1c40f', fmt='%.3f')

        self.canvas.draw()

    def reset_sim(self):
        self.running = False
        self.total_time = 0.0
        self.current_state = 1
        self.time_spent = {1: 0.0, 2: 0.0, 3: 0.0}
        self.history_time = [0.0]
        self.history_states = [1]
        self.reset_csv()
        
        self.time_label.configure(text="Модельное время: 0")
        self.state_label.configure(text="СБРОШЕНО")
        self.start_btn.configure(text="Запустить симуляцию", fg_color="#2ecc71")
        
        self.ax_line.clear()
        self.ax_bar.clear()
        self.canvas.draw()

if __name__ == "__main__":
    app = WeatherContinuousMarkov()
    app.mainloop()