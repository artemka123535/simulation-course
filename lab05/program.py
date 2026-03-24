import tkinter as tk
from tkinter import ttk, scrolledtext
from collections import Counter
import random

MAGIC_ANSWERS = [
    "Определённо да",
    "Скорее всего",
    "Хорошие шансы",
    "Возможно, да",
    "Не могу сказать сейчас",
    "Спроси позже",
    "Шансы невелики",
    "Скорее нет",
    "Мой ответ — нет",
    "Очень сомнительно"
]

def get_yes_no(p=0.5):
    return "ДА" if random.random() < p else "НЕТ"

def get_magic_8_ball():
    m = len(MAGIC_ANSWERS)
    p_i = 1.0 / m  

    alpha = random.random()

    cumulative_p = 0.0

    for k in range(m):
        cumulative_p += p_i

        if alpha < cumulative_p:
            return MAGIC_ANSWERS[k]
        
    return MAGIC_ANSWERS[-1]

def run_simulation(mode, n_trials=10000):
    if mode == "YESNO":
        results = [get_yes_no() for _ in range(n_trials)]
        keys = ["ДА", "НЕТ"]
        theoretical_p = 0.5
    else:  # MAGIC
        results = [get_magic_8_ball() for _ in range(n_trials)]
        keys = MAGIC_ANSWERS
        theoretical_p = 1.0 / len(MAGIC_ANSWERS)

    counts = Counter(results)
    stats = []
    for key in keys:
        n_k = counts.get(key, 0)
        p_hat = n_k / n_trials
        stats.append({
            'answer': key,
            'count': n_k,
            'p_empirical': p_hat,
            'p_theory': theoretical_p,
            'diff': abs(p_hat - theoretical_p)
        })
    return stats, n_trials

class PredictionApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Предсказания")
        self.root.geometry("600x500")
        self.root.resizable(True, True)

        self.notebook = ttk.Notebook(root)
        self.notebook.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

        self.frame_yesno = ttk.Frame(self.notebook)
        self.notebook.add(self.frame_yesno, text="Да / Нет")
        self.setup_yesno_tab()

        self.frame_magic = ttk.Frame(self.notebook)
        self.notebook.add(self.frame_magic, text="Шар предсказаний")
        self.setup_magic_tab()
        
        self.frame_sim = ttk.Frame(self.notebook)
        self.notebook.add(self.frame_sim, text="Статистика")
        self.setup_sim_tab()

    def setup_yesno_tab(self):
        tk.Label(self.frame_yesno, text="Задайте вопрос:", font=("Arial", 12)).pack(pady=10)
        self.question_yesno = tk.Entry(self.frame_yesno, width=50)
        self.question_yesno.pack(pady=5)

        btn_ask = tk.Button(self.frame_yesno, text="Спросить", command=self.ask_yesno,
                            font=("Arial", 12))
        btn_ask.pack(pady=10)

        self.label_answer_yesno = tk.Label(self.frame_yesno, text="", font=("Arial", 24, "bold"), fg="blue")
        self.label_answer_yesno.pack(pady=20)

        btn_clear = tk.Button(self.frame_yesno, text="Очистить вопрос", command=self.clear_question_yesno)
        btn_clear.pack()

    def setup_magic_tab(self):
        tk.Label(self.frame_magic, text="Задайте вопрос магическому шару:", font=("Arial", 12)).pack(pady=10)
        self.question_magic = tk.Entry(self.frame_magic, width=50)
        self.question_magic.pack(pady=5)

        btn_ask = tk.Button(self.frame_magic, text="Трясти шар", command=self.ask_magic,
                            font=("Arial", 12))
        btn_ask.pack(pady=10)

        self.label_answer_magic = tk.Label(self.frame_magic, text="", font=("Arial", 20, "bold"), fg="purple")
        self.label_answer_magic.pack(pady=20)

        btn_clear = tk.Button(self.frame_magic, text="Очистить вопрос", command=self.clear_question_magic)
        btn_clear.pack()

    def setup_sim_tab(self):
        tk.Label(self.frame_sim, text="Режим симуляции:", font=("Arial", 12)).pack(pady=5)
        self.sim_mode = tk.StringVar(value="YESNO")
        mode_frame = ttk.Frame(self.frame_sim)
        mode_frame.pack()
        ttk.Radiobutton(mode_frame, text="Да / Нет", variable=self.sim_mode, value="YESNO").pack(side=tk.LEFT, padx=10)
        ttk.Radiobutton(mode_frame, text="Шар предсказаний", variable=self.sim_mode, value="MAGIC").pack(side=tk.LEFT, padx=10)

        tk.Label(self.frame_sim, text="Количество испытаний:", font=("Arial", 12)).pack(pady=5)
        self.trials_entry = tk.Entry(self.frame_sim, width=15)
        self.trials_entry.insert(0, "10000")
        self.trials_entry.pack()

        btn_run = tk.Button(self.frame_sim, text="Запустить симуляцию", command=self.run_simulation_gui,
                            font=("Arial", 12))
        btn_run.pack(pady=10)

        self.sim_output = scrolledtext.ScrolledText(self.frame_sim, width=70, height=20, wrap=tk.WORD)
        self.sim_output.pack(padx=5, pady=5, fill=tk.BOTH, expand=True)

    def ask_yesno(self):
        answer = get_yes_no()
        self.label_answer_yesno.config(text=answer)

    def ask_magic(self):
        answer = get_magic_8_ball()
        self.label_answer_magic.config(text=answer)

    def clear_question_yesno(self):
        self.question_yesno.delete(0, tk.END)

    def clear_question_magic(self):
        self.question_magic.delete(0, tk.END)

    def run_simulation_gui(self):
        try:
            n_trials = int(self.trials_entry.get())
            if n_trials <= 0:
                raise ValueError
        except ValueError:
            self.sim_output.insert(tk.END, "Ошибка: введите положительное целое число.\n")
            return

        mode = self.sim_mode.get()
        stats, total = run_simulation(mode, n_trials)

        self.sim_output.delete(1.0, tk.END)
        self.sim_output.insert(tk.END, f"Результаты симуляции (n = {total})\n")
        self.sim_output.insert(tk.END, f"Режим: {'Да/Нет' if mode == 'YESNO' else 'Шар предсказаний'}\n")
        self.sim_output.insert(tk.END, "-" * 60 + "\n")
        self.sim_output.insert(tk.END, f"{'Ответ':<20} {'Частота':<10} {'Эмп. вер.':<12} {'Теор. вер.':<12} {'Откл.':<10}\n")
        self.sim_output.insert(tk.END, "-" * 60 + "\n")
        for row in stats:
            ans = row['answer']
            if len(ans) > 18:
                ans = ans[:15] + "..."
            self.sim_output.insert(tk.END,
                f"{ans:<20} {row['count']:<10} {row['p_empirical']:<12.5f} {row['p_theory']:<12.5f} {row['diff']:<10.5f}\n")
        self.sim_output.insert(tk.END, "-" * 60 + "\n")
        self.sim_output.see(tk.END)

if __name__ == "__main__":
    root = tk.Tk()
    app = PredictionApp(root)
    root.mainloop()