import tkinter as tk
from tkinter import ttk
import numpy as np
import random

EMPTY = 0
ASH = 1
TREE_YOUNG = 2
TREE_MEDIUM = 3
TREE_OLD = 4
BURNING = 5
WATER = 6

COLORS = {
    EMPTY: "#744000",       
    ASH: "#4a4a4a",         
    TREE_YOUNG: "#00cc00",  
    TREE_MEDIUM: "#008200", 
    TREE_OLD: "#004F00",    
    BURNING: "#ff4500",     
    WATER: "#1e90ff"       
}

class ForestFireWindModel:
    def __init__(self, root):
        self.root = root
        self.root.title("Клеточный автомат")
        
        self.size = 50
        self.grid = np.zeros((self.size, self.size), dtype=int)
        
        self.prob_age = 0.005
        self.prob_ash_clear = 0.05
        self.wind_dir = (0, 0)
        
        self.running = False
        self.setup_ui()
        
        self.canvas.bind("<Button-1>", self.add_water)
        self.canvas.bind("<B1-Motion>", self.add_water)

    def setup_ui(self):
        controls = ttk.Frame(self.root, padding="10")
        controls.pack(side=tk.RIGHT, fill=tk.Y)

        self.btn_start = ttk.Button(controls, text="Старт", command=self.toggle)
        self.btn_start.pack(fill=tk.X, pady=5)
        
        ttk.Button(controls, text="Очистить поле", command=self.clear_all).pack(fill=tk.X)
        
        ttk.Label(controls, text="\nВероятность роста (p):").pack()
        self.grow_slider = ttk.Scale(controls, from_=0, to=0.05, value=0.01)
        self.grow_slider.pack(fill=tk.X)

        ttk.Label(controls, text="Вероятность молнии (f):").pack()
        self.f_slider = ttk.Scale(controls, from_=0, to=0.005, value=0.0005)
        self.f_slider.pack(fill=tk.X)

        ttk.Label(controls, text="\nНаправление ветра:", font=('Arial', 10, 'bold')).pack(pady=5)
        self.wind_var = tk.StringVar(value="None")
        directions = {
            "Нет": (0, 0),
            "Север ↑": (-1, 0),
            "Юг ↓": (1, 0),
            "Восток →": (0, 1),
            "Запад ←": (0, -1)
        }
        for text, vec in directions.items():
            ttk.Radiobutton(controls, text=text, variable=self.wind_var, 
                            value=text, command=lambda v=vec: self.set_wind(v)).pack(anchor=tk.W)

        self.c_size = 600
        self.cell_w = self.c_size // self.size
        self.canvas = tk.Canvas(self.root, width=self.c_size, height=self.c_size, bg="#744000")
        self.canvas.pack(side=tk.LEFT)
        self.draw_grid()

    def set_wind(self, vec):
        self.wind_dir = vec

    def add_water(self, event):
        x, y = event.x // self.cell_w, event.y // self.cell_w
        if 0 <= x < self.size and 0 <= y < self.size:
            self.grid[y, x] = WATER
            self.draw_cell(x, y)

    def toggle(self):
        self.running = not self.running
        self.btn_start.config(text="Стоп" if self.running else "Старт")
        if self.running: self.loop()

    def clear_all(self):
        self.grid.fill(EMPTY)
        self.draw_grid()

    def draw_cell(self, x, y):
        self.canvas.create_rectangle(x*self.cell_w, y*self.cell_w, 
                                   (x+1)*self.cell_w, (y+1)*self.cell_w, 
                                   fill=COLORS[self.grid[y, x]], outline="")

    def draw_grid(self):
        self.canvas.delete("all")
        for y in range(self.size):
            for x in range(self.size):
                if self.grid[y, x] != EMPTY:
                    self.draw_cell(x, y)

    def loop(self):
        if not self.running: return
        
        new_grid = self.grid.copy()
        p_grow = self.grow_slider.get()
        f_lightning = self.f_slider.get()
        wy, wx = self.wind_dir

        for y in range(self.size):
            for x in range(self.size):
                state = self.grid[y, x]
                
                if state == ASH:
                    if random.random() < self.prob_ash_clear:
                        new_grid[y, x] = EMPTY
                
                elif state == EMPTY:
                    if random.random() < p_grow:
                        new_grid[y, x] = TREE_YOUNG
                
                elif state in [TREE_YOUNG, TREE_MEDIUM, TREE_OLD]:
                    if state == TREE_YOUNG and random.random() < self.prob_age:
                        new_grid[y, x] = TREE_MEDIUM
                    elif state == TREE_MEDIUM and random.random() < self.prob_age:
                        new_grid[y, x] = TREE_OLD
                    
                    if random.random() < f_lightning:
                        new_grid[y, x] = BURNING
                    else:
                        is_ignited = False
                        for dy in [-1, 0, 1]:
                            for dx in [-1, 0, 1]:
                                if dy == 0 and dx == 0: continue
                                ny, nx = y + dy, x + dx
                                
                                if 0 <= ny < self.size and 0 <= nx < self.size:
                                    if self.grid[ny, nx] == BURNING:
                                        base_chance = 0.2 if state == TREE_YOUNG else (0.5 if state == TREE_MEDIUM else 0.8)
                                        
                                        if (wy == -dy and wy != 0) or (wx == -dx and wx != 0):
                                            base_chance += 0.5
                                        elif (wy == dy and wy != 0) or (wx == dx and wx != 0):
                                            base_chance -= 0.15
                                            
                                        if random.random() < base_chance:
                                            is_ignited = True
                                            break
                            if is_ignited: break
                        
                        if is_ignited:
                            new_grid[y, x] = BURNING
                
                elif state == BURNING:
                    new_grid[y, x] = ASH

        self.grid = new_grid
        self.draw_grid()
        self.root.after(100, self.loop)

if __name__ == "__main__":
    root = tk.Tk()
    app = ForestFireWindModel(root)
    root.mainloop()