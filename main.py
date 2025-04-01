# main.py
import tkinter as tk
from frontend.gui import PrimerDesignerApp

def main():
    root = tk.Tk()
    try:
        root.iconbitmap("data/primer_icon.ico")
    except Exception as e:
        print(f"Ошибка загрузки иконки: {e}")
    app = PrimerDesignerApp(root)
    root.mainloop()

if __name__ == "__main__":
    main()
