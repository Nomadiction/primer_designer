# frontend/gui.py
import tkinter as tk
from tkinter import filedialog, messagebox, scrolledtext, simpledialog, ttk
import os
import time
import numpy as np
from Bio import SeqIO
from itertools import combinations
from backend.ncbi_fetch import download_pUC18
from backend.vector_analysis import find_restriction_sites
from backend.primer_tools import generate_primers, check_tm_difference
import matplotlib.pyplot as plt

# Каталог с данными
DATA_FOLDER = "data"

class PrimerDesignerApp:
    def __init__(self, root):
        self.root = root
        self.root.title("pUC18 Primer Designer")
        self.root.geometry("700x500")
        
        # Панель кнопок верхнего уровня
        top_frame = tk.Frame(root)
        top_frame.pack(pady=5)
        
        # Кнопка для отображения информации о векторе
        self.info_button = tk.Button(top_frame, text="Информация о векторе", command=self.show_vector_info)
        self.info_button.pack(side=tk.LEFT, padx=5)
        
        # Метка для выбора файла вставки
        tk.Label(root, text="Выберите последовательность вставки:").pack(pady=5)

        # Фрейм для кнопок загрузки и анализа
        button_frame = tk.Frame(root)
        button_frame.pack(pady=5)

        self.upload_button = tk.Button(button_frame, text="Загрузить файл", command=self.load_sequence)
        self.upload_button.pack(side=tk.LEFT, padx=5)

        self.auto_load_button = tk.Button(button_frame, text="Загрузить из папки", command=self.auto_load_sequence)
        self.auto_load_button.pack(side=tk.LEFT, padx=5)

        self.analyze_button = tk.Button(button_frame, text="Анализировать", command=self.analyze, state=tk.DISABLED)
        self.analyze_button.pack(side=tk.LEFT, padx=5)

        # Текстовое поле для вывода результатов анализа
        self.result_text = scrolledtext.ScrolledText(root, width=80, height=15)
        self.result_text.pack(pady=5)

        self.plot_button = tk.Button(root, text="Показать карту вектора", command=self.plot_vector, state=tk.DISABLED)
        self.plot_button.pack(pady=5)

        self.progress = ttk.Progressbar(root, orient="horizontal", length=300, mode="determinate")
        self.progress.pack(pady=10)

        # Переменные для хранения загруженной последовательности и выбранных параметров
        self.insert_seq = None
        self.selected_sites = {}
        # Сохранённая выбранная пара рестрикционных сайтов и координаты вставки для рекомбинантного вектора
        self.chosen_pair = None
        self.insertion_start = None
        self.insertion_end = None

    def show_vector_info(self):
            info_text = (
                "Информация о векторе pUC18:\n\n"
                "• pUC18 – классический плазмидный вектор, предназначенный для клонирования в бактериях (E. coli).\n"
                "• Используется для клонирования генов или фрагментов генома, имеет высокую копию и систему blue/white screening.\n"
                "• Вставляемые фрагменты получают дополнительные последовательности:\n"
                "   – Область гомологии (минимум 16 нуклеотидов),\n"
                "   – Рестрикционный сайт (для векторного клонирования),\n"
                "   – Навеска (1–5 нуклеотидов) для корректного распознавания рестриктазой.\n\n"
                "Перед клонированием обязательно проверяйте, чтобы выбранные рестрикционные сайты отсутствовали во вставляемой последовательности."
            )
            messagebox.showinfo("Информация о векторе pUC18", info_text)

    def slow_print(self, text):
        """Постепенно выводит текст в текстовое поле для имитации «живого» вывода."""
        for char in text:
            self.result_text.insert(tk.END, char)
            self.result_text.see(tk.END)
            self.root.update()
            time.sleep(0.01)

    def find_fasta_in_data(self):
        """
        Ищет файлы с расширением .fasta в каталоге DATA_FOLDER.
        Если найден один файл – возвращает его путь, иначе предлагает выбрать.
        """
        if not os.path.exists(DATA_FOLDER):
            os.makedirs(DATA_FOLDER)
        fasta_files = [f for f in os.listdir(DATA_FOLDER) if f.endswith(".fasta")]
        if not fasta_files:
            messagebox.showerror("Ошибка", "В папке 'data' нет файлов .fasta!")
            return None
        if len(fasta_files) == 1:
            return os.path.join(DATA_FOLDER, fasta_files[0])
        choice = simpledialog.askstring("Выбор файла",
                                        f"Найдено несколько файлов:\n{', '.join(fasta_files)}\nВведите имя файла:")
        if choice in fasta_files:
            return os.path.join(DATA_FOLDER, choice)
        else:
            messagebox.showerror("Ошибка", "Файл не найден!")
            return None

    def auto_load_sequence(self):
        """Автоматически загружает последовательность вставки из каталога DATA_FOLDER."""
        file_path = self.find_fasta_in_data()
        if file_path:
            self.load_sequence(file_path)

    def load_sequence(self, file_path=None):
        """
        Загружает последовательность вставки из выбранного файла в формате FASTA.
        При успешной загрузке активируется кнопка анализа.
        """
        if file_path is None:
            file_path = filedialog.askopenfilename(filetypes=[("FASTA files", "*.fasta"), ("All files", "*.*")])
        if file_path:
            try:
                self.insert_seq = SeqIO.read(file_path, "fasta").seq
                self.result_text.delete(1.0, tk.END)  # Очищаем предыдущее сообщение
                self.selected_sites = {}  # Сброс выбранных рестрикционных сайтов
                self.slow_print(f"🔬 Загружена последовательность: {len(self.insert_seq)} bp\n")
                self.analyze_button.config(state=tk.NORMAL)
            except Exception:
                messagebox.showerror("Ошибка", "Неверный формат файла!")

    def select_restriction_sites(self, unique_sites):
        """
        Открывает дополнительное окно для выбора рестрикционных сайтов.
        Пользователь может отметить отдельные ферменты или выбрать все сразу.
        """
        self.selected_sites = {}
        select_window = tk.Toplevel(self.root)
        select_window.title("Выберите сайты рестрикции")
        select_window.geometry("300x400")
        check_vars = {enzyme: tk.BooleanVar() for enzyme in unique_sites}

        def toggle_all():
            """Устанавливает или сбрасывает все чекбоксы."""
            state = select_all_var.get()
            for var in check_vars.values():
                var.set(state)

        def save_selection():
            """Сохраняет выбранные рестрикционные сайты и закрывает окно."""
            self.selected_sites = {enzyme: unique_sites[enzyme] for enzyme, var in check_vars.items() if var.get()}
            select_window.destroy()

        # Чекбокс "Выбрать всё"
        select_all_var = tk.BooleanVar()
        select_all_cb = tk.Checkbutton(select_window, text="Выбрать всё", variable=select_all_var, command=toggle_all)
        select_all_cb.pack(anchor="w", padx=10, pady=5)

        # Создаём чекбоксы для каждого фермента
        for enzyme, var in check_vars.items():
            chk = tk.Checkbutton(select_window, text=enzyme, variable=var)
            chk.pack(anchor='w', padx=20)
        btn_save = tk.Button(select_window, text="Выбрать", command=save_selection)
        btn_save.pack(pady=10)
        select_window.wait_window()

    def analyze(self):
        """
        Основной метод анализа:
         1. Загружает последовательность pUC18 (вектор) из NCBI.
         2. Находит рестрикционные сайты в векторе.
         3. Позволяет пользователю выбрать необходимые сайты.
         4. Формирует пары сайтов, удовлетворяющие условию расстояния.
         5. Генерирует праймеры и отбирает пары с допустимой разницей Tm.
        """
        if not self.insert_seq:
            messagebox.showerror("Ошибка", "Сначала загрузите последовательность!")
            return

        self.progress["value"] = 0
        self.root.update()
        self.progress["value"] = 25
        self.root.update()

        pUC18_record = download_pUC18()
        if not pUC18_record:
            messagebox.showerror("Ошибка", "Не удалось загрузить pUC18!")
            return

        self.progress["value"] = 45
        self.root.update()

        unique_sites = find_restriction_sites(pUC18_record.seq)
        self.select_restriction_sites(unique_sites)
        if not self.selected_sites:
            self.slow_print("\n❌ Вы не выбрали сайты рестрикции!\n")
            return

        self.progress["value"] = 65
        self.root.update()

        # Отбираем пары сайтов, расстояние между которыми находится в диапазоне 40–200 bp
        selected_pairs = [
            ((enzyme1, site1), (enzyme2, site2), abs(site1 - site2))
            for (enzyme1, site1), (enzyme2, site2) in combinations(self.selected_sites.items(), 2)
            if 40 <= abs(site1 - site2) <= 200
        ]
        self.slow_print("\n🔬 Анализ pUC18 завершён.\n")
        if selected_pairs:
            self.slow_print("\n✅ Выбранные сайты рестрикции:\n")
            for (enzyme1, site1), (enzyme2, site2), distance in selected_pairs:
                self.slow_print(f"🔹 {enzyme1} на позиции {site1}\n")
                self.slow_print(f"🔹 {enzyme2} на позиции {site2}\n")
                self.slow_print(f"📏 Расстояние между сайтами: {distance} bp\n\n")
            self.progress["value"] = 85
            self.root.update()

            primers = generate_primers(self.insert_seq, selected_pairs)
            valid_primers = check_tm_difference(primers, max_tm_diff=3)

            self.slow_print("\n✅ Праймеры с допустимой разницей Tm:\n")
            for enzyme1, forward, tm1, enzyme2, reverse, tm2 in valid_primers:
                self.slow_print(f"🔹 Прямой ({enzyme1}): {forward} (Tm: {tm1}°C)\n")
                self.slow_print(f"🔹 Обратный ({enzyme2}): {reverse} (Tm: {tm2}°C)\n")
            self.plot_button.config(state=tk.NORMAL)
        else:
            self.slow_print("\n❌ Нет подходящих сайтов рестрикции!\n")
            messagebox.showerror("Ошибка", "Не удалось найти подходящую пару сайтов.")

        self.progress["value"] = 100
        self.root.update()

    def plot_vector(self):
        """
        Строит полярную карту вектора pUC18.
        Если загружена вставка и выбраны сайты, формируется рекомбинантный вектор:
          участок между выбранными позициями заменяется на вставку.
        Если после рекомбинации сайты рестрикции не найдены, используется исходный вектор.
        Вставленный фрагмент выделяется красной дугой.
        """
        pUC18_record = download_pUC18()
        if not pUC18_record:
            messagebox.showerror("Ошибка", "Не удалось загрузить pUC18!")
            return

        vector_seq = pUC18_record.seq

        if self.insert_seq and self.selected_sites:
            chosen_sites = list(self.selected_sites.values())
            insertion_start = min(chosen_sites)
            insertion_end = max(chosen_sites)
            new_vector = vector_seq[:insertion_start] + self.insert_seq + vector_seq[insertion_end:]
            insert_start_new = insertion_start
            insert_end_new = insertion_start + len(self.insert_seq)
        else:
            new_vector = vector_seq
            insert_start_new = None

        # Пробуем найти рестрикционные сайты в новом векторе
        unique_sites = find_restriction_sites(new_vector)
        if not unique_sites:
            # Если не найдены, используем сайты из исходного вектора
            unique_sites = find_restriction_sites(vector_seq)
            if not unique_sites:
                messagebox.showerror("Ошибка", "Сайты рестрикции не найдены!")
                return

        sequence_length = len(new_vector)
        fig, ax = plt.subplots(figsize=(7, 7), subplot_kw={'projection': 'polar'})
        theta = np.linspace(0, 2 * np.pi, 100)
        ax.plot(theta, [1] * 100, color="black", linewidth=2)
        colors = plt.cm.get_cmap("tab10", len(unique_sites))
        for i, (enzyme, site) in enumerate(sorted(unique_sites.items(), key=lambda x: x[1])):
            angle = (site / sequence_length) * 2 * np.pi
            ax.scatter(angle, 1, color=colors(i), s=100, label=f"{enzyme}")

        if insert_start_new is not None:
            angle_start = (insert_start_new / sequence_length) * 2 * np.pi
            angle_end = (insert_end_new / sequence_length) * 2 * np.pi
            ax.plot([angle_start, angle_end], [1, 1], color="red", linewidth=6, label="Вставка")

        ax.legend(loc="lower center", bbox_to_anchor=(0.5, -0.2), ncol=3, fontsize=10)
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax.set_yticks([])
        ax.spines['polar'].set_visible(False)
        plt.title("Карта вектора pUC18", fontsize=14)
        plt.show()

if __name__ == "__main__":
    root = tk.Tk()
    app = PrimerDesignerApp(root)
    root.mainloop()

