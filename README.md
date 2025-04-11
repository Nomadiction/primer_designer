<<<<<<< HEAD
# pUC18 Primer Designer

**pUC18 Primer Designer** — это программное обеспечение для автоматизированного подбора праймеров с рестрикционными сайтами для вектора pUC18. Приложение позволяет:
- Загружать последовательность вектора pUC18 из NCBI (с возможностью использования локального файла).
- Загружать последовательности вставок в формате FASTA.
- Анализировать вектор на наличие рестрикционных сайтов.
- Выбирать пары рестрикционных сайтов с учетом расстояния между ними (40–200 bp).
- Генерировать праймеры с навесками, включающими рестрикционные сайты, область гомологии и дополнительные нуклеотиды для корректной работы рестриктазы.
- Выполнять базовую оценку параметров праймеров: температуру плавления (Tm), GC-состав, наличие ошибок в 3'-конце, образование шпилек и димеров.
- Строить полярную карту вектора с обозначением позиций рестрикционных сайтов и вставленного фрагмента.

## Содержание

- [Требования](#требования)
- [Установка](#установка)
- [Структура проекта](#структура-проекта)
- [Описание модулей](#описание-модулей)
  - [backend/ncbi_fetch.py](#backendncbi_fetchpy)
  - [backend/primer_tools.py](#backendprimer_toolspy)
  - [backend/vector_analysis.py](#backendvector_analysispy)
  - [frontend/gui.py](#frontendguipy)
  - [main.py](#mainpy)
- [Использование](#использование)

## Требования

- Python 3.12+
- [Biopython](https://biopython.org/)
- [Tkinter](https://docs.python.org/3/library/tkinter.html)
- [Matplotlib](https://matplotlib.org/)
- [NumPy](https://numpy.org/)

Список зависимостей приведён в файле `requirements.txt`.

## Установка

1. Склонируйте репозиторий:
   ```bash
   git clone https://your-repo-link.git
   cd primer_designer
   ```

2. Установите зависимости:
   ```bash
   pip install -r requirements.txt
   ```

3. Убедитесь, что каталог `data` содержит:
   - Файл генбанк для pUC18 (`pUC18.gb`) (при отсутствии, он будет загружен автоматически).
   - Файл с последовательностью вставки (например, `example_insert.fasta`).

## Структура проекта

```
primer_designer/
│
├── main.py                  # Основной файл для запуска приложения.
├── README.md                # Документация проекта.
├── requirements.txt         # Файл со списком зависимостей.
├── structure.txt            # Описание структуры проекта (при наличии дополнительной информации).
│
├── backend/                 # Бэкэнд-логика проекта.
│   ├── ncbi_fetch.py        # Модуль для загрузки последовательности pUC18 из NCBI.
│   ├── primer_tools.py      # Модуль с инструментами для генерации и проверки праймеров.
│   └── vector_analysis.py   # Модуль для поиска рестрикционных сайтов в последовательности.
│
├── data/                    # Каталог с данными.
│   ├── example_insert.fasta # Пример последовательности вставки.
│   ├── primer_icon.ico      # Иконка для приложения.
│   └── pUC18.gb             # Файл с последовательностью вектора (будет создан/обновлён автоматически).
│
└── frontend/                # Фронтенд-приложение.
    └── gui.py             # Графический интерфейс, реализованный с помощью Tkinter.
```

## Описание модулей

### backend/ncbi_fetch.py

- **Функциональность:**  
  Загружает последовательность вектора pUC18 из NCBI или использует локально сохранённый файл.  
- **Ключевые функции:**
  - `download_pUC18()`: Проверяет наличие локального файла и при его отсутствии загружает последовательность с использованием Entrez API из Biopython.

### backend/primer_tools.py

- **Функциональность:**  
  Генерирует праймеры для вставки в вектор с учетом добавления навесок, рестрикционных сайтов и эффективной области гомологии.
- **Ключевые функции:**
  - `generate_primers(insert_seq, chosen_sites, min_homology=20, logger=print)`: Генерирует набор праймеров по заданной последовательности вставки и выбранным рестрикционным сайтам.
  - `calculate_tm(effective_seq)`: Вычисляет температуру плавления праймера, используя формулу для коротких последовательностей или метод NN для длинных.
  - `check_secondary_structure(...)`: Проверяет вторичную структуру праймеров (3'-концевой нуклеотид, образование шпилек, димеров).

### backend/vector_analysis.py

- **Функциональность:**  
  Ищет рестрикционные сайты в последовательности вектора с использованием набора ферментов (например, EcoRI, HindIII, BamHI, и др.).
- **Ключевые функции:**
  - `find_restriction_sites(sequence, allow_multiple=False)`: Ищет и возвращает позиции рестрикционных сайтов.

### frontend/gui.py

- **Функциональность:**  
  Реализует графический интерфейс приложения с использованием Tkinter. Позволяет пользователю:
  - Загружать последовательность вставки.
  - Выбирать рестрикционные сайты из найденных в векторе pUC18.
  - Анализировать и генерировать праймеры.
  - Отображать карту вектора с позицией вставки и рестрикционных сайтов.
- **Ключевые компоненты:**
  - Кнопки для загрузки файлов и анализа.
  - Текстовое поле для вывода результатов.
  - Интерактивный выбор рестрикционных сайтов.
  - Построение карты вектора с помощью Matplotlib.

### main.py

- **Функциональность:**  
  Точка входа в приложение. Запускает графический интерфейс, устанавливает иконку и инициализирует класс `PrimerDesignerApp`.

## Использование

1. **Запуск приложения:**
   ```bash
   python main.py
   ```

2. **Интерфейс:**
   - В главном окне вы можете:
     - Нажать кнопку **«Информация о векторе»** для получения справки о векторе pUC18.
     - Загрузить последовательность вставки из файла или из каталога `data`.
     - Нажать **«Анализировать»** для поиска рестрикционных сайтов и генерации праймеров.
     - После анализа можно нажать **«Показать карту вектора»** для визуального отображения позиций рестрикционных сайтов и вставленного фрагмента.

3. **Выбор рестрикционных сайтов:**
   - После загрузки последовательности вектора появляется окно для выбора рестрикционных сайтов.
   - Выберите нужные ферменты или используйте опцию «Выбрать всё».

4. **Генерация праймеров:**
   - Приложение автоматически сгенерирует праймеры, учитывая навески, эффективные области, GC-состав и разницу температур плавления.

5. **Отображение карты:**
   - По завершении анализа отображается полярная карта, где позиции рестрикционных сайтов помечены, а вставленный фрагмент выделен красной дугой.
=======
# primer_designer
>>>>>>> 72c0f83d1375143ea38ee105309f943e7b3ebafa
