import sys
import os
import streamlit as st
from Bio.Restriction import EcoRI, BamHI, XhoI, SalI

# Добавление пути к модулям
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from backend.restriction_sites import restriction_sites
from backend.primer_tools import generate_primers, check_tm_difference
from frontend.utils.load_insert import load_insert
from frontend.utils.report_generator import generate_pdf_report
from frontend.utils.visualize_insert import visualize_insert
from frontend.utils.visualize_vector import visualize_vector

# НАСТРОЙКА СТРАНИЦЫ
st.set_page_config(page_title="🔬 Помощник клонирования", layout="wide")

# ВЫБОР ТЕМЫ
theme = st.radio("🎨 Выберите тему", options=["Светлая", "Тёмная"], horizontal=True)

if theme == "Светлая":
    st.markdown("""
<style>
        .main { background-color: #111827; }
        .block-container { padding-top: 2rem; padding-bottom: 2rem; }
        h1, h2, h3, p {
            font-family: 'Segoe UI', sans-serif;
            color: #f9fafb;
        }
        .stButton>button {
            background-color: #10b981;
            color: white;
            border-radius: 10px;
            padding: 0.6em 1.4em;
            font-weight: 600;
            border: none;
        }
        .stButton>button:hover {
            background-color: #059669;
        }
        .stDownloadButton>button {
            background-color: #3b82f6;
            color: white;
            border-radius: 8px;
            padding: 0.6em 1.4em;
            font-weight: 600;
        }
        .stDownloadButton>button:hover {
            background-color: #2563eb;
        }
        .primer-box {
            background-color: #1f2937;
            border-radius: 8px;
            padding: 1em;
            margin-bottom: 1em;
            border-left: 4px solid #10b981;
            color: #f9fafb;
        }
    </style>
    """, unsafe_allow_html=True)
else:
    st.markdown("""
    <style>
        .main { background-color: #f9fafb; }
        .block-container { padding-top: 2rem; padding-bottom: 2rem; }
        h1, h2, h3, p {
            font-family: 'Segoe UI', sans-serif;
            color: #1f2937;
        }
        .stButton>button {
            background-color: #10b981;
            color: white;
            border-radius: 10px;
            padding: 0.6em 1.4em;
            font-weight: 600;
            border: none;
        }
        .stButton>button:hover {
            background-color: #059669;
        }
        .stDownloadButton>button {
            background-color: #3b82f6;
            color: white;
            border-radius: 8px;
            padding: 0.6em 1.4em;
            font-weight: 600;
        }
        .stDownloadButton>button:hover {
            background-color: #2563eb;
        }
        .primer-box {
            background-color: #f3f4f6;
            border-radius: 8px;
            padding: 1em;
            margin-bottom: 1em;
            border-left: 4px solid #10b981;
        }
    </style>
    """, unsafe_allow_html=True)

# ЗАГОЛОВОК
st.markdown("""
<h1>🔬 Помощник клонирования</h1>
<p style="font-size: 1.1em;">
    Умный инструмент для подбора праймеров, визуализации и генерации отчёта. Просто загрузите FASTA-файл и выберите рестриктазы!
</p>
""", unsafe_allow_html=True)

# ИНФОРМАЦИЯ О ПРАЙМЕРЕ
with st.expander("ℹ️ Информация о векторе pUC18 и правила подбора праймеров"):
    st.markdown(
        """
**🔬 О векторе pUC18**

`pUC18` — это синтетический плазмидный вектор, разработанный на основе вектора pBR322. Он широко используется в молекулярном клонировании благодаря своей высокой эффективности, простоте в работе и наличию полезных генетических элементов.

**Основные характеристики pUC18:**
- 📏 Размер: около **2686 пар оснований (bp)**
- 🧬 Репликон: **ColE1** — обеспечивает высокое число копий в *E. coli*
- 💊 **ampR**: ген устойчивости к **ампициллину** — используется для селекции трансформантов
- 🔷 **lacZα**: участок гена β-галактозидазы для **синего/белого скрининга** (вставка нарушает рамку считывания, колонии — белые)
- 🧩 **MCS (Multiple Cloning Site)**: мультиклонинг-сайт, включающий более 10 уникальных сайтов рестрикции, вставленный в рамку lacZα

**Преимущества использования pUC18:**
- 🔝 Высокое количество копий плазмиды (до 500–700 копий на клетку)
- 🔬 Удобство контроля трансформации (через лак-систему)
- ⚙️ Универсальность для различных методов клонирования
- 🧪 Совместимость с множеством рестриктаз и стандартных праймеров (например, M13-форвард/реверс)

---

**🧬 Структура праймера:**

Каждый праймер должен включать:
1. 🧷 Навеску (1–5 нт) — для эффективного связывания рестриктазы
2. ✂️ Сайт рестрикции (например, `GAATTC` для EcoRI)
3. 🔗 Гомологичный участок к вставке (минимум 16 нт, предпочтительно 20)

---

**🧠 Требования к вставке:**
- ❌ Вставка **не должна содержать сайты рестрикции**, используемые в праймерах
- ✔️ Выбранные сайты должны быть **уникальны в векторе** (не дублируются)
- 🌡️ Разница Tm между праймерами — **не более 3°C**
- 🔄 Обязательная проверка на:
  - Hairpin-структуры (внутримолекулярные петли)
  - Primer-dimer (взаимная комплементарность 3'-концов)
  - Нежелательные T на 3'-конце

---

**🧪 Поддерживаемые ферменты (примеры):**
`EcoRI`, `BamHI`, `XhoI`, `SalI`, `PstI`, `KpnI`, `XbaI`, `SmaI` и многие другие.  
Полный список доступен в модуле `restriction_sites.py`.

---

**📄 Генерируемые материалы:**
- 🖼️ Карта вставки с аннотацией праймеров (линейная и круговая)
- 📎 Карта вектора pUC18 с сайтами рестрикции
- 📄 PDF-отчёт с температурой, последовательностями, визуализацией
- 📧 Отправка email с вложениями отчёта (при необходимости)
        """,
        unsafe_allow_html=True
    )

# ЗАГРУЗКА ПОСЛЕДОВАТЕЛЬНОСТИ
st.subheader("📂 Загрузка последовательности")
uploaded_file = st.file_uploader("Загрузите файл вставки в формате FASTA", type=["fasta", "fa"])
express_insert = st.button("📂 Быстрая вставка и анализ")

insert_seq = None
if uploaded_file:
    with open("data/user_insert.fasta", "wb") as f:
        f.write(uploaded_file.getvalue())
    insert_seq = load_insert("data/user_insert.fasta")
    st.success("Файл вставки успешно загружен.")
elif express_insert:
    insert_seq = load_insert("data/example_insert.fasta")
    st.info("Загружен базовый пример вставки.")

# ВЫБОР РЕСТРИКТАЗ
st.subheader("🧪 Выбор рестриктаз")
enzymes_options = {
    'EcoRI': EcoRI,
    'BamHI': BamHI,
    'XhoI': XhoI,
    'SalI': SalI,
}
selected_enzymes = st.multiselect(
    "Выберите ферменты для клонирования",
    options=list(enzymes_options.keys()),
    default=['EcoRI', 'BamHI']
)

# АНАЛИЗ
if insert_seq and selected_enzymes:
    selected_enzymes_objs = [enzymes_options[enzyme] for enzyme in selected_enzymes]

    for enzyme in selected_enzymes_objs:
        name = enzyme.__name__ if hasattr(enzyme, "__name__") else str(enzyme)
        if name in restriction_sites:
            enzyme.site = restriction_sites[name]['sequence']
            enzyme.cut_position = restriction_sites[name]['cut_position']
            enzyme.overhang = restriction_sites[name]['overhang']
            enzyme.type = restriction_sites[name]['type']

    restriction_pairs = [
        (selected_enzymes_objs[i], selected_enzymes_objs[j])
        for i in range(len(selected_enzymes_objs))
        for j in range(i + 1, len(selected_enzymes_objs))
    ]

    st.subheader("🧬 Генерация праймеров")

    with st.spinner("Генерация праймеров и проверка температуры..."):
        raw_primers = generate_primers(insert_seq, restriction_pairs, min_homology=20)
        valid = check_tm_difference(raw_primers, max_tm_diff=5)

    if valid:
        with st.expander("📋 Валидные пары праймеров"):
            for enzyme1, forward, tm1, enzyme2, reverse, tm2 in valid:
                st.markdown(f"""
                    <div class=\"primer-box\">
                        <strong>{enzyme1.__name__} / {enzyme2.__name__}</strong><br>
                        🔹 Forward: <code>{forward}</code> (Tm: {tm1:.2f}°C)<br>
                        🔸 Reverse: <code>{reverse}</code> (Tm: {tm2:.2f}°C)
                    </div>
                """, unsafe_allow_html=True)

        st.subheader("🖼️ Визуализация")
        used_enzymes = {e.__name__ for e1, _, _, e2, _, _ in valid for e in (e1, e2)}
        insert_fig = visualize_insert(insert_seq, valid, output_path="primer_map_insert_both.png")
        vector_fig = visualize_vector(insert_seq=insert_seq, used_enzymes=used_enzymes, output_path="primer_map_vector_both.png")

        col1, col2 = st.columns(2)
        with col1:
            st.image("primer_map_insert_both.png", caption="🧬 Вставка с праймерами", use_column_width=True)
        with col2:
            st.image("primer_map_vector_both.png", caption="🧪 Вектор с сайтами рестрикции", use_column_width=True)

        st.subheader("📄 Генерация отчёта")
        generate_pdf_report(insert_seq, valid, insert_fig, vector_fig, output_pdf="cloning_report.pdf")
        with open("cloning_report.pdf", "rb") as pdf_file:
            st.download_button(
                label="📥 Скачать PDF-отчёт",
                data=pdf_file,
                file_name="cloning_report.pdf",
                mime="application/pdf",
            )
    else:
        st.warning("⚠️ Не удалось найти валидные пары праймеров. Попробуйте изменить вставку или рестриктазы.")
elif insert_seq and not selected_enzymes:
    st.info("Пожалуйста, выберите как минимум два фермента.")
else:
    st.info("Загрузите FASTA-файл или воспользуйтесь демонстрационной вставкой.")
