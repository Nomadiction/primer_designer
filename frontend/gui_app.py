import sys
import os
import streamlit as st
from Bio.Restriction import EcoRI, BamHI, XhoI, SalI

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from backend.restriction_sites import restriction_sites
from backend.primer_tools import generate_primers, check_tm_difference
from frontend.utils.load_insert import load_insert
from frontend.utils.report_generator import generate_pdf_report
from frontend.utils.visualize_insert import visualize_insert
from frontend.utils.visualize_vector import visualize_vector


# НАСТРОЙКА СТРАНИЦЫ
st.set_page_config(page_title="🔬 Cloning Assistant", layout="wide")

# ПЕРЕКЛЮЧАТЕЛЬ ТЕМЫ
theme = st.radio("🎨 Select Theme", options=["Light", "Dark"], horizontal=True)

if theme == "Light":
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
<h1>🔬 Cloning Assistant</h1>
<p style="font-size: 1.1em;">
    A smart tool for primer selection, visualization and report generation for your insert. Simply upload your FASTA file and select your enzymes!
</p>
""", unsafe_allow_html=True)

# ЗАГРУЗКА ПОСЛЕДОВАТЕЛЬНОСТИ
st.subheader("📂 Sequence loading")
uploaded_file = st.file_uploader("Download the FASTA insert file", type=["fasta", "fa"])
express_insert = st.button("📂 Express Insert")

insert_seq = None
if uploaded_file:
    with open("data/user_insert.fasta", "wb") as f:
        f.write(uploaded_file.getvalue())
    insert_seq = load_insert("data/user_insert.fasta")
    st.success("The insert file has been successfully uploaded.")
elif express_insert:
    insert_seq = load_insert("data/example_insert.fasta")
    st.info("Uploaded a quick data example")

# ВЫБОР РЕСТРИКТАЗ
st.subheader("🧪 Selection of restrictionases")
enzymes_options = {
    'EcoRI': EcoRI,
    'BamHI': BamHI,
    'XhoI': XhoI,
    'SalI': SalI,
}
selected_enzymes = st.multiselect(
    "Select the enzymes to be cloned",
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

    st.subheader("🧬 Primer generation")

    with st.spinner("Generating primers and validating Tm..."):
        raw_primers = generate_primers(insert_seq, restriction_pairs, min_homology=20)
        valid = check_tm_difference(raw_primers, max_tm_diff=5)

    if valid:
        with st.expander("📋 Valid primer pairs"):
            for enzyme1, forward, tm1, enzyme2, reverse, tm2 in valid:
                st.markdown(f"""
                    <div class="primer-box">
                        <strong>{enzyme1.__name__} / {enzyme2.__name__}</strong><br>
                        🔹 Forward: <code>{forward}</code> (Tm: {tm1:.2f}°C)<br>
                        🔸 Reverse: <code>{reverse}</code> (Tm: {tm2:.2f}°C)
                    </div>
                """, unsafe_allow_html=True)

        st.subheader("🖼️ Visualization")
        used_enzymes = {e.__name__ for e1, _, _, e2, _, _ in valid for e in (e1, e2)}
        insert_fig = visualize_insert(insert_seq, valid, output_path="primer_map_insert_both.png")
        vector_fig = visualize_vector(insert_seq=insert_seq, used_enzymes=used_enzymes, output_path="primer_map_vector_both.png")

        col1, col2 = st.columns(2)
        with col1:
            st.image("primer_map_insert_both.png", caption="🧬 Insert with primers", use_column_width=True)
        with col2:
            st.image("primer_map_vector_both.png", caption="🧪 Vector with restriction sites", use_column_width=True)

        st.subheader("📄 Report generation")
        generate_pdf_report(insert_seq, valid, insert_fig, vector_fig, output_pdf="cloning_report.pdf")
        with open("cloning_report.pdf", "rb") as pdf_file:
            st.download_button(
                label="📥 Download PDF report",
                data=pdf_file,
                file_name="cloning_report.pdf",
                mime="application/pdf",
            )
    else:
        st.warning("⚠️ No valid primer pairs were found. Try changing the insert or enzymes.")
elif insert_seq and not selected_enzymes:
    st.info("Please select at least two restriction enzymes.")
else:
    st.info("Download the FASTA file or use Express Insert.")
