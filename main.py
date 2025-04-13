from Bio import SeqIO
from Bio.Restriction import EcoRI, BamHI, XhoI, SalI
from backend.primer_tools import generate_primers, check_tm_difference
from backend.ncbi_fetch import download_pUC18
from backend.vector_analysis import find_restriction_sites
from backend.restriction_sites import restriction_sites
from dna_features_viewer import GraphicFeature, GraphicRecord, CircularGraphicRecord
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import math

def load_insert(filepath):
    record = SeqIO.read(filepath, "fasta")
    seq = str(record.seq).upper()
    print(f"Загружена вставка ({len(seq)} bp): {seq[:30]}...")
    return seq

def visualize_insert(insert_seq, valid_pairs, output_path="primer_map_insert.png"):
    features = []
    insert_length = len(insert_seq)
    effective_len = 20

    features.append(GraphicFeature(start=0, end=insert_length, strand=0,
                                   color="#eeeeee", label="Insert"))

    color_map = {
        "EcoRI": "#1f77b4",
        "BamHI": "#ff7f0e",
        "XhoI": "#2ca02c",
        "SalI": "#d62728"
    }

    for idx, pair in enumerate(valid_pairs):
        enzyme1, forward, tm1, enzyme2, reverse, tm2 = pair
        enzyme1_name = enzyme1.__name__ if hasattr(enzyme1, "__name__") else str(enzyme1)
        enzyme2_name = enzyme2.__name__ if hasattr(enzyme2, "__name__") else str(enzyme2)

        forward_start = 10 * (idx + 1)
        forward_end = forward_start + effective_len
        reverse_end = insert_length - 10 * (idx + 1)
        reverse_start = reverse_end - effective_len

        features.append(GraphicFeature(start=forward_start, end=forward_end, strand=+1, level=idx + 1,
                                       color=color_map.get(enzyme1_name, "#cccccc"),
                                       label=f"F: {enzyme1_name} (Tm {tm1:.0f}°C)"))
        features.append(GraphicFeature(start=reverse_start, end=reverse_end, strand=-1, level=-(idx + 1),
                                       color=color_map.get(enzyme2_name, "#cccccc"),
                                       label=f"R: {enzyme2_name} (Tm {tm2:.0f}°C)"))

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))
    GraphicRecord(sequence_length=insert_length, features=features).plot(ax=ax1, with_ruler=True)
    ax1.set_title("Линейная карта вставки")
    CircularGraphicRecord(sequence_length=insert_length, features=features).plot(ax=ax2)
    ax2.set_title("Круговая карта вставки")
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    print(f"🧬 Визуализация вставки сохранена как {output_path}")
    return fig

def visualize_vector(insert_seq=None, used_enzymes=None, output_path="primer_map_vector.png"):
    pUC18_record = download_pUC18()
    if not pUC18_record:
        print("Ошибка: не удалось загрузить pUC18!")
        return

    vector_seq = pUC18_record.seq
    vector_length = len(vector_seq)
    features = []

    features.append(GraphicFeature(start=0, end=vector_length, strand=0,
                                   color="#eeeeee", label="pUC18"))

    if insert_seq:
        insertion_start = 100
        insertion_end = insertion_start + len(insert_seq)
        features.append(GraphicFeature(start=insertion_start, end=insertion_end, strand=+1,
                                       color="#ff9999", label="Вставка"))

    unique_sites = find_restriction_sites(vector_seq)
    for enzyme, sites in unique_sites.items():
        # Фильтруем по используемым рестриктазам (MCS)
        if used_enzymes and enzyme not in used_enzymes:
            continue
        # Извлекаем первое значение из списка позиций
        site = sites[0] if isinstance(sites, list) and sites else None
        if site is not None:
            # Здесь выводится только имя фермента, а данные из restriction_sites используются в работе
            features.append(GraphicFeature(start=site, end=site + 1, strand=0,
                                           color="#6699cc", label=str(enzyme)))

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))
    GraphicRecord(sequence_length=vector_length, features=features).plot(ax=ax1, with_ruler=True)
    ax1.set_title("Линейная карта вектора pUC18")
    CircularGraphicRecord(sequence_length=vector_length, features=features).plot(ax=ax2)
    ax2.set_title("Круговая карта вектора pUC18")
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    print(f"🧬 Визуализация вектора сохранена как {output_path}")
    return fig

def generate_pdf_report(insert_seq, valid_pairs, insert_fig, vector_fig, output_pdf="cloning_report.pdf"):
    with PdfPages(output_pdf) as pdf:
        plt.figure(figsize=(11.69, 8.27))
        plt.axis('off')
        plt.text(0.5, 0.75, "Отчёт по клонированию", ha='center', fontsize=24, weight='bold')
        plt.text(0.5, 0.6, f"Длина вставки: {len(insert_seq)} bp", ha='center', fontsize=14)
        plt.text(0.5, 0.5, f"Количество валидных пар праймеров: {len(valid_pairs)}", ha='center', fontsize=14)
        pdf.savefig()
        plt.close()

        pairs_per_page = 5
        num_pages = math.ceil(len(valid_pairs) / pairs_per_page)
        for page in range(num_pages):
            plt.figure(figsize=(11.69, 8.27))
            plt.axis('off')
            start = page * pairs_per_page
            end = min((page + 1) * pairs_per_page, len(valid_pairs))
            y = 0.95
            plt.text(0.05, y, "Валидные праймерные пары:", fontsize=14, weight='bold')
            y -= 0.05
            for idx in range(start, end):
                enzyme1, forward, tm1, enzyme2, reverse, tm2 = valid_pairs[idx]
                enzyme1_name = enzyme1.__name__ if hasattr(enzyme1, "__name__") else str(enzyme1)
                enzyme2_name = enzyme2.__name__ if hasattr(enzyme2, "__name__") else str(enzyme2)
                text = (
                    f"{enzyme1_name} / {enzyme2_name}\n"
                    f"  Forward: {forward} (Tm: {tm1:.1f}°C)\n"
                    f"  Reverse: {reverse} (Tm: {tm2:.1f}°C)\n"
                )
                plt.text(0.05, y, text, fontsize=10, family='monospace', va='top')
                y -= 0.18
            pdf.savefig()
            plt.close()

        pdf.savefig(insert_fig)
        plt.close(insert_fig)

        pdf.savefig(vector_fig)
        plt.close(vector_fig)

    print(f"✅ PDF-отчёт успешно сохранён как {output_pdf}")

def main():
    insert_seq = load_insert("data/example_insert.fasta")
    enzymes = [EcoRI, BamHI, XhoI, SalI]

    # Обновляем данные рестриктаз из словаря restriction_sites для использования при вставке
    for enzyme in enzymes:
        name = enzyme.__name__ if hasattr(enzyme, "__name__") else str(enzyme)
        if name in restriction_sites:
            enzyme.site = restriction_sites[name]['sequence']
            enzyme.cut_position = restriction_sites[name]['cut_position']
            enzyme.overhang = restriction_sites[name]['overhang']
            enzyme.type = restriction_sites[name]['type']

    restriction_pairs = [(enzymes[i], enzymes[j]) for i in range(len(enzymes)) for j in range(i + 1, len(enzymes))]

    raw_primers = generate_primers(insert_seq, restriction_pairs, min_homology=20)
    valid = check_tm_difference(raw_primers, max_tm_diff=5)

    print("\n=== Валидные праймерные пары ===\n")
    if valid:
        for enzyme1, forward, tm1, enzyme2, reverse, tm2 in valid:
            enzyme1_name = enzyme1.__name__ if hasattr(enzyme1, "__name__") else str(enzyme1)
            enzyme2_name = enzyme2.__name__ if hasattr(enzyme2, "__name__") else str(enzyme2)
            print(f"{enzyme1_name} - {enzyme2_name}")
            print(f"  Forward: {forward} (Tm: {tm1:.2f}°C)")
            print(f"  Reverse: {reverse} (Tm: {tm2:.2f}°C)\n")
    else:
        print("❗Нет валидных пар праймеров для визуализации вставки.")

    # Собираем перечень используемых рестриктаз для формирования MCS
    used_enzymes = set()
    if valid:
        for enzyme1, forward, tm1, enzyme2, reverse, tm2 in valid:
            used_enzymes.add(enzyme1.__name__ if hasattr(enzyme1, "__name__") else str(enzyme1))
            used_enzymes.add(enzyme2.__name__ if hasattr(enzyme2, "__name__") else str(enzyme2))

    insert_fig = visualize_insert(insert_seq, valid, output_path="primer_map_insert_both.png") if valid else None
    vector_fig = visualize_vector(insert_seq, used_enzymes=used_enzymes, output_path="primer_map_vector_both.png")

    if insert_fig and vector_fig:
        generate_pdf_report(insert_seq, valid, insert_fig, vector_fig, output_pdf="cloning_report.pdf")

if __name__ == "__main__":
    main()
