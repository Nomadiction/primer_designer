# main.py
from Bio.Restriction import EcoRI, BamHI, XhoI, SalI
from backend.restriction_sites import restriction_sites
from backend.primer_tools import generate_primers, check_tm_difference
from frontend.utils.load_insert import load_insert
from frontend.utils.report_generator import generate_pdf_report
from frontend.utils.email_sender import send_ya_mail
from frontend.utils.visualize_insert import visualize_insert
from frontend.utils.visualize_vector import visualize_vector
from frontend.utils.email_sender import send_ya_mail


def main():
    insert_seq = load_insert("data/example_insert.fasta")
    enzymes = [EcoRI, BamHI, XhoI, SalI]

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

    used_enzymes = set()
    if valid:
        for enzyme1, _, _, enzyme2, _, _ in valid:
            used_enzymes.add(enzyme1.__name__)
            used_enzymes.add(enzyme2.__name__)

    insert_fig = visualize_insert(insert_seq, valid, output_path="primer_map_insert_both.png") if valid else None
    vector_fig = visualize_vector(insert_seq, used_enzymes=used_enzymes, output_path="primer_map_vector_both.png")
    
    # Генерация PDF
    generate_pdf_report(
        insert_seq, 
        valid, 
        insert_fig, 
        vector_fig, 
        output_pdf="cloning_report.pdf"
    )

    # Отправка письма
    send_ya_mail(
        recipients_emails="nik.volkov.vnv@gmail.com",
        insert_seq=insert_seq,
        valid_pairs=valid,
        insert_img_path="primer_map_insert_both.png",
        vector_img_path="primer_map_vector_both.png",
        file_path="cloning_report.pdf",
        attachment_filename="cloning_report.pdf"
    )

if __name__ == "__main__":
    main()
