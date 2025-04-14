# frontend/utils/report_generator.py

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import math
import base64
from io import BytesIO
from frontend.utils.email_sender import send_ya_mail

def generate_pdf_report(
        insert_seq, 
        valid_pairs, 
        insert_fig=None, 
        vector_fig=None, 
        output_pdf="cloning_report.pdf",
        recipients_emails=None
):
    with PdfPages(output_pdf) as pdf:
        plt.figure(figsize=(11.69, 8.27))
        plt.axis('off')
        plt.text(0.5, 0.75, "Отчёт по клонированию", ha='center', fontsize=24, weight='bold')
        plt.text(0.5, 0.6, f"Длина вставки: {len(insert_seq)} bp", ha='center', fontsize=14)
        plt.text(0.5, 0.5, f"Количество валидных пар праймеров: {len(valid_pairs)}", ha='center', fontsize=14)
        pdf.savefig()
        plt.close()

        if valid_pairs:
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
                    enzyme1_name = enzyme1.__name__
                    enzyme2_name = enzyme2.__name__
                    text = (
                        f"{enzyme1_name} / {enzyme2_name}\n"
                        f"  Forward: {forward} (Tm: {tm1:.1f}°C)\n"
                        f"  Reverse: {reverse} (Tm: {tm2:.1f}°C)\n"
                    )
                    plt.text(0.05, y, text, fontsize=10, family='monospace', va='top')
                    y -= 0.18
                pdf.savefig()
                plt.close()

        if insert_fig:
            pdf.savefig(insert_fig)
            plt.close(insert_fig)

        if vector_fig:
            pdf.savefig(vector_fig)
            plt.close(vector_fig)

    send_ya_mail(
        recipients_emails="nik.volkov.vnv@gmail.com",
        insert_seq=insert_seq,
        valid_pairs=valid_pairs,
        insert_img_path="primer_map_insert_both.png",
        vector_img_path="primer_map_vector_both.png",
        file_path="cloning_report.pdf",
        attachment_filename="cloning_report.pdf"
    )

    print(f"✅ PDF-отчёт успешно сохранён как {output_pdf}")


def fig_to_base64(fig):
    buf = BytesIO()
    fig.savefig(buf, format='png', bbox_inches='tight')
    buf.seek(0)
    return base64.b64encode(buf.read()).decode('utf-8')

def create_html_report(insert_seq, valid_pairs, insert_img_data=None, vector_img_data=None):
    html = f"""
    <html>
    <head>
        <meta charset="utf-8">
        <style>
            body {{ font-family: sans-serif; padding: 20px; line-height: 1.5; }}
            h1 {{ text-align: center; }}
            pre {{ font-family: monospace; background: #f5f5f5; padding: 10px; border-radius: 5px; }}
            .pair {{ margin-bottom: 20px; }}
        </style>
    </head>
    <body>
        <h1>Отчёт по клонированию</h1>
        <p><strong>Длина вставки:</strong> {len(insert_seq)} п.н.</p>
        <p><strong>Количество валидных пар праймеров:</strong> {len(valid_pairs)}</p>
    """

    if valid_pairs:
        html += "<h2>Валидные праймерные пары</h2>"
        for enzyme1, forward, tm1, enzyme2, reverse, tm2 in valid_pairs:
            html += f"""
            <div class="pair">
                <strong>{enzyme1.__name__} / {enzyme2.__name__}</strong><br>
                <pre>Forward: {forward} (Tm: {tm1:.1f}°C)
Reverse: {reverse} (Tm: {tm2:.1f}°C)</pre>
            </div>
            """

    if insert_img_data:
        html += "<h2>Карта вставки</h2>"
        html += f'<img src="data:image/png;base64,{insert_img_data}" style="max-width:100%;"/>'

    if vector_img_data:
        html += "<h2>Карта вектора</h2>"
        html += f'<img src="data:image/png;base64,{vector_img_data}" style="max-width:100%;"/>'

    html += "</body></html>"
    return html
