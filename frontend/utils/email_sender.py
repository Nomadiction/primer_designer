import smtplib
import base64
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
from email.mime.base import MIMEBase
from email import encoders
from email.header import Header
from datetime import datetime

def encode_image_base64(path):
    with open(path, "rb") as img_file:
        encoded = base64.b64encode(img_file.read())
        return encoded.decode("ascii") 

def send_ya_mail(
    recipients_emails,
    insert_seq,
    valid_pairs,
    insert_img_path,
    vector_img_path,
    file_path,
    attachment_filename
):
    from email.mime.image import MIMEImage

    login = 'francesco.margaretti@yandex.com'
    password = 'ygfnzatiugdfbptt'

    msg = MIMEMultipart('related')
    msg['From'] = login
    msg['To'] = ', '.join(recipients_emails)
    msg['Subject'] = Header("🧬 Отчёт по клонированию", 'utf-8')

    # Альтернативный multipart — plain и html
    alt_part = MIMEMultipart('alternative')
    msg.attach(alt_part)

    # Текст
    alt_part.attach(MIMEText("Во вложении PDF-отчёт. Ниже HTML-версия с графикой.", 'plain', 'utf-8'))

    # CID для inline-картинок
    insert_cid = 'insert_map'
    vector_cid = 'vector_map'

    # HTML отчёт с cid-картинками
    html_report = create_html_report_with_cid(insert_seq, valid_pairs, insert_cid, vector_cid)
    alt_part.attach(MIMEText(html_report, 'html', 'utf-8'))

    # Вставка картинок как MIMEImage + CID
    if insert_img_path:
        with open(insert_img_path, 'rb') as f:
            img = MIMEImage(f.read(), _subtype="png")
            img.add_header('Content-ID', f'<{insert_cid}>')
            msg.attach(img)

    if vector_img_path:
        with open(vector_img_path, 'rb') as f:
            img = MIMEImage(f.read(), _subtype="png")
            img.add_header('Content-ID', f'<{vector_cid}>')
            msg.attach(img)

    # PDF вложение
    with open(file_path, 'rb') as attachment:
        part = MIMEBase('application', 'pdf')
        part.set_payload(attachment.read())
        encoders.encode_base64(part)
        part.add_header('Content-Disposition', f'attachment; filename="{attachment_filename}"')
        msg.attach(part)

    # SMTP отправка
    try:
        with smtplib.SMTP_SSL('smtp.yandex.ru', 465) as server:
            server.login(login, password)
            server.sendmail(login, recipients_emails, msg.as_string())
        print("✅ Письмо отправлено с inline-картинками.")
    except Exception as e:
        print("❌ Ошибка при отправке письма:", str(e))

def create_html_report_with_cid(insert_seq, valid_pairs, insert_cid, vector_cid):
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    html = f"""\
<html>
<head>
    <meta charset="utf-8">
    <style>
        body {{ font-family: Arial, sans-serif; padding: 20px; line-height: 1.5; }}
        h1 {{ text-align: center; }}
        .primer {{ font-family: monospace; background: #f5f5f5; padding: 10px; border-radius: 5px; }}
        .pair {{ margin-bottom: 20px; }}
    </style>
</head>
<body>
    <h1>🧬 Отчёт по клонированию</h1>
    <p><strong>📅 Дата:</strong> {timestamp}</p>
    <p><strong>📏 Длина вставки:</strong> {len(insert_seq)} bp</p>
    <p><strong>🔢 Валидных пар праймеров:</strong> {len(valid_pairs)}</p>
    <hr/>
"""

    if valid_pairs:
        html += "<h2>🔬 Валидные праймерные пары:</h2>"
        for enzyme1, forward, tm1, enzyme2, reverse, tm2 in valid_pairs:
            html += f"""
            <div class="pair">
                <div class="primer">
                    <b>{enzyme1.__name__}</b> / <b>{enzyme2.__name__}</b><br/>
                    Forward: {forward} (Tm: {tm1:.1f}°C)<br/>
                    Reverse: {reverse} (Tm: {tm2:.1f}°C)
                </div>
            </div>
            """

    html += f"""
    <h2>📌 Карта вставки</h2>
    <img src="cid:{insert_cid}" style="max-width:100%; border:1px solid #ccc;"/>

    <h2>📌 Карта вектора</h2>
    <img src="cid:{vector_cid}" style="max-width:100%; border:1px solid #ccc;"/>
</body></html>
"""
    return html
