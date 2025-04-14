# backend/ncbi_fetch.py
import os
import logging
from Bio import Entrez, SeqIO

# Настройка логирования
logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")

Entrez.email = "reichsfurer1933@gmail.com"
DATA_FOLDER = "data"
PUC18_FILE = os.path.join(DATA_FOLDER, "pUC18.gb")

def download_pUC18():
    """
    Загружаем последовательность pUC18 из NCBI или используем локальный файл,
    если он уже сохранён.
    """
    genbank_id = "L08752.1"

    # Проверяем наличие локального файла
    if os.path.exists(PUC18_FILE):
        try:
            record = SeqIO.read(PUC18_FILE, "genbank")
            logging.info("Загружен локальный файл pUC18.")
            return record
        except Exception as e:
            logging.error(f"Ошибка при чтении локального файла pUC18: {e}. Будет произведена загрузка из NCBI.")

    # Если файла нет или произошла ошибка, выполняем загрузку из NCBI
    try:
        with Entrez.efetch(db="nucleotide", id=genbank_id, rettype="gb", retmode="text") as handle:
            pUC18_record = SeqIO.read(handle, "genbank")
        # Создаём каталог, если его нет
        if not os.path.exists(DATA_FOLDER):
            os.makedirs(DATA_FOLDER)
        SeqIO.write(pUC18_record, PUC18_FILE, "genbank")
        logging.info("Файл pUC18 успешно загружен и сохранён.")
        return pUC18_record
    except Exception as e:
        logging.error(f"Ошибка загрузки pUC18 из NCBI: {e}")
        return None
