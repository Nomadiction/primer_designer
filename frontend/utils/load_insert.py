# frontend/utils/load_insert.py
from Bio import SeqIO
import os

def load_insert(filepath):
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"Файл вставки не найден: {filepath}")
    record = SeqIO.read(filepath, "fasta")
    seq = str(record.seq).upper()
    print(f"Загружена вставка ({len(seq)} bp): {seq[:30]}...")
    return seq