# backend/vector_analysis.py

from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt
from Bio.SeqUtils import gc_fraction
from backend.restriction_sites import restriction_sites  
import re


def find_restriction_sites(sequence, allow_multiple=False, blunt_only=False):
    """
    Ищет сайты рестрикции по пользовательскому словарю.
    Возвращает словарь: {фермент: [позиции]}
    """
    results = {}

    for enzyme, props in restriction_sites.items():
        if blunt_only and props["type"] != "blunt":
            continue

        pattern = props["sequence"]
        positions = [m.start() for m in re.finditer(f"(?={pattern})", str(sequence))]
        if not positions:
            continue

        if not allow_multiple and len(positions) > 1:
            continue

        results[enzyme] = positions

    return results


def find_unique_enzyme_pairs(vector_seq, insert_seq, blunt_only=False):
    """
    Находит пары рестриктаз, у которых по одному уникальному сайту в векторе
    и отсутствуют сайты в вставке
    """
    vector_sites = find_restriction_sites(vector_seq, allow_multiple=False, blunt_only=blunt_only)
    insert_sites = find_restriction_sites(insert_seq, allow_multiple=True, blunt_only=blunt_only)

    insert_enzymes = set(insert_sites.keys())
    usable_pairs = []

    enzymes = list(vector_sites.keys())
    for i in range(len(enzymes)):
        for j in range(i + 1, len(enzymes)):
            e1, e2 = enzymes[i], enzymes[j]
            if e1 not in insert_enzymes and e2 not in insert_enzymes:
                usable_pairs.append((e1, e2))

    return usable_pairs


def get_enzyme_cut_details(enzyme_name):
    """
    Возвращает подробности о ферменте: позиция разреза, тип конца
    """
    return restriction_sites.get(enzyme_name, None)


def calculate_tm(primer_seq, na_conc=50, mg_conc=1.5, dNTPs=0.2):
    """
    Рассчитывает температуру отжига праймера по формуле ближайшего соседа
    """
    seq = Seq(primer_seq)
    return round(mt.Tm_NN(seq, Na=na_conc, Mg=mg_conc, dNTPs=dNTPs), 2)


def gc_content(primer_seq):
    """
    Возвращает GC-состав (%) праймера
    """
    return round(gc_fraction(primer_seq) * 100, 2)
