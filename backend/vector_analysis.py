# backend/vector_analysis.py

from Bio.Restriction import AllEnzymes, RestrictionBatch
from Bio.SeqUtils import MeltingTemp as mt
from Bio.SeqUtils import gc_fraction
from Bio.Seq import Seq

def get_supported_enzymes(blunt_only=False):
    """
    Возвращает список поддерживаемых рестриктаз.
    Если blunt_only=True, возвращаются только ферменты с тупыми концами.
    """
    enzymes = [e for e in AllEnzymes]
    if blunt_only:
        enzymes = [e for e in enzymes if hasattr(e, 'is_blunt') and e.is_blunt()]
    return enzymes


def find_restriction_sites(sequence, allow_multiple=False, blunt_only=False):
    """
    Ищет сайты рестрикции в последовательности.
    
    Параметры:
    - allow_multiple: если False, оставляет только ферменты с 1 сайтом.
    - blunt_only: если True, фильтрует по тупым концам.
    
    Возвращает словарь: {'EcoRI': [позиции]}
    """
    enzymes = RestrictionBatch(get_supported_enzymes(blunt_only))
    enzyme_sites = {enzyme.__name__: enzyme.search(sequence) for enzyme in enzymes}
    enzyme_sites = {name: sites for name, sites in enzyme_sites.items() if sites}

    if not allow_multiple:
        return {name: sites for name, sites in enzyme_sites.items() if len(sites) == 1}
    
    return enzyme_sites


def find_unique_enzyme_pairs(vector_seq, insert_seq, blunt_only=False):
    """
    Находит пары рестриктаз, у которых по одному уникальному сайту в векторе
    и отсутствуют сайты в вставке.

    Возвращает: [('EcoRI', 'HindIII'), ...]
    """
    vector_sites = find_restriction_sites(vector_seq, allow_multiple=False, blunt_only=blunt_only)
    insert_sites = find_restriction_sites(insert_seq, allow_multiple=True, blunt_only=blunt_only)

    insert_enzymes = set(insert_sites.keys())
    usable = []

    enzymes = list(vector_sites.keys())
    for i in range(len(enzymes)):
        for j in range(i + 1, len(enzymes)):
            e1, e2 = enzymes[i], enzymes[j]
            if e1 not in insert_enzymes and e2 not in insert_enzymes:
                usable.append((e1, e2))
    
    return usable


def calculate_tm(primer_seq, na_conc=50, mg_conc=1.5, dNTPs=0.2):
    """
    Рассчитывает температуру отжига праймера.
    Используется Wallace rule (по умолчанию).
    """
    seq = Seq(primer_seq)
    return round(mt.Tm_NN(seq, Na=na_conc, Mg=mg_conc, dNTPs=dNTPs), 2)


def gc_content(primer_seq):
    """
    Возвращает GC-состав (%) праймера.
    """
    # Умножение результата на 100 для преобразования в проценты
    return round(gc_fraction(primer_seq) * 100, 2)

