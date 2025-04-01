# backend/vector_analysis.py
from Bio.Restriction import RestrictionBatch, EcoRI, HindIII, BamHI, PstI, SalI, SphI, XmaI

def find_restriction_sites(sequence, allow_multiple=False):
    """
    Ищет сайты рестрикции в заданной последовательности.
    Используются ферменты: EcoRI, HindIII, BamHI, PstI, SalI, SphI, XmaI.
    
    Если allow_multiple=False, возвращает для каждого фермента только первый сайт,
    при условии, что он встречается единожды.
    Если True – возвращает все найденные сайты.
    """
    restriction_enzymes = RestrictionBatch([EcoRI, HindIII, BamHI, PstI, SalI, SphI, XmaI])
    enzyme_sites = {enzyme: enzyme.search(sequence) for enzyme in restriction_enzymes}
    # Отфильтровываем ферменты, где сайты найдены
    enzyme_sites = {enzyme: sites for enzyme, sites in enzyme_sites.items() if sites}

    if not allow_multiple:
        return {enzyme: sites[0] for enzyme, sites in enzyme_sites.items() if len(sites) == 1}
    
    return enzyme_sites
