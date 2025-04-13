# restriction_sites.py

restriction_sites = {
    'EcoRI': {
        'sequence': 'GAATTC',
        'cut_position': (1, 5),
        'type': 'sticky',
        'overhang': 'AATT'
    },
    'Eco53kI': {
        'sequence': 'GAGCTC',
        'cut_position': (5, 1),
        'type': 'sticky',
        'overhang': 'GAGC'
    },
    'SacI': {
        'sequence': 'GAGCTC',
        'cut_position': (1, 5),
        'type': 'sticky',
        'overhang': 'GAGC'
    },
    'Acc65I': {
        'sequence': 'GGTACC',
        'cut_position': (1, 5),
        'type': 'sticky',
        'overhang': 'GGTAC'
    },
    'AvaI': {
        'sequence': 'CCCGGG',
        'cut_position': (1, 5),
        'type': 'sticky',
        'overhang': 'CCCGG'
    },
    'KpnI': {
        'sequence': 'GGTACC',
        'cut_position': (5, 1),
        'type': 'sticky',
        'overhang': 'GGTAC'
    },
    'XmaI': {
        'sequence': 'CCCGGG',
        'cut_position': (1, 5),
        'type': 'sticky',
        'overhang': 'CCCGG'
    },
    'BsoBI': {
        'sequence': 'CCCGGG',
        'cut_position': (2, 3),
        'type': 'sticky',
        'overhang': 'CC'
    },
    'TspMI': {
        'sequence': 'CCCGGG',
        'cut_position': (1, 5),
        'type': 'sticky',
        'overhang': 'CCCGG'
    },
    'SmaI': {
        'sequence': 'CCCGGG',
        'cut_position': (3, 3),
        'type': 'blunt',
        'overhang': ''
    },
    'BamHI': {
        'sequence': 'GGATCC',
        'cut_position': (1, 5),
        'type': 'sticky',
        'overhang': 'GATC'
    },
    'XbaI': {
        'sequence': 'TCTAGA',
        'cut_position': (1, 5),
        'type': 'sticky',
        'overhang': 'TCTAG'
    },
    'SalI': {
        'sequence': 'GTCGAC',
        'cut_position': (1, 5),
        'type': 'sticky',
        'overhang': 'GTCGAC'
    },
    'AccI': {
        'sequence': 'GTCGAC',
        'cut_position': (1, 5),
        'type': 'sticky',
        'overhang': 'GTCGAC'
    },
    'HincII': {
        'sequence': 'GTCGAC',
        'cut_position': (3, 3),
        'type': 'blunt',
        'overhang': ''
    },
    'PstI': {
        'sequence': 'CTGCAG',
        'cut_position': (5, 1),
        'type': 'sticky',
        'overhang': 'TGCAG'
    },
    'SbfI': {
        'sequence': 'CCTGCAGG',
        'cut_position': (6, 2),
        'type': 'sticky',
        'overhang': 'CCTGCA'
    },
    'SphI': {
        'sequence': 'GCATGC',
        'cut_position': (1, 5),
        'type': 'sticky',
        'overhang': 'GCATG'
    },
    'HindIII': {
        'sequence': 'AAGCTT',
        'cut_position': (1, 5),
        'type': 'sticky',
        'overhang': 'AAGCT'
    }
}
