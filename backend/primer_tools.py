# backend/primer_tools.py

from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt
import random

def generate_primers(insert_seq, chosen_sites, min_homology=20, logger=print):
    """
    Генерирует праймеры для заданной вставки и пар сайтов рестрикции.
    Возвращает список кортежей:
    (enzyme1, forward_primer, forward_effective, enzyme2, reverse_primer, reverse_effective)
    """
    primers = []
    effective_length = min(len(insert_seq), min_homology)

    for enzyme1, enzyme2 in chosen_sites:
        recognition_seq1 = enzyme1.site
        recognition_seq2 = enzyme2.site

        # Пропускаем, если сайты уже есть во вставке
        if recognition_seq1 in insert_seq or recognition_seq2 in insert_seq:
            logger(f"Пропуск пары {enzyme1} и {enzyme2}: сайт уже есть во вставке.")
            continue

        forward_padding = ''.join(random.choices("AGCT", k=random.randint(1, 5)))
        reverse_padding = ''.join(random.choices("AGCT", k=random.randint(1, 5)))

        forward_effective = str(insert_seq[:effective_length])
        reverse_effective = str(Seq(insert_seq[-effective_length:]).reverse_complement())

        # Добавляем GC в начало при низком содержании GC
        if (forward_effective.count('G') + forward_effective.count('C')) < effective_length // 2:
            forward_effective = 'GC' + forward_effective
        if (reverse_effective.count('G') + reverse_effective.count('C')) < effective_length // 2:
            reverse_effective = 'GC' + reverse_effective

        forward_primer = f'{forward_padding}{recognition_seq1}{forward_effective}'
        reverse_primer = f'{reverse_padding}{recognition_seq2}{reverse_effective}'

        logger(f"Пара {enzyme1} - {enzyme2}:")
        logger(f"  Forward: padding={forward_padding}, site={recognition_seq1}, effective={forward_effective}")
        logger(f"  Reverse: padding={reverse_padding}, site={recognition_seq2}, effective={reverse_effective}")

        primers.append((enzyme1, forward_primer, forward_effective,
                        enzyme2, reverse_primer, reverse_effective))
    return primers


def calculate_tm(effective_seq):
    """
    Расчёт температуры плавления для эффективной части праймера.
    Использует упрощённую формулу для коротких последовательностей и
    Nearest-Neighbor метод для длинных.
    """
    effective_seq = Seq(effective_seq)
    length = len(effective_seq)
    if length <= 20:
        return 4 * (effective_seq.count('G') + effective_seq.count('C')) + 2 * (effective_seq.count('A') + effective_seq.count('T'))
    else:
        return mt.Tm_NN(effective_seq, dnac1=25, dnac2=25, Na=50, saltcorr=5)


def check_3prime_end(primer):
    """
    Проверка, что 3'-конец не заканчивается на 'T'
    """
    return primer[-1].upper() != 'T'


def check_hairpin(primer, min_stem=4, loop_size=3):
    """
    Проверка шпильки (hairpin) внутри праймера.
    Поиск комплементарных участков, разделённых петлёй (loop).
    """
    seq = primer.upper()
    L = len(seq)
    for i in range(L - (2 * min_stem + loop_size) + 1):
        max_stem = (L - i - loop_size) // 2
        for stem_len in range(min_stem, max_stem + 1):
            left = seq[i:i + stem_len]
            right_start = i + stem_len + loop_size
            right = seq[right_start:right_start + stem_len]
            if left == str(Seq(right).reverse_complement()):
                return True
    return False


def check_primer_dimer(forward, reverse, dimer_length=4):
    """
    Проверка димера между forward и reverse праймерами.
    Ищет комплементарность на 3'-концах.
    """
    f_end = forward[-dimer_length:].upper()
    r_end = reverse[-dimer_length:].upper()
    r_end_rc = str(Seq(r_end).reverse_complement())
    return f_end in r_end_rc


def check_secondary_structure(forward, reverse, logger=print):
    """
    Проверка на вторичные структуры:
    - Шпильки (hairpins)
    - Димеры между праймерами
    - T на 3'-конце (только предупреждение)
    """
    issues = False
    if not check_3prime_end(forward):
        logger("⚠️ Forward праймер заканчивается на T.")
    if not check_3prime_end(reverse):
        logger("⚠️ Reverse праймер заканчивается на T.")
    if check_hairpin(forward):
        logger("❌ Forward праймер образует шпильку.")
        issues = True
    if check_hairpin(reverse):
        logger("❌ Reverse праймер образует шпильку.")
        issues = True
    if check_primer_dimer(forward, reverse):
        logger("❌ Обнаружен праймерный димер между forward и reverse.")
        issues = True
    return not issues


def check_tm_difference(primers, max_tm_diff=5, logger=print):
    """
    Отбор праймеров с допустимой разницей температур плавления.
    Возвращает список только тех пар, где разница Tm не превышает max_tm_diff
    и отсутствуют вторичные структуры.
    """
    valid_pairs = []
    for enzyme1, forward, forward_eff, enzyme2, reverse, reverse_eff in primers:
        tm1 = calculate_tm(forward_eff)
        tm2 = calculate_tm(reverse_eff)

        logger(f"Пара {enzyme1} - {enzyme2}: Tm_forward = {tm1:.0f}°C, Tm_reverse = {tm2:.0f}°C")

        if abs(tm1 - tm2) <= max_tm_diff and check_secondary_structure(forward, reverse, logger):
            valid_pairs.append((enzyme1, forward, tm1, enzyme2, reverse, tm2))
        else:
            logger(f"❌ Пара {enzyme1}-{enzyme2} отклонена: ΔTm={abs(tm1 - tm2):.1f}°C или проблемы во вторичной структуре.")
    return valid_pairs
