# backend/primer_tools.py

from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt
import random

def generate_primers(insert_seq, chosen_sites, min_homology=20, logger=print, min_len=18, max_len=35):
    """
    Генерация праймеров с сайтами рестрикции и контрольными структурами.
    Возвращает: [(enzyme1, forward_primer, forward_effective, enzyme2, reverse_primer, reverse_effective), ...]
    """
    primers = []
    effective_length = min(len(insert_seq), min_homology)

    for enzyme1, enzyme2 in chosen_sites:
        recognition_seq1 = enzyme1.site
        recognition_seq2 = enzyme2.site

        # Пропустить, если сайт уже присутствует во вставке
        if recognition_seq1 in insert_seq or recognition_seq2 in insert_seq:
            logger(f"⚠️ Пропущена пара {enzyme1} и {enzyme2}: сайт уже найден во вставке.")
            continue

        # Padding добавляет несколько случайных нуклеотидов для лучшей отрезки ферментом
        forward_padding = random_padding_with_gc()
        reverse_padding = random_padding_with_gc()

        forward_effective = str(insert_seq[:effective_length])
        reverse_effective = str(Seq(insert_seq[-effective_length:]).reverse_complement())

        # Добавление GC в начало эффективной части при недостатке GC
        forward_effective = balance_gc(forward_effective)
        reverse_effective = balance_gc(reverse_effective)

        forward_primer = f'{forward_padding}{recognition_seq1}{forward_effective}'
        reverse_primer = f'{reverse_padding}{recognition_seq2}{reverse_effective}'

        # Контроль длины
        if not (min_len <= len(forward_primer) <= max_len and min_len <= len(reverse_primer) <= max_len):
            logger(f"❌ Пропущена пара {enzyme1}-{enzyme2}: длина праймера вне допустимого диапазона.")
            continue

        logger(f"🔹 Пара {enzyme1} - {enzyme2}:")
        logger(f"  ➤ Forward: {forward_primer}")
        logger(f"  ➤ Reverse: {reverse_primer}")

        primers.append((enzyme1, forward_primer, forward_effective,
                        enzyme2, reverse_primer, reverse_effective))
    return primers


def random_padding_with_gc():
    """
    Генерирует случайную вставку из 1-5 нуклеотидов, обязательно включая хотя бы одну G или C
    """
    length = random.randint(1, 5)
    padding = ''.join(random.choices("AGCT", k=length))
    if 'G' not in padding and 'C' not in padding:
        padding = 'G' + padding[1:]
    return padding


def balance_gc(seq):
    """
    Добавляет GC в начало, если содержание GC ниже 50%
    """
    gc_count = seq.count('G') + seq.count('C')
    if gc_count < len(seq) // 2:
        return 'GC' + seq
    return seq


def calculate_tm(effective_seq):
    """
    Температура плавления эффективной части праймера.
    Короткие — упрощённая формула, длинные — Nearest Neighbor.
    """
    effective_seq = Seq(effective_seq)
    length = len(effective_seq)
    if length <= 20:
        return 4 * (effective_seq.count('G') + effective_seq.count('C')) + 2 * (effective_seq.count('A') + effective_seq.count('T'))
    else:
        return mt.Tm_NN(effective_seq, dnac1=25, dnac2=25, Na=50, saltcorr=5)


def check_3prime_end(primer):
    """
    Проверка, что 3'-конец НЕ заканчивается на T
    """
    return primer[-1].upper() != 'T'


def check_hairpin(primer, min_stem=4, loop_size=3):
    """
    Поиск шпильки (hairpin) в последовательности: обратнокомплементарные участки с петлёй
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
    Поиск димеров на 3’-концах forward и reverse праймеров
    """
    f_end = forward[-dimer_length:].upper()
    r_end = reverse[-dimer_length:].upper()
    r_end_rc = str(Seq(r_end).reverse_complement())
    return f_end in r_end_rc


def check_secondary_structure(forward, reverse, logger=print):
    """
    Проверка на вторичные структуры:
    - Hairpins
    - Primer-dimers
    - T на 3'-конце (не критично, но предупреждение)
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
    Отбор праймеров с допустимой ΔTm и без вторичных структур.
    Возвращает только валидные пары.
    """
    valid_pairs = []
    for enzyme1, forward, forward_eff, enzyme2, reverse, reverse_eff in primers:
        tm1 = calculate_tm(forward_eff)
        tm2 = calculate_tm(reverse_eff)
        delta_tm = abs(tm1 - tm2)

        logger(f"🔍 Пара {enzyme1} - {enzyme2}: Tm_forward = {tm1:.1f}°C, Tm_reverse = {tm2:.1f}°C")

        if delta_tm <= max_tm_diff and check_secondary_structure(forward, reverse, logger):
            logger("✅ Пара прошла проверку.\n")
            valid_pairs.append((enzyme1, forward, tm1, enzyme2, reverse, tm2))
        else:
            logger(f"❌ Пара отклонена: ΔTm={delta_tm:.1f}°C или проблемы структуры.\n")
    return valid_pairs
