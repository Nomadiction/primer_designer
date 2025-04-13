# frontend/utils/visualize_vector.py
from dna_features_viewer import GraphicFeature, GraphicRecord, CircularGraphicRecord
import matplotlib.pyplot as plt
from backend.ncbi_fetch import download_pUC18
from backend.vector_analysis import find_restriction_sites


def visualize_vector(insert_seq=None, chosen_site=None, used_enzymes=None, output_path="primer_map_vector.png"):
    """
    Визуализирует карту вектора pUC18 с учетом вставки и пересчета позиций элементов.
    
    Если указан chosen_site и найден соответствующий сайт рестрикции, 
    то вставка производится по его позиции, иначе вставка происходит в начало области MCS (по умолчанию 396 bp).
    
    Все элементы, находящиеся после позиции вставки, сдвигаются на длину вставляемой последовательности.
    При этом, если область MCS пересекается с точкой вставки, она делится на две части.
    """
    # Загружаем последовательность вектора pUC18
    pUC18_record = download_pUC18()
    if not pUC18_record:
        print("Ошибка: не удалось загрузить pUC18!")
        return

    vector_seq = pUC18_record.seq
    vector_length = len(vector_seq)
    new_length = vector_length  # если вставки нет, длина не меняется

    # Определяем область MCS для pUC18 
    mcs_start = 396
    mcs_end = 471

    # Получаем уникальные сайты рестрикции
    unique_sites = find_restriction_sites(vector_seq)

    # Если есть вставка, определяем позицию вставки и корректируем длину вектора
    if insert_seq:
        insert_length = len(insert_seq)
        if chosen_site and chosen_site in unique_sites and unique_sites[chosen_site]:
            # Если выбран сайт рестрикции найден, берем его первую позицию
            insertion_position = unique_sites[chosen_site][0]
        else:
            # По умолчанию вставляем в начало области MCS
            insertion_position = mcs_start

        new_length = vector_length + insert_length
        features = []

        # Добавляем первую часть вектора до вставки
        features.append(GraphicFeature(start=0, end=insertion_position, strand=0,
                                       color="#eeeeee", label="pUC18"))
        # Добавляем вторую часть вектора после вставки (с учетом сдвига)
        features.append(GraphicFeature(start=insertion_position + insert_length, end=new_length, strand=0,
                                       color="#eeeeee", label="pUC18"))
        # Добавляем саму вставку
        features.append(GraphicFeature(start=insertion_position,
                                       end=insertion_position + insert_length, strand=+1,
                                       color="#ff9999", label="Insertion"))
        # Если точка вставки попадает внутрь области MCS, делим MCS на две части
        if mcs_start < insertion_position < mcs_end:
            # Первая часть MCS до вставки
            features.append(GraphicFeature(start=mcs_start, end=insertion_position, strand=0,
                                           color="lightblue", label="MCS"))
            # Вторая часть MCS после вставки (сдвинутая)
            features.append(GraphicFeature(start=insertion_position + insert_length, end=mcs_end + insert_length, strand=0,
                                           color="lightblue"))
    else:
        # Если вставки нет, отображаем весь вектор
        features = [GraphicFeature(start=0, end=vector_length, strand=0, color="#eeeeee", label="pUC18")]

    # Добавляем сайты рестрикции, корректируя их позицию, если они находятся после вставки
    for enzyme, sites in unique_sites.items():
        if used_enzymes and enzyme not in used_enzymes:
            continue
        site = sites[0] if isinstance(sites, list) and sites else None
        if site is not None:
            if insert_seq and site >= insertion_position:
                adjusted_site = site + insert_length
            else:
                adjusted_site = site
            features.append(GraphicFeature(start=adjusted_site, end=adjusted_site + 1, strand=0,
                                           color="#6699cc", label=str(enzyme)))

    # Создаем два подграфика для визуализации: линейную и круговую карты
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))
    # Линейная карта
    GraphicRecord(sequence_length=new_length, features=features).plot(ax=ax1, with_ruler=True)
    ax1.set_title("Linear map of vector pUC18")
    # Круговая карта
    CircularGraphicRecord(sequence_length=new_length, features=features).plot(ax=ax2)
    ax2.set_title("Circular map of vector pUC18")
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    print(f"🧬 Визуализация вектора сохранена как {output_path}")
    return fig


