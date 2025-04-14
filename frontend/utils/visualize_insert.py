# frontend/utils/visualize_insert.py

from dna_features_viewer import GraphicFeature, GraphicRecord, CircularGraphicRecord
import matplotlib.pyplot as plt

# Функция визуализации вставки с отображением праймеров
def visualize_insert(insert_seq, valid_pairs, output_path="primer_map_insert.png"):
    features = []   # Список объектов графических элементов 
    insert_length = len(insert_seq) # Длина вставляемой последовательности
    effective_len = 20  # Длина отображаемой части праймера 

    # Добавляем основную вставку как фоновую область
    features.append(GraphicFeature(
        start=0, end=insert_length, strand=0,  # Вся последовательность, без ориентации
        color="#eeeeee", label="Insert"       
    ))

    # Цвета для ферментов рестрикции
    color_map = {
        "EcoRI": "#1f77b4",
        "BamHI": "#ff7f0e",
        "XhoI": "#2ca02c",
        "SalI": "#d62728"
    }

    # Перебираем пары валидных праймеров и добавляем графические элементы
    for idx, pair in enumerate(valid_pairs):
        enzyme1, forward, tm1, enzyme2, reverse, tm2 = pair
        enzyme1_name = enzyme1.__name__  # Название фермента (функции) прямого праймера
        enzyme2_name = enzyme2.__name__  # Название фермента для обратного праймера

        # Смещённые координаты, чтобы праймеры не наслаивались
        forward_start = 10 * (idx + 1)
        forward_end = forward_start + effective_len
        reverse_end = insert_length - 10 * (idx + 1)
        reverse_start = reverse_end - effective_len

        # Добавляем прямой праймер
        features.append(GraphicFeature(
            start=forward_start, end=forward_end, strand=+1, level=idx + 1,
            color=color_map.get(enzyme1_name, "#cccccc"),  
            label=f"F: {enzyme1_name} (Tm {tm1:.0f}°C)"     # Метка с температурой отжига
        ))

        # Добавляем обратный праймер
        features.append(GraphicFeature(
            start=reverse_start, end=reverse_end, strand=-1, level=-(idx + 1),
            color=color_map.get(enzyme2_name, "#cccccc"),
            label=f"R: {enzyme2_name} (Tm {tm2:.0f}°C)"
        ))

    # Создаём фигуру с двумя графиками — линейный и круговой формат
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))

    # Линейная карта
    GraphicRecord(sequence_length=insert_length, features=features).plot(ax=ax1, with_ruler=True)
    ax1.set_title("Линейная карта вставки")

    # Круговая карта
    CircularGraphicRecord(sequence_length=insert_length, features=features).plot(ax=ax2)
    ax2.set_title("Круговая карта вставки")

    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    print(f"🧬 Визуализация вставки сохранена как {output_path}")

    return fig  