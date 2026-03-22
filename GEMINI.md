# Yleaf Pipeline — GEMINI.md
> Last Updated: 2026-03-22

## Обзор

Yleaf — инструмент предсказания Y-хромосомных гаплогрупп из WGS VCF/BAM.
Модифицированная версия Erasmus MC v3.2.1 с обогащением из YFull/FTDNA/YBrowse.

**Директория:** `/home/valalav/_dna/wgs/yleaf/Yleaf`

---

## Архитектура Данных

### Источники SNP маркеров

| Источник | Файл | SNPs | Обновление |
|----------|------|------|------------|
| **YFull** | `temp_update/current_tree.json` | ~427K | GitHub releases |
| **FTDNA** | `~/ystr-matcher/ftdna_haplo/data/get.json` | ~786K | FTDNA API |
| **YBrowse** | `~/wgs/ybrowse/snps_hg38.csv` | ~3.1M | Еженедельно |

### ⚠️ Деревья YFull и FTDNA НЕ совпадают по топологии

### Файлы данных Yleaf

| Файл | Назначение | Размер |
|------|-----------|--------|
| `yleaf/data/hg38/new_positions.txt` | SNP позиции hg38 | ~852K строк |
| `yleaf/data/hg_prediction_tables/tree.json` | Дерево гаплогрупп | ~192K nodes |

---

## Запуск

```bash
# BAM (рекомендуется — все позиции)
python3 -m yleaf.Yleaf -bam /path/to/sample.bam -o /tmp/output -rg hg38

# VCF (только варианты — меньше маркеров!)
python3 -m yleaf.Yleaf -vcf /path/to/sample.vcf.gz -o /tmp/output -rg hg38
```

> **ПРАВИЛО: НИКОГДА не сравнивать BAM и VCF прогоны!**
> BAM: ~324K маркеров, VCF: ~497. Разный input = разный результат.

---

## ✅ YFull Deepener — Post-processing углубление (2026-03-22)

### Архитектура

`yleaf/yfull_deepener.py` — post-processing модуль, углубляющий Yleaf предсказания
через YFull дерево + прямую проверку SNP в BAM (samtools mpileup).

```
Yleaf prediction → Deepener → Check SNPs in BAM → Deeper haplogroup
```

### Принцип работы

1. Получает Yleaf prediction (напр. `G-Y4464`)
2. Находит **все** ветки-потомки в YFull дереве (без лимита depth)
3. Для каждой ветки ищет адреса SNP из **3 источников** (по приоритету):
   - YFull CSV (`yfull_snp_p1-p658.csv`, confirmed=Yes) — ~515K SNP
   - Yleaf positions (`new_positions.txt`) — fallback, ~852K SNP
   - *(planned)* YBrowse CSV (`snps_hg38.csv`) — ~3.1M SNP
4. Проверяет позиции в BAM через `samtools mpileup`
5. Определяет derived/ancestral статус каждой ветки
6. **Итеративно** повторяет от найденной deeper ветки до терминальной
7. Transitive deepening: если потомок derived → все предки тоже derived

### Запуск

```bash
python3 yleaf/yfull_deepener.py \
  --prediction "G-Y4464" \
  --bam /path/to/sample.bam \
  --yfull-tree temp_update/current_tree.json \
  --csv /path/to/yfull_snp_p1-p658.csv \
  --positions yleaf/data/hg38/new_positions.txt
```

### Батч-результаты (76 образцов, 2026-03-22)

| Метрика | Значение |
|---------|----------|
| Углублено | **43/76 (56%)** |
| Без deeper | 33 |

### Известные ограничения

1. **T2T-only SNP**: ветки с SNP из T2T reference (CHM13) не имеют hg38 координат →
   невозможно проверить в hg38 BAM. Пример: R-Y518906 (4 SNP, все без hg38 addr)
2. **Гетерозиготные SNP**: смесь der/anc reads — считаются по большинству (>50%)
3. **Покрытие**: позиции с depth < min_depth (default=3) пропускаются

---

## ✅ YFull CSV SNP интеграция (2026-03-20/21)

- `yleaf/yfull_integrator.py` — фильтрация + интеграция YFull CSV SNP
- +82,030 маркеров (852K → 934K), +120 нод дерева
- BAM→BAM прогон: 3/4 образцов SAME, 1 CHANGED (не из-за YFull)

```bash
python3 -m yleaf.yfull_integrator \
    ~/wgs/yfull/data/yfull_snp_p1-p658.csv \
    --positions yleaf/data/hg38/new_positions.txt \
    --tree yleaf/data/hg_prediction_tables/tree.json \
    --yfull-tree ~/wgs/yfull/data/current_tree.json
```

> **YFull интеграция БЕЗОПАСНА.** +65K проверяемых маркеров, не ломает предсказания.

---

## Модифицированные файлы

| Файл | Изменение |
|------|-----------|
| `yleaf/predict_haplogroup.py` | Soft-fallback QC: если QC1 low, но QC2+QC3 ok → accept |
| `yleaf/Yleaf.py` | VCF index skip: 21с → 6.6с/образец |
| `yleaf/tree.py` | Итеративный DFS вместо рекурсии |
| `yleaf/updater.py` | Dual mapping: YFull primary + YBrowse fallback |
| `enrich_with_paths.py` | Post-processing: dual-tree paths через haplo server (порт 9003) |
| `yleaf/yfull_deepener.py` | **NEW**: YFull deepening через BAM |
| `yleaf/yfull_integrator.py` | **NEW**: YFull CSV SNP интеграция |

---

## Батч-результаты (852K, 2026-03-19)

| Метрика | Значение |
|---------|----------|
| Образцов | 76 |
| Определены | 66 (86.8%) |
| Low_Y_Signal | 10 |
| Среднее время | 6.6с |
| Распределение | G=33, J=17, R=9, E=2, I=2, C=2, Q=1 |

### Данные
- **BAM:** `/mnt/truenas-data/.../BAM/{SAMPLE}/{SAMPLE}.MGI.cutadapt.bwa.MarkDuplicates.bam`
- **VCF:** `/mnt/truenas-data/.../VCF/{SAMPLE}/{SAMPLE}.MGI.cutadapt.bwa.MarkDuplicates.DeepVariant.vcf.gz`

### Бэкапы
- `new_positions.txt.bak_before_yfull` — позиции до интеграции
- `tree.json.bak_before_yfull` — дерево до интеграции
