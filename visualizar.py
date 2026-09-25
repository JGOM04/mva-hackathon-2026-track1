import gzip
from collections import Counter
import matplotlib.pyplot as plt

vcf_path = "data/WGS_EX2312012_HGWCNDSX7.vcf.gz"

chrom_counts = Counter()
tipo_counts = Counter()       # SNV / insercion / delecion
zigosidad_counts = Counter()  # heterocigota / homocigota_alt
filtro_counts = Counter()     # PASS vs. filtros de calidad fallidos
total_variants = 0
columnas = None

with gzip.open(vcf_path, "rt") as f:
    for linea in f:
        if linea.startswith("##"):
            continue
        if linea.startswith("#CHROM"):
            columnas = linea.strip().split("\t")
            continue

        campos = linea.rstrip("\n").split("\t")
        chrom, pos, vid, ref, alt, qual, filtro, info, fmt, muestra = campos[:10]

        chrom_counts[chrom] += 1
        total_variants += 1

        # PASS o el nombre del filtro que fallo (ej. MQ40)
        filtro_counts[filtro] += 1

        # Tipo de variante por longitud de REF/ALT (asume un solo ALT)
        primer_alt = alt.split(",")[0]
        if len(ref) == len(primer_alt) == 1:
            tipo_counts["SNV"] += 1
        elif len(ref) < len(primer_alt):
            tipo_counts["Insercion"] += 1
        else:
            tipo_counts["Delecion"] += 1

        # Zygosity a partir del genotipo (primer campo de FORMAT/muestra)
        gt = muestra.split(":")[0]
        alelos = gt.replace("|", "/").split("/")
        if len(set(alelos)) == 1:
            zigosidad_counts["Homocigota"] += 1
        else:
            zigosidad_counts["Heterocigota"] += 1

print(f"Total de variantes: {total_variants:,}")
print(f"Columnas del VCF: {columnas}")

print("\nVariantes por cromosoma (top 10):")
for chrom, count in chrom_counts.most_common(10):
    print(f"  {chrom}: {count:,}")

print("\nVariantes por tipo:")
for tipo, count in tipo_counts.most_common():
    print(f"  {tipo}: {count:,}")

print("\nZygosity:")
for z, count in zigosidad_counts.most_common():
    print(f"  {z}: {count:,}")

print("\nEstado de filtro (PASS = paso los filtros de calidad de GATK):")
for filtro, count in filtro_counts.most_common(10):
    print(f"  {filtro}: {count:,}")

# Solo graficamos cromosomas "normales" (1-22, X, Y, M)
cromosomas_principales = [str(i) for i in range(1, 23)] + ["X", "Y", "M"]
datos_chrom = {c: chrom_counts[c] for c in cromosomas_principales if c in chrom_counts}

fig, axes = plt.subplots(2, 2, figsize=(14, 10))

axes[0, 0].bar(datos_chrom.keys(), datos_chrom.values(), color="#4C72B0")
axes[0, 0].set_title("Variantes por cromosoma")
axes[0, 0].set_xlabel("Cromosoma")
axes[0, 0].set_ylabel("Numero de variantes")
axes[0, 0].tick_params(axis="x", rotation=90)

axes[0, 1].bar(tipo_counts.keys(), tipo_counts.values(), color="#DD8452")
axes[0, 1].set_title("Variantes por tipo")
axes[0, 1].set_ylabel("Numero de variantes")

axes[1, 0].bar(zigosidad_counts.keys(), zigosidad_counts.values(), color="#55A868")
axes[1, 0].set_title("Zygosity")
axes[1, 0].set_ylabel("Numero de variantes")

top_filtros = dict(filtro_counts.most_common(6))
axes[1, 1].bar(top_filtros.keys(), top_filtros.values(), color="#C44E52")
axes[1, 1].set_title("Estado de filtro GATK (top 6)")
axes[1, 1].set_ylabel("Numero de variantes")
axes[1, 1].tick_params(axis="x", rotation=45)

plt.tight_layout()
plt.savefig("variantes_resumen.png")
plt.show()
