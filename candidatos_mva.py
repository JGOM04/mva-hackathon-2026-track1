import gzip
import requests

vcf_path = "data/WGS_EX2312012_HGWCNDSX7.vcf.gz"

# Coordenadas GRCh38 (Ensembl), +/-5kb para capturar splice cercano al gen.
# Genes causales conocidos del sindrome de Aneuploidia Mosaico Variegada (MVA),
# que encaja con el fenotipo: restriccion de crecimiento + predisposicion a
# rabdomiosarcoma/tumor de Wilms + inestabilidad cromosomica.
GENES_CANDIDATOS = {
    # Sindrome MVA (herencia recesiva)
    "BUB1B": ("15", 40_160_984 - 5_000, 40_221_137 + 5_000),
    "CEP57": ("11", 95_789_965 - 5_000, 95_837_070 + 5_000),
    "TRIP13": ("5", 892_849 - 5_000, 919_357 + 5_000),
    # Diferencial: otros trastornos de restriccion de crecimiento +
    # inestabilidad cromosomica / microcefalia primordial
    "ATR": ("3", 142_449_007 - 5_000, 142_578_868 + 5_000),       # Seckel syndrome
    "ATRIP": ("3", 48_446_710 - 5_000, 48_467_645 + 5_000),       # Seckel syndrome 2
    "NBN": ("8", 89_924_515 - 5_000, 90_003_228 + 5_000),         # Nijmegen breakage syndrome
    "PCNT": ("21", 46_324_082 - 5_000, 46_445_769 + 5_000),       # MOPD II
    "MCM4": ("8", 47_960_185 - 5_000, 47_979_181 + 5_000),        # MVA syndrome 2
    "GINS1": ("20", 25_391_008 - 5_000, 25_452_700 + 5_000),      # NK cell + growth deficiency
    "CDC6": ("17", 40_287_831 - 5_000, 40_304_657 + 5_000),       # Meier-Gorlin syndrome
    "ORC1": ("1", 52_370_994 - 5_000, 52_404_492 + 5_000),        # Meier-Gorlin syndrome
    "ORC4": ("2", 147_930_396 - 5_000, 148_021_784 + 5_000),      # Meier-Gorlin syndrome
    "ORC6": ("16", 46_689_301 - 5_000, 46_698_399 + 5_000),       # Meier-Gorlin syndrome
    "CDT1": ("16", 88_803_777 - 5_000, 88_809_258 + 5_000),       # Meier-Gorlin syndrome
}

# Consecuencias que de verdad importan (alteran proteina o splicing)
IMPACTO_RELEVANTE = {"HIGH", "MODERATE"}

candidatas = []  # (gen, campos_vcf)

with gzip.open(vcf_path, "rt") as f:
    for linea in f:
        if linea.startswith("#"):
            continue
        campos = linea.rstrip("\n").split("\t")
        chrom, pos = campos[0], int(campos[1])

        for gen, (g_chrom, g_start, g_end) in GENES_CANDIDATOS.items():
            if chrom == g_chrom and g_start <= pos <= g_end:
                candidatas.append((gen, campos))

print(f"Variantes en ventana de los genes candidatos: {len(candidatas)}")

# Anotar con VEP (Ensembl) para quedarnos solo con lo funcionalmente relevante
variants_str = []
for gen, campos in candidatas:
    chrom, pos, vid, ref, alt = campos[0], campos[1], campos[2], campos[3], campos[4]
    primer_alt = alt.split(",")[0]
    variants_str.append(f"{chrom} {pos} . {ref} {primer_alt} . . .")

anotaciones = {}
LOTE = 200  # limite del endpoint POST de VEP
for i in range(0, len(variants_str), LOTE):
    lote = variants_str[i:i + LOTE]
    resp = requests.post(
        "https://rest.ensembl.org/vep/human/region",
        headers={"Content-Type": "application/json", "Accept": "application/json"},
        json={"variants": lote, "hgvs": 1, "canonical": 1},
        timeout=60,
    )
    resp.raise_for_status()
    for a in resp.json():
        anotaciones[a["input"]] = a

print("\n=== Variantes con impacto funcional relevante (HIGH/MODERATE) ===\n")
hallazgos = []
for (gen, campos), var_str in zip(candidatas, variants_str):
    a = anotaciones.get(var_str)
    if not a:
        continue

    chrom, pos, vid, ref, alt, qual, filtro, info, fmt, muestra = campos[:10]
    gt = muestra.split(":")[0]

    tcs = a.get("transcript_consequences", [])
    # Prioriza el transcript canonico; si no esta anotado, cae al primero que aplique
    tcs_gen = [tc for tc in tcs if tc.get("gene_symbol") == gen]
    tcs_gen.sort(key=lambda tc: tc.get("canonical") != 1)
    for tc in tcs_gen:
        if tc.get("impact") not in IMPACTO_RELEVANTE:
            continue

        # Frecuencia poblacional (gnomAD genomes, global) si esta anotada
        gnomad_af = None
        for cv in a.get("colocated_variants", []):
            freqs = cv.get("frequencies", {})
            for allele_freqs in freqs.values():
                if "gnomadg" in allele_freqs:
                    gnomad_af = allele_freqs["gnomadg"]

        hallazgos.append({
            "gen": gen,
            "pos": f"chr{chrom}:{pos}",
            "cambio": f"{ref}>{alt}",
            "gt": gt,
            "qual": qual,
            "filter": filtro,
            "consecuencia": ",".join(tc.get("consequence_terms", [])),
            "impacto": tc.get("impact"),
            "hgvsc": tc.get("hgvsc", "-"),
            "hgvsp": tc.get("hgvsp", "-"),
            "gnomad_af": gnomad_af,
        })
        break  # una anotacion representativa por variante basta

for h in hallazgos:
    rareza = "RARA/no vista" if h["gnomad_af"] is None else f"AF gnomAD={h['gnomad_af']}"
    print(f"[{h['gen']}] {h['pos']} {h['cambio']}  GT={h['gt']}  QUAL={h['qual']}  FILTER={h['filter']}")
    print(f"    {h['consecuencia']} ({h['impacto']})  {h['hgvsc']} {h['hgvsp']}  {rareza}")

if not hallazgos:
    print("Ninguna variante codificante/splice de impacto alto-moderado en estos 3 genes.")
    print("Sugerencia: revisar otros genes de la ruta (ATR, ATM, NBN, PCNT, MCM4) o bajar el umbral de impacto a LOW.")
