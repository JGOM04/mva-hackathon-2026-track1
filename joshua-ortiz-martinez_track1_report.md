# Methods Report — MVA Hackathon 2026, Track 1: Variant Prediction

**Submitted by:** Joshua Gabriel Ortiz Martinez (individual submission, no team)
**Model:** 1 — `phenotype-vep-candidate-gene`
**Predictions file:** [`joshua-ortiz-martinez_phenotype-vep-candidate-gene.csv`](joshua-ortiz-martinez_phenotype-vep-candidate-gene.csv)


## 1. Summary

A phenotype-first, candidate-gene approach was used to narrow a proband
whole-genome VCF (~5.0M calls) down to a small, ranked set of variants. The
clinical phenotype (rhabdomyosarcoma, nephrocalcinosis, short stature, failure to
thrive, extreme prematurity/IUGR, and recurrent parental miscarriages) was matched
to a known inherited-cancer-predisposition/chromosomal-instability syndrome, a
panel of causal and differential genes for that syndrome was defined, and only the
variants overlapping those genes were pulled out of the VCF and functionally
annotated. This produced two high-confidence candidate variants in **BUB1B**,
submitted as a proposed compound-heterozygous pair (EPCR 0.85), plus one
lower-confidence secondary finding in **GINS1** (EPCR 0.15).

## 2. Phenotype interpretation

The clinical document (`data/Challenge_Clinical_Phenotype_1.docx`) lists, with HPO
terms:

| Feature | HPO ID |
|---|---|
| Rhabdomyosarcoma | HP:0002859 |
| Nephrocalcinosis | HP:0000121 |
| Short stature | HP:0004322 |
| Failure to thrive / skeletal muscle atrophy | HP:0001508 / HP:0003202 |
| Premature birth (32 weeks) | HP:0001622 |
| Very low birth weight (~1 kg, IUGR) | HP:0001518 |
| Parental recurrent miscarriages | HP:0200067 |

The document itself flags the combination of cancer predisposition + growth
restriction + adverse perinatal history + parental reproductive loss as pointing
to a **chromosomal instability disorder** rather than any single feature in
isolation. That cluster — early-onset rhabdomyosarcoma together with severe IUGR
and a mosaic/constitutional aneuploidy signature — is the classic presentation of
**Mosaic Variegated Aneuploidy syndrome (MVA)**, an autosomal-recessive condition
caused by biallelic loss-of-function variants in **BUB1B** (MVA type 1), and less
commonly in **CEP57** or **TRIP13**. (The dataset repository name itself,
`mva-hackathon-2026-data`, was read as a corroborating hint.) Recurrent parental
miscarriage is consistent with a heritable mechanism affecting chromosome
segregation.

## 3. Candidate gene panel

Because MVA is genetically heterogeneous and overlaps clinically with other
growth-restriction / genome-instability syndromes, the search was widened beyond
the three primary MVA genes to include the main differential diagnoses:

| Category | Genes |
|---|---|
| MVA syndrome (primary) | BUB1B, CEP57, TRIP13 |
| Seckel syndrome | ATR, ATRIP |
| Nijmegen breakage syndrome | NBN |
| Microcephalic osteodysplastic primordial dwarfism II | PCNT |
| MVA syndrome 2 / pre-RC gene disorders | MCM4, GINS1 |
| Meier-Gorlin syndrome | CDC6, ORC1, ORC4, ORC6, CDT1 |

GRCh38 gene coordinates were pulled from the Ensembl REST API
(`lookup/symbol/homo_sapiens/<GENE>`) and padded ±5 kb to catch splice/regulatory
variants near exon boundaries.

## 4. Pipeline

All scripts are in the repository root and run against the downloaded VCF
(`data/*.vcf.gz`, GRCh38, GATK HaplotypeCaller output, single proband sample,
`WGS_EX2312012`).

1. **`descargar_completo.py` / `instructions.py`** — download the challenge
   dataset from the SageBio Hugging Face repo (`SageBio/mva-hackathon-2026-data`).
2. **`visualizar.py`** — whole-VCF QC/overview pass: parses every record with a
   streaming `gzip` reader (no external VCF library), and reports:
   - total variant count and per-chromosome distribution
   - variant type (SNV / insertion / deletion), classified by REF/ALT length
   - zygosity (het/hom), parsed from the `GT` field
   - GATK hard-filter status (`PASS` vs. `MQ40`, `QD2`, `LowQual`, ...)

   Output: `variantes_resumen.png`.
3. **`candidatos_mva.py`** — targeted extraction + annotation:
   - streams the VCF once, keeping only records whose position falls inside a
     candidate-gene window
   - batches the surviving variants (up to 200 per call, the Ensembl VEP REST
     limit) to `POST https://rest.ensembl.org/vep/human/region` with
     `hgvs=1&canonical=1`, to get consequence, impact, canonical-transcript HGVS
     notation, and colocated-variant population frequencies (gnomAD genomes) in
     one pass, without installing a local VEP/ANNOVAR cache
   - for each variant, the canonical transcript's annotation for the target gene
     is used (falling back to the first available transcript if the canonical one
     isn't annotated for that gene), discarding neighboring genes incidentally
     captured by the ±5 kb window
   - keeps only `HIGH`/`MODERATE`-impact consequences (protein-truncating,
     missense, splice) — synonymous/intronic/UTR noise is dropped
4. **`visualizar_candidatos.py`** — evidence figure for the top candidates:
   allele-depth balance (supports genotype call quality), patient allele
   frequency vs. gnomAD population frequency (log scale, rarity check), and
   variant position on the BUB1B protein (needle/lollipop plot, 1050 aa,
   `ENSP00000287598`). Output: `variantes_candidatas_bub1b.png`.

EPCR values in the final submission were assigned manually from these outputs
(impact severity, population rarity, allele-depth balance, and whether a variant
plausibly completes a biallelic hit) — the pipeline itself is fully automated and
reproducible, but the final probability estimate is an expert judgment call, not
a fitted score.

## 5. Results

**Cohort-level QC** (`visualizar.py`, full VCF):

- 5,012,204 total variant records
- By type: 4,084,915 SNVs, 481,164 deletions, 446,125 insertions
- By zygosity: 3,113,937 heterozygous, 1,898,267 homozygous-alt
- By filter: 4,740,790 PASS (94.6%); remaining fail on mapping quality (`MQ40`,
  177,522), quality-by-depth (`QD2`, 35,402), or combinations thereof

**Candidate-gene search**: 718 variants fell inside the 13-gene panel window; 21
had `HIGH`/`MODERATE` functional impact. Of those, 18 were common population
polymorphisms (gnomAD AF 0.02–0.99 — far too frequent to explain a rare disease)
and were discarded. Three rare/novel variants remained:

| EPCR | Gene | Variant | Consequence | Zygosity | gnomAD AF | Finding type |
|---|---|---|---|---|---|---|
| 0.85 | BUB1B | chr15:40209701 T>G — c.2210T>G, p.Leu737Ter | stop_gained (HIGH) | het | 3.29e-05 | primary |
| 0.85 | BUB1B | chr15:40220612 T>G — c.3006T>G, p.Asn1002Lys | missense (MODERATE) | het | not in gnomAD | primary |
| 0.15 | GINS1 | chr20:25452451 TC>T — c.278del, p.Ser93fs*6 | frameshift (HIGH) | het | 1.05e-04 | secondary |

Full record in [`joshua-ortiz-martinez_phenotype-vep-candidate-gene.csv`](joshua-ortiz-martinez_phenotype-vep-candidate-gene.csv)
(the two BUB1B rows are submitted as a single compound-het row per the submission
template's `chrom_2/pos_2/ref_2/alt_2` columns).

### Primary call: BUB1B compound heterozygote

Both BUB1B variants are `PASS`, have well-balanced allele depth for a het call
(21 REF/25 ALT reads and 15 REF/13 ALT reads respectively — see
`variantes_candidatas_bub1b.png`, panel A), and are rare-to-absent in gnomAD
(panel B). A premature stop at residue 737 of a 1050-aa protein
(`p.Leu737Ter`, panel C) is expected to trigger nonsense-mediated decay or produce
a non-functional truncated protein; the second variant is a rare missense change
39 residues from the C-terminus. Two different loss-of-function/likely-damaging
alleles in the same recessive-disease gene, in a patient whose phenotype matches
that gene's associated syndrome, is the expected molecular signature of a
biallelic (compound-heterozygous) hit.

### Secondary/incidental finding: GINS1

A single heterozygous frameshift variant in GINS1 (also an autosomal-recessive
growth-restriction gene in the same differential) was found, but no second rare
variant in GINS1 was identified, so this is reported as a lower-confidence,
secondary finding — most plausibly incidental heterozygous carrier status rather
than a causal hit — per the challenge's allowance for secondary/incidental
findings.

## 6. Limitations

- **No parental/trio data.** The VCF contains a single proband sample, so the two
  BUB1B variants cannot be formally phased. They are ~11 kb apart — well beyond
  short-read insert size — so read-backed phasing from the FASTQ is not possible
  either. The compound-heterozygous interpretation is a strong hypothesis based
  on gene, consequence severity, and rarity, not a phased confirmation.
  Statistical phasing against a reference panel (e.g. Eagle2/SHAPEIT) or parental
  sequencing would resolve this.
- **Candidate-gene search, not a genome-wide scan.** `candidatos_mva.py` only
  inspects a 13-gene panel chosen from the phenotype. A causal variant outside
  this panel would not be found by this script; `visualizar.py`'s whole-VCF QC
  pass is the only genome-wide step, and it does not annotate consequence. This
  is a targeted approach chosen for time-boxed hackathon scope, not an unbiased
  genome-wide association — it trades sensitivity to variants outside the
  phenotype-driven gene panel for speed and precision on the genes most likely to
  explain this specific phenotype.
- **Annotation via live API, not a pinned local cache.** VEP annotation was
  pulled from the public Ensembl REST API at analysis time rather than a
  versioned local VEP cache, so results depend on Ensembl's current annotation
  build and are not perfectly pinned for long-term reproducibility.
- **No FASTQ-level QC in this pipeline.** Read/alignment quality (mapping rate,
  coverage uniformity, contamination) was not assessed from the raw FASTQ files;
  the analysis trusts the provided VCF's GATK quality fields.

## 7. Reproducibility

```
pip install -r requirements.txt
python descargar_completo.py      # downloads data/ from Hugging Face (needs HF_TOKEN in .env)
python visualizar.py              # whole-VCF QC -> variantes_resumen.png
python candidatos_mva.py          # candidate-gene search + VEP annotation -> console output
python visualizar_candidatos.py   # evidence figure -> variantes_candidatas_bub1b.png
```

`HF_TOKEN` must be set in a local `.env` file (see `.gitignore` — it is not
committed).

## 8. Method abstract (submission form Q13)

This submission uses a phenotype-first, candidate-gene strategy rather than an
unbiased genome-wide variant-calling/prioritization pipeline. The proband's
clinical features (early-onset rhabdomyosarcoma, severe intrauterine growth
restriction, short stature, nephrocalcinosis, and a family history of recurrent
miscarriage) were interpreted as a chromosomal-instability disorder, most
consistent with Mosaic Variegated Aneuploidy (MVA) syndrome. A panel of 13 genes
was assembled from MVA's primary causal genes (BUB1B, CEP57, TRIP13) and its main
clinical differentials (Seckel syndrome, Nijmegen breakage syndrome, MOPD II,
Meier-Gorlin syndrome). GRCh38 coordinates for each gene were retrieved from the
Ensembl REST API, and the proband's ~5-million-record VCF was streamed once to
extract the ~700 variants overlapping those gene windows. Each candidate was
annotated in batch via the Ensembl VEP REST API (consequence, impact, canonical
HGVS, and gnomAD population frequency), and variants were filtered to HIGH/MODERATE
protein-level impact and population rarity. This reduced 718 windowed variants to
3 plausible candidates. Two heterozygous BUB1B variants — a stop-gained
(p.Leu737Ter) and a rare missense (p.Asn1002Lys) — are proposed as a compound
heterozygous pair and submitted as the primary finding (EPCR 0.85); the gene,
consequence severity, rarity, and phenotype match are all consistent with a
biallelic BUB1B hit causing MVA1, though phase could not be confirmed from
proband-only short-read WGS. A single heterozygous GINS1 frameshift variant is
submitted as a lower-confidence secondary/incidental finding (EPCR 0.15), since
no second rare GINS1 hit was found to complete a biallelic call.

**Strengths:** the approach is fast (whole pipeline runs in minutes on a laptop,
no paid compute), fully scripted/reproducible, explicitly designed to output
compound-heterozygous pairs (not just single variants), and grounds every
candidate in orthogonal evidence (impact + rarity + allele-depth balance +
phenotype-gene match) rather than a single filter. **Limitations:** it is a
targeted, hypothesis-driven search rather than an unbiased genome-wide scan, so a
causal variant outside the 13-gene panel would be missed by the targeted script
(though `visualizar.py` still QCs the full VCF); phasing of the BUB1B pair is
inferred, not molecularly confirmed, for lack of trio or long-read data; and
annotation depends on a live third-party API rather than a version-pinned local
cache.
