# LIFE748 Assessment 2 — Analysis Scripts

Benchmarking genome assembly, annotation, machine learning and structural
prediction tools on a clinical *Escherichia coli* isolate (sample **GN3**).

**Author:** \<Faiza Khan\> (\<201936907\>)
**Module:** LIFE748.1 Scientific Report, Assessment 2 — August 2026

---

## What is here

| File | Part | What it does |
|---|---|---|
| `Part1_code.sh` | 1 | Assembly and annotation pipeline as executed: seqkit QC → Flye → minimap2 + Racon → SPAdes → Prokka → Barrnap |
| `assembly_stats.py` | 1 | Reference-free contiguity metrics (contigs, N50/L50/N90, GC, ambiguous bases) |
| `make_figure_part1.py` | 1 | Builds Figure A, the six-panel benchmarking figure |
| `Part2_code.qmd` | 2 | Quarto document: DESeq2 differential expression, clustering, and LR/LDA/SVM classification with nested feature selection |
| `gene_names.csv` | 2 | b-number → gene symbol map for *E. coli* K-12 MG1655 |
| `Part3_code.py` | 3 | AlphaFold model confidence analysis and structure-based DNA-binding site prediction |
| `compile_report.py` | — | Assembles the final report, embedding figures and result tables |

---

## Environment

Run on a MacBook Air (Apple Silicon, 8 GB RAM, 4 threads). `CONDA_SUBDIR=osx-64`
is required on Apple Silicon — bioconda publishes no native arm64 builds for
these tools, so they run under Rosetta 2.

```bash
CONDA_SUBDIR=osx-64 conda create -y -n asm \
    flye=2.9.3 spades=4.0.0 seqkit minimap2=2.26 racon=1.5.0

CONDA_SUBDIR=osx-64 conda create -y -n annot \
    prokka=1.14.6 barrnap=0.9

conda create -y -n rstats -c conda-forge -c bioconda \
    bioconductor-deseq2 r-ggplot2 r-ggrepel r-pheatmap r-caret \
    r-e1071 r-glmnet r-cluster r-mclust r-rcolorbrewer r-knitr r-kernlab

python3 -m pip install biotite numpy matplotlib
```

**Versions used:** Flye 2.9.3-b1797 · SPAdes 4.0.0 · minimap2 2.26 · Racon 1.5.0 ·
seqkit 2.8.2 · Prokka 1.14.6 · Barrnap 0.9 · R 4.5.3 · DESeq2 1.50.2 · Biotite 1.7.1

---

## Running the analysis

```bash
# Part 1 — assembly and annotation  (~25 min)
bash Part1_code.sh

# Part 2 — differential expression and machine learning  (~3 min)
conda activate rstats
Rscript -e 'knitr::purl("Part2_code.qmd", output="Part2_code.R", documentation=1)'
Rscript Part2_code.R

# Part 3 — structural bioinformatics  (~2 min)
# Requires the five AlphaFold models in structures/
for id in P0ACP1 Q57083 Q47129 P36673 P77743; do
  curl -s -O "https://alphafold.ebi.ac.uk/files/AF-${id}-F1-model_v6.pdb"
done
python3 Part3_code.py
```

**Note on paths:** SPAdes fails with a `std::filesystem::path` conversion error
if any input or output path contains a space. Work in a directory without them.

---

## Inputs

- `GN3_hifix30.fastq` — PacBio HiFi reads, 8,390 reads / 150.0 Mb / 32.3× coverage (not redistributed here)
- `dataset_LIFE748_1.csv` — count matrix, 4,464 genes × 50 samples (not redistributed here)
- AlphaFold DB v6 models: P0ACP1 (Cra), Q47129 (FeaR), P36673 (TreR), Q57083 (PerR), P77743 (PrpR)

---

## Key results

- **Assembly.** Flye: 4 circular contigs, 5,243,266 bp, N50 4.95 Mb, 3 min 18 s, 0.90 GB. SPAdes: 18 contigs, 8,028,811 bp — 1.53× the genome actually present.
- **Polishing.** Racon changed the assembly by 10 bp and the annotation not at all (4,865 → 4,848 CDS). HiFi consensus needs no external polishing.
- **Annotation.** 4,848 CDS, 90 tRNA, 22 rRNA (seven *rrn* operons resolved), 74.8% functional coverage.
- **Machine learning.** 77 DEGs. Hierarchical clustering 88% agreement. Ridge logistic regression 78.0 ± 9.2% under nested feature selection, versus 86.9 ± 7.4% when features are chosen from all samples — an 8.9-point selection bias.
- **Structure.** HTH motifs recovered 2.2–4.1× above chance across five transcription factors.

---

## Generative AI declaration

Generative AI (Anthropic Claude, model `claude-opus-5`, August 2026) was used to
assist in drafting, structuring and commenting these scripts, and to diagnose
runtime failures. All analyses were executed and verified by the author, and
every reported value derives from those runs. The full declaration and prompt
history is submitted with the report as `GAI_declaration_and_prompts.md`.
