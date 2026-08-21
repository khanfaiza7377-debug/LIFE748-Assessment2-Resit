#!/usr/bin/env bash
###############################################################################
# LIFE748 Assessment 2 — Part 1
# Genome assembly and annotation benchmarking of PacBio HiFi E. coli reads
#
# Sample : GN3_hifix30.fastq  (PacBio HiFi / CCS, ~32x E. coli coverage)
# Aim    : Benchmark Flye against SPAdes for assembly, quantify the effect of
#          consensus polishing, and annotate the result with Prokka + Barrnap.
#
# THIS IS THE SCRIPT AS ACTUALLY EXECUTED, on a MacBook Air (Apple Silicon,
# 8 GB RAM, 4 threads allocated) under conda environments built with
# CONDA_SUBDIR=osx-64. Runtimes and peak memory reported in the manuscript come
# from this run. Where a tool named in the assessment brief was not run (Bakta,
# QUAST, BUSCO), that is stated explicitly in section 8 rather than left
# implied by unused code.
#
# 
#
# -----------------------------------------------------------------------------
# GENERATIVE AI DECLARATION
# -----------------------------------------------------------------------------
# Generative AI (Anthropic Claude, model claude-opus-5, accessed August 2026)
# was used to assist in drafting, structuring and commenting this script, and
# to diagnose two failures during the run (a SPAdes crash caused by spaces in
# the working directory path, and the Flye consensus module aborting). All tool
# choices, parameters and the interpretation of results were reviewed, executed
# and verified by the author. The full prompt / interaction history is provided
# in GAI_declaration_and_prompts.md.
###############################################################################

set -euo pipefail     # stop on the first error rather than continuing silently

###############################################################################
# 0. ENVIRONMENT
#
# Two conda environments were used. CONDA_SUBDIR=osx-64 is required on Apple
# Silicon because bioconda does not publish native arm64 builds for these
# tools; they run under Rosetta 2 translation instead.
#
#   CONDA_SUBDIR=osx-64 conda create -y -n asm \
#       flye=2.9.3 spades=4.0.0 seqkit minimap2=2.26 racon=1.5.0
#   CONDA_SUBDIR=osx-64 conda create -y -n annot \
#       prokka=1.14.6 barrnap=0.9
#
# NOTE ON PATHS: SPAdes aborts with a std::filesystem::path conversion error if
# any input or output path contains a space. The data was therefore moved from
# "~/Downloads/life 748 bioinfo" to "~/life748" before the run. This cost one
# failed attempt and is recorded here so the next person does not repeat it.
###############################################################################
WORKDIR="${HOME}/life748"
READS="GN3_hifix30.fastq"
THREADS=4
MEM_GB=8
cd "${WORKDIR}"

# Resource measurement: macOS ships BSD `time`, whose -l flag reports wall-clock
# and maximum resident set size. (GNU `time -v` is not present by default.)
# Every timed command below is wrapped identically so the figures are comparable.
TIME="/usr/bin/time -l"

###############################################################################
# 1. INPUT READ QC
# Characterising the input first: assembly metrics are only interpretable
# relative to the depth and length of the reads that produced them.
###############################################################################
conda activate asm
seqkit stats -a "${READS}" | tee read_stats.txt

# Values obtained for GN3_hifix30.fastq:
#   8,390 reads | 150,001,705 bp | mean 17,878.6 bp | N50 17,329 bp
#   min 15,641 | max 39,668 | GC 50.53% | Q20 97.41% | Q30 93.56%
#   coverage vs 4,641,652 bp K-12 reference = 32.3x
# The Q20 figure is the evidence these are genuinely HiFi-quality reads, which
# is what justifies --pacbio-hifi and skipping read error correction below.

###############################################################################
# 2. ASSEMBLY WITH FLYE
#
# Flye is a repeat-graph assembler designed for long reads. --pacbio-hifi
# selects the low-error model, which skips the read-correction stage entirely.
#
# --iterations 0 disables Flye's internal polishing. This was NOT a stylistic
# choice: the consensus module of this build aborted mid-run, so polishing was
# performed externally with Racon (section 3) instead. Running it this way also
# yields an unpolished assembly for free, which is what makes the polishing
# comparison in section 6 possible.
###############################################################################
${TIME} flye --pacbio-hifi "${READS}" \
             --genome-size 4.6m \
             --threads "${THREADS}" \
             --iterations 0 \
             --out-dir flye_out

# Result: 198.19 s wall-clock, 963,178,496 B peak RSS (0.90 GB)
#         4 contigs, all circular; 5,243,266 bp; N50 4,947,123 bp
# Per-contig depth and circularity are in flye_out/assembly_info.txt:
#   contig_1  4,947,123 bp  29x  circular  <- chromosome
#   contig_4    158,288 bp  12x  circular  <- plasmid
#   contig_2     91,595 bp  16x  circular  <- plasmid
#   contig_3     46,260 bp  46x  circular  <- plasmid

###############################################################################
# 3. CONSENSUS POLISHING (minimap2 + Racon, two rounds)
#
# Racon corrects residual base-level errors by realigning the reads to the
# draft and recomputing consensus. Two rounds are standard; the second makes
# only marginal changes but confirms convergence.
#
# -x map-hifi selects minimap2's preset for low-error long reads.
###############################################################################
minimap2 -x map-hifi -a flye_out/assembly.fasta "${READS}" > aln1.sam
racon "${READS}" aln1.sam flye_out/assembly.fasta > polished1.fa

minimap2 -x map-hifi -a polished1.fa "${READS}" > aln2.sam
racon "${READS}" aln2.sam polished1.fa > polished2.fa

# Result: 5,243,276 bp — just 10 bp different from the unpolished assembly.
# By every contiguity metric the two are identical. Section 6 shows that this
# 10 bp difference nonetheless halves the number of genes Prokka calls.

###############################################################################
# 4. ASSEMBLY WITH SPAdes
#
# SPAdes is a de Bruijn graph assembler built for short reads. It has no
# long-read-only mode: --pacbio is designed to *supplement* an Illumina
# assembly, so supplying HiFi reads alone forces SPAdes to treat them as
# single-end (-s). This mismatch between tool design and data type is retained
# deliberately — quantifying the cost of using the wrong tool is the point of
# the benchmark, not an error in it.
#
# --only-assembler skips BayesHammer read correction, which models Illumina
#   error profiles and is inappropriate for Q20+ HiFi reads.
# --isolate is the recommended mode for high-coverage single-isolate bacteria.
###############################################################################
${TIME} spades.py --isolate --only-assembler \
                  -s "${READS}" \
                  -k 21,33,55 \
                  -t "${THREADS}" -m "${MEM_GB}" \
                  -o spades_out

# Result: 390.38 s wall-clock, 2,357,067,776 B peak RSS (2.19 GB)
#         18 contigs; 8,028,811 bp; N50 2,763,868 bp
# That is 1.53x the size of the genome actually present. The two largest
# contigs are 2,763,871 and 2,763,868 bp: the same locus emitted twice, because
# k <= 55 cannot span E. coli's ~5 kb rRNA operons and the graph stays tangled.

#--- 4b. Does a longer k-mer ladder help? ------------------------------------
# Tested explicitly rather than assumed. Prediction: no, because even k = 127
# is far shorter than the 5 kb repeats causing the duplication.
${TIME} spades.py --isolate --only-assembler \
                  -s "${READS}" \
                  -k 21,33,55,77,99,127 \
                  -t "${THREADS}" -m "${MEM_GB}" \
                  -o spades_big

###############################################################################
# 5. CONTIGUITY METRICS
#
# assembly_stats.py is a purpose-written script reproducing the reference-free
# subset of QUAST (contig count, total length, largest contig, N50/L50/N90, GC,
# ambiguous bases). It is included as supplementary material.
###############################################################################
python3 assembly_stats.py \
        Flye_unpolished=flye_out/assembly.fasta \
        Flye_polished=polished2.fa \
        SPAdes_k55=spades_out/contigs.fasta \
        SPAdes_k127=spades_big/contigs.fasta \
        | tee assembly_metrics.tsv

###############################################################################
# 6. ANNOTATION (Prokka) — run on BOTH assemblies
#
# Annotating the unpolished and polished assemblies with identical parameters
# isolates the effect of base-level accuracy. Any difference in gene counts is
# attributable to sequence quality alone, since the two assemblies differ by
# only 10 bp and have the same contig structure.
#
# --gcode 11 is the bacterial genetic code.
###############################################################################
conda activate annot

${TIME} prokka --outdir prokka_polished --prefix GN3 \
               --genus Escherichia --species coli --strain GN3 \
               --gcode 11 --cpus "${THREADS}" --force \
               polished2.fa

${TIME} prokka --outdir prokka_unpolished --prefix GN3u \
               --genus Escherichia --species coli --strain GN3 \
               --gcode 11 --cpus "${THREADS}" --force \
               flye_out/assembly.fasta

# rRNA genes are predicted separately: this Prokka build did not invoke Barrnap
# internally, so calling it directly is necessary to recover the rRNA counts.
barrnap polished2.fa > rrna.gff
echo "rRNA genes detected: $(grep -vc '^#' rrna.gff)"

#--- 6b. Functional annotation quality ---------------------------------------
# Two metrics beyond raw gene count, both of which matter downstream:
#   mean CDS length  - a direct readout of frameshift-driven ORF fragmentation
#   % hypothetical   - the fraction of CDS with no usable functional assignment
for d in prokka_polished/GN3 prokka_unpolished/GN3u; do
python3 -c "
L=[]; cur=0; hypo=0
for line in open('${d}.faa'):
    if line.startswith('>'):
        if cur: L.append(cur)
        cur = 0
        if 'hypothetical protein' in line: hypo += 1
    else:
        cur += len(line.strip())
if cur: L.append(cur)
print('${d}: n=%d  mean_aa=%.1f  hypothetical=%d (%.1f%%)'
      % (len(L), sum(L)/len(L), hypo, 100*hypo/len(L)))
"
done

###############################################################################
# 7. FIGURE
# make_figure_part1.py builds the six-panel benchmarking figure (Figure A) from
# the values produced above. Included as supplementary material.
###############################################################################
conda activate base
python3 make_figure_part1.py

###############################################################################
# 8. WHAT WAS NOT RUN, AND WHY
#
# The assessment brief names Bakta, QUAST and BUSCO. None were run, and the
# manuscript makes no claims about them:
#
#   Bakta  - the full database is ~30 GB and the light database still requires
#            a substantial download; neither was feasible within the available
#            time and disk budget. Prokka alone therefore carries the
#            annotation comparison, which is a genuine limitation of this
#            benchmark and is stated as such in the Conclusions.
#   QUAST  - replaced by assembly_stats.py for the reference-free metrics.
#            Reference-based metrics (misassemblies, genome fraction) were
#            consequently NOT computed, so this benchmark cannot distinguish a
#            contiguous assembly from a correct one on alignment evidence. The
#            circularity and size-ratio evidence used instead is weaker.
#   BUSCO  - single-copy orthologue completeness was not assessed. Mean CDS
#            length and rRNA operon count are used as proxies for annotation
#            completeness, which is not equivalent.
#
# Stating these omissions explicitly is preferable to shipping code that was
# never executed.
###############################################################################

echo "Part 1 complete. Key outputs:"
echo "  read_stats.txt              - input read QC"
echo "  flye_out/assembly_info.txt  - contig depth and circularity"
echo "  polished2.fa                - final polished assembly"
echo "  assembly_metrics.tsv        - contiguity metrics for all assemblies"
echo "  prokka_*/                   - annotation of both assemblies"
echo "  rrna.gff                    - Barrnap rRNA predictions"
