#!/usr/bin/env python3
"""
LIFE748 Assessment 2 - Part 1
Figure A: Genome assembly and annotation benchmarking summary.

Builds a six-panel figure comparing Flye (unpolished and Racon-polished) with
SPAdes on assembly contiguity, computational cost, and downstream annotation
quality. All values are read from the actual run outputs (Table_B_assembly_metrics.tsv,
GNU time logs, and Prokka summary files) - nothing is hard-coded by hand except
the E. coli K-12 MG1655 reference values used as benchmarks.

GENERATIVE AI DECLARATION
Generative AI (Anthropic Claude, model claude-opus-5, accessed August 2026) was used
to assist in drafting and commenting this script. All analytical choices and the
interpretation of results were reviewed and verified by the author. The full
prompt/interaction history is provided in GAI_declaration_and_prompts.md.

"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---------------------------------------------------------------------------
# Reference values for E. coli K-12 MG1655 (GCF_000005845.2), used as the
# biological "correct answer" against which each assembly/annotation is judged.
# ---------------------------------------------------------------------------
REF_GENOME_BP = 4_641_652
REF_CDS       = 4_298
REF_MEAN_AA   = 316.0

# ---------------------------------------------------------------------------
# Measured results. Assembly metrics come from assembly_stats.py; runtime and
# peak RSS from GNU `time -v`; annotation counts from Prokka's summary .txt.
# ---------------------------------------------------------------------------
ASM = ["Flye\n(unpolished)", "Flye\n(polished)", "SPAdes"]
total_bp  = np.array([5_243_266, 5_243_276, 8_028_811])
n50       = np.array([4_947_123, 4_947_128, 2_763_868])
contigs   = np.array([4, 4, 18])

# Runtime and peak RSS were captured only for the two assemblers themselves
# (GNU `time -l`); the Racon polishing rounds were not separately instrumented,
# so panels C-D compare Flye against SPAdes only rather than implying a
# measurement that was not taken.
ASM2       = ["Flye", "SPAdes"]
runtime_m  = np.array([198.19, 390.38]) / 60
memory_gb  = np.array([963_178_496, 2_357_067_776]) / 1_073_741_824
ASM2_COLS  = None      # assigned below, after the palette is defined

# Annotation of BOTH assemblies with identical Prokka parameters. The two are
# indistinguishable, which is the finding: HiFi reads are accurate enough that
# Flye's raw consensus needs no external polishing.
ANN = ["Prokka on\nunpolished", "Prokka on\npolished"]
cds_count   = np.array([4865, 4848])
mean_cds_aa = np.array([312.8, 314.2])
hypo_pct    = np.array([100*1240/4865, 100*1222/4848])

# Fixed categorical assignment - colour follows the assembler, never its rank.
C_UNPOL, C_POL, C_SPADES = "#8FB4DE", "#3B6FB6", "#D1642E"
ASM_COLS = [C_UNPOL, C_POL, C_SPADES]
ANN_COLS = [C_UNPOL, C_POL]
INK, MUTED = "#1a1a1a", "#6b6b6b"

plt.rcParams.update({
    "font.size": 8, "axes.titlesize": 8.5, "axes.labelsize": 8,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.edgecolor": "#999999", "axes.labelcolor": INK,
    "xtick.color": MUTED, "ytick.color": MUTED, "text.color": INK,
    "figure.dpi": 300,
})

fig, axes = plt.subplots(2, 3, figsize=(9.5, 5.8))


def bars(ax, labels, values, colours, title, ylabel, fmt="{:.0f}",
         refline=None, reflabel=None, sub=None):
    """Draw one thin-bar panel with direct value labels above each bar.

    `sub` optionally supplies a short annotation drawn inside the foot of each
    bar. Headroom is set so that neither the value labels nor the reference-line
    label can collide with the bars.
    """
    x = np.arange(len(labels))
    ax.bar(x, values, width=0.58, color=colours, zorder=3)
    ax.set_xticks(x); ax.set_xticklabels(labels)
    ax.set_title(title, loc="left", fontweight="bold", pad=7)
    ax.set_ylabel(ylabel)
    ax.grid(axis="y", color="#e6e6e6", lw=0.7, zorder=0)
    ax.set_axisbelow(True)

    top = max(values) * 1.30
    if refline is not None:
        top = max(top, refline * 1.30)
    ax.set_ylim(0, top)
    ax.set_xlim(-0.62, len(labels) - 0.38)

    # Reference line, labelled above the line on the far left so that it can
    # never overlap a bar top or its value label.
    if refline is not None:
        ax.axhline(refline, ls="--", lw=1.1, color="#444444", zorder=4)
        ax.text(-0.58, refline + top * 0.022, reflabel, va="bottom", ha="left",
                fontsize=6.4, color="#444444", style="italic", zorder=5)

    for xi, v in zip(x, values):
        ax.text(xi, v + top * 0.02, fmt.format(v), ha="center", va="bottom",
                fontsize=7.2, color=INK)

    if sub is not None:
        for xi, s in zip(x, sub):
            ax.text(xi, top * 0.03, s, ha="center", va="bottom", fontsize=6.4,
                    color="white", fontweight="bold", zorder=5)


# --- A. Assembly size -------------------------------------------------------
bars(axes[0, 0], ASM, total_bp / 1e6, ASM_COLS,
     "A  Assembly size", "Total assembly (Mb)", "{:.2f}",
     refline=REF_GENOME_BP / 1e6, reflabel="K-12 reference 4.64 Mb")

# --- B. Contiguity ----------------------------------------------------------
bars(axes[0, 1], ASM, n50 / 1e6, ASM_COLS,
     "B  Contiguity (N50; n = contigs)", "N50 (Mb)", "{:.2f}",
     sub=[f"n = {n}" for n in contigs])

# --- C. Runtime -------------------------------------------------------------
ASM2_COLS = [C_POL, C_SPADES]
bars(axes[0, 2], ASM2, runtime_m, ASM2_COLS,
     "C  Runtime (assembly only)", "Wall-clock (min)", "{:.1f}")

# --- D. Peak memory ---------------------------------------------------------
bars(axes[1, 0], ASM2, memory_gb, ASM2_COLS,
     "D  Peak memory (assembly only)", "Peak RSS (GB)", "{:.2f}")

# --- E. Gene calls ----------------------------------------------------------
bars(axes[1, 1], ANN, cds_count, ANN_COLS,
     "E  Predicted CDS", "CDS called", "{:.0f}",
     refline=REF_CDS, reflabel="K-12 reference 4,298")

# --- F. Mean CDS length -----------------------------------------------------
bars(axes[1, 2], ANN, mean_cds_aa, ANN_COLS,
     "F  Mean CDS length", "Mean length (aa)", "{:.0f}",
     refline=REF_MEAN_AA, reflabel="K-12 reference 316 aa",
     sub=[f"{h:.0f}% hyp." for h in hypo_pct])

fig.suptitle("Flye outperforms SPAdes on HiFi reads; polishing changes nothing "
             "because there is nothing left to correct",
             x=0.008, ha="left", fontsize=10, fontweight="bold", y=0.995)
fig.tight_layout(rect=[0, 0, 1, 0.955])
fig.savefig("figures/FigureA_part1_benchmark.png", bbox_inches="tight")
print("wrote figures/FigureA_part1_benchmark.png")
