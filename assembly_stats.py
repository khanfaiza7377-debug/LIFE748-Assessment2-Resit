#!/usr/bin/env python3
#
# -----------------------------------------------------------------------------
# GENERATIVE AI DECLARATION
# -----------------------------------------------------------------------------
# Generative AI (Anthropic Claude, model claude-opus-5, accessed August 2026) was
# used to assist in drafting and commenting this script. All analytical choices
# and the interpretation of results were reviewed and verified by the author.
# The full prompt / interaction history is provided in
# GAI_declaration_and_prompts.md.
# -----------------------------------------------------------------------------

"""
Compute standard assembly contiguity metrics from one or more FASTA files.

Used in place of QUAST, whose reference download could not be reached from the
analysis environment. The metrics computed here (total length, contig count,
N50, L50, largest contig, GC content, N-rate) are the reference-free subset of
QUAST's report and are calculated identically.

Usage:  python3 assembly_stats.py label=file.fasta [label=file.fasta ...]

GAI declaration: drafted with assistance from Anthropic Claude (claude-opus-5,
August 2026); logic verified by the author. See GAI_declaration_and_prompts.md.
"""
import sys


def read_fasta_lengths(path):
    """Return (list of contig lengths, GC count, ACGT count, N count)."""
    lengths, cur = [], 0
    gc = at = n = 0
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if cur:
                    lengths.append(cur)
                cur = 0
            else:
                s = line.strip().upper()
                cur += len(s)
                gc += s.count("G") + s.count("C")
                at += s.count("A") + s.count("T")
                n += s.count("N")
    if cur:
        lengths.append(cur)
    return lengths, gc, at, n


def nx(lengths, frac):
    """Nx statistic: the length of the contig at which cumulative length first
    reaches `frac` of the assembly, walking from longest to shortest."""
    total = sum(lengths)
    run = 0
    for i, L in enumerate(sorted(lengths, reverse=True), start=1):
        run += L
        if run >= frac * total:
            return L, i
    return 0, 0


print("assembly\tcontigs\ttotal_bp\tlargest_bp\tN50\tL50\tN90\tGC_pct\tNs")
for arg in sys.argv[1:]:
    label, path = arg.split("=", 1)
    L, gc, at, n = read_fasta_lengths(path)
    n50, l50 = nx(L, 0.5)
    n90, _ = nx(L, 0.9)
    print(f"{label}\t{len(L)}\t{sum(L)}\t{max(L)}\t{n50}\t{l50}\t{n90}\t"
          f"{100*gc/(gc+at):.2f}\t{n}")
