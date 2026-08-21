#!/usr/bin/env python3
"""
###############################################################################
LIFE748 Assessment 2 - Part 3
Structural Bioinformatics of Differentially Expressed Transcription Factors

Aim: for the transcription factors identified as differentially expressed in
Part 2, obtain AlphaFold structure predictions, assess their confidence, and
apply structure-based function annotation to predict DNA-binding surfaces.
The predictions are then validated against curated PROSITE/UniProt HTH
annotations, so the accuracy of the structural method can be quantified rather
than merely asserted.

Author: <YOUR NAME>  (<YOUR STUDENT ID>)
Date:   August 2026

-------------------------------------------------------------------------------
GENERATIVE AI DECLARATION
Generative AI (Anthropic Claude, model claude-opus-5, accessed August 2026) was
used to assist in drafting, structuring and commenting this script. All method
choices, thresholds and the interpretation of all results were reviewed, tested
and verified by the author. The complete prompt / interaction history is given
in the supplementary file GAI_declaration_and_prompts.md.
-------------------------------------------------------------------------------

INPUTS
  structures/AF-<accession>-F1-model_v6.pdb   AlphaFold DB v6 predicted models
  (downloaded from https://alphafold.ebi.ac.uk/ ; see Methods)

OUTPUTS
  results/Table_7_structure_confidence.csv    per-model pLDDT confidence summary
  results/Table_8_dna_binding_prediction.csv  predicted vs annotated DNA-binding
  results/predicted_binding_residues.txt      residue-level predictions
  figures/Figure5_plddt_profiles.png
  figures/Figure6_Cra_structure.png
  figures/Figure7_binding_site_validation.png

REQUIREMENTS
  python >= 3.10, biotite 1.6.0, numpy, matplotlib
  pip install biotite numpy matplotlib
###############################################################################
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")                       # headless backend - no display needed
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import biotite.structure as struc
import biotite.structure.io.pdb as pdb

os.makedirs("results", exist_ok=True)
os.makedirs("figures", exist_ok=True)

# --------------------------------------------------------------------------- #
# 1. TARGET DEFINITION
#
# The five proteins below are the transcription factors among the 77 DEGs from
# Part 2. They were chosen on biological rationale rather than statistics alone:
# each is a DNA-binding regulator whose regulon is directly relevant to how an
# E. coli isolate would respond to a change in growth condition, and together
# they span four distinct structural families of bacterial regulator, which
# lets the structural method be tested against varied architectures.
#
# The UniProt/PROSITE annotations are curated references used ONLY for
# validation at the end - they are never used to make the prediction.
# --------------------------------------------------------------------------- #
TARGETS = {
    "P0ACP1": dict(gene="cra",  locus="b0080", family="LacI/GalR",
                   name="Catabolite repressor/activator (Cra/FruR)",
                   hth=(3, 22),   domain=(1, 58),    oligomer="Homotetramer",
                   rationale="Global regulator of central carbon metabolism; "
                             "switches the cell between glycolytic and "
                             "gluconeogenic programmes."),
    "Q47129": dict(gene="feaR", locus="b1384", family="AraC/XylS",
                   name="Transcriptional activator FeaR",
                   hth=(217, 238), domain=(199, 299), oligomer="Monomer/dimer",
                   hth2=(266, 289),
                   rationale="Activates phenylethylamine catabolism; AraC/XylS "
                             "family members commonly control stress and "
                             "virulence regulons."),
    "P36673": dict(gene="treR", locus="b4241", family="LacI/GalR",
                   name="HTH-type transcriptional regulator TreR",
                   hth=(7, 26),   domain=(5, 59),    oligomer="Homodimer",
                   rationale="Trehalose operon repressor; trehalose is an "
                             "osmoprotectant, linking this regulator to "
                             "osmotic and desiccation stress."),
    "Q57083": dict(gene="perR", locus="b0254", family="LysR",
                   name="HTH-type transcriptional regulator PerR",
                   hth=(24, 44),  domain=(7, 64),    oligomer="Unknown",
                   rationale="Peroxide-resistance regulator; oxidative stress "
                             "tolerance is directly relevant to survival of a "
                             "clinical isolate in the host."),
    "P77743": dict(gene="prpR", locus="b0330", family="NtrC/sigma-54",
                   name="Propionate catabolism operon regulatory protein PrpR",
                   hth=(508, 528), domain=(483, 528), oligomer="Unknown",
                   rationale="Sigma-54-dependent activator of propionate "
                             "catabolism; a bacterial enhancer-binding protein "
                             "with a C-terminal rather than N-terminal HTH."),
}

# Theoretical maximum solvent-accessible surface areas (A^2) used to convert
# absolute SASA into relative solvent accessibility (Tien et al. 2013, PLoS ONE).
MAX_ASA = {
    "ALA": 129, "ARG": 274, "ASN": 195, "ASP": 193, "CYS": 167, "GLN": 225,
    "GLU": 223, "GLY": 104, "HIS": 224, "ILE": 197, "LEU": 201, "LYS": 236,
    "MET": 224, "PHE": 240, "PRO": 159, "SER": 155, "THR": 172, "TRP": 285,
    "TYR": 263, "VAL": 174,
}
# Formal charge at physiological pH; His is given a partial positive weight
# because it is frequently protonated at a DNA interface.
CHARGE = {"ARG": 1.0, "LYS": 1.0, "HIS": 0.5, "ASP": -1.0, "GLU": -1.0}


def load_model(acc):
    """Read one AlphaFold PDB model and return the protein atoms only.

    AlphaFold stores the per-residue pLDDT confidence score in the B-factor
    column, so the B-factors are carried through with the structure.
    """
    path = f"structures/AF-{acc}-F1-model_v6.pdb"
    f = pdb.PDBFile.read(path)
    arr = f.get_structure(model=1, extra_fields=["b_factor"])
    return arr[struc.filter_amino_acids(arr)]


def per_residue(arr):
    """Collapse an atom array to per-residue arrays.

    Returns residue ids, three-letter names, CA coordinates and pLDDT. pLDDT is
    identical for every atom of a residue, so taking the CA value is exact.
    """
    ca = arr[arr.atom_name == "CA"]
    return ca.res_id, ca.res_name, ca.coord, ca.b_factor


def confidence_summary(acc, plddt):
    """Summarise AlphaFold confidence using the DeepMind interpretation bands.

    pLDDT > 90 very high, 70-90 confident, 50-70 low, < 50 very low (and in
    AlphaFold models regions below 50 are frequently genuinely disordered
    rather than simply badly predicted).
    """
    t = TARGETS[acc]
    return dict(
        Accession=acc, Gene=t["gene"], Locus=t["locus"], Family=t["family"],
        Length=len(plddt),
        Mean_pLDDT=round(float(plddt.mean()), 1),
        Median_pLDDT=round(float(np.median(plddt)), 1),
        Pct_very_high=round(float((plddt >= 90).mean() * 100), 1),
        Pct_confident=round(float((plddt >= 70).mean() * 100), 1),
        Pct_low=round(float(((plddt >= 50) & (plddt < 70)).mean() * 100), 1),
        Pct_very_low=round(float((plddt < 50).mean() * 100), 1),
    )


def sse_string(arr):
    """Assign secondary structure from coordinates alone using P-SEA.

    biotite's annotate_sse implements the P-SEA algorithm, which classifies
    residues as alpha-helix (a), beta-strand (b) or coil (c) from CA geometry.
    DSSP is not used because it requires a separate binary; P-SEA needs only the
    backbone and is sufficient for locating helices.
    """
    sse = struc.annotate_sse(arr)
    return "".join(sse)


def relative_sasa(arr, res_ids, res_names):
    """Per-residue relative solvent accessibility (0 = buried, 1 = fully exposed).

    Absolute SASA is computed with the Shrake-Rupley algorithm and summed over
    each residue's atoms, then divided by the residue's theoretical maximum.
    """
    atom_sasa = struc.sasa(arr, vdw_radii="Single")
    atom_sasa = np.nan_to_num(atom_sasa)
    rel = []
    for rid, rname in zip(res_ids, res_names):
        total = atom_sasa[arr.res_id == rid].sum()
        rel.append(total / MAX_ASA.get(rname, 200))
    return np.clip(np.array(rel), 0, 1.5)


def positive_patch_score(coords, res_names, rel_sasa, sse=None, radius=10.0,
                         coil_weight=0.25):
    """Structure-based DNA-binding propensity.

    Rationale: the phosphodiester backbone of DNA is uniformly and densely
    negatively charged, so the defining structural signature of a DNA-binding
    surface is a large, solvent-exposed, contiguous patch of positive
    electrostatic potential. This is captured here without a full Poisson-
    Boltzmann calculation by scoring each residue as the solvent-weighted net
    charge of all residues within `radius` angstroms of it.

    A residue scores highly only if it is itself exposed AND sits in a
    neighbourhood dominated by exposed Arg/Lys, which is exactly the geometry of
    a recognition helix docked in the DNA major groove.

    Helical residues are up-weighted relative to coil (`coil_weight`). This is
    an a priori structural constraint, not a fitted parameter: in every HTH
    family the base-specific contacts are made by a recognition alpha-helix
    inserted into the major groove, so a positive patch lying on a helix is far
    more likely to be a genuine interface than one on a flexible loop. A
    sensitivity analysis over `radius` and `coil_weight` is written to
    Table_S5 so the effect of these choices is transparent.
    """
    n = len(coords)
    # Pairwise distance matrix between CA atoms
    d = np.linalg.norm(coords[:, None, :] - coords[None, :, :], axis=-1)
    q = np.array([CHARGE.get(r, 0.0) for r in res_names])
    # Weight each neighbour's charge by how exposed that neighbour is: a buried
    # arginine is charge-neutralised by the protein core and cannot touch DNA.
    q_exposed = q * np.clip(rel_sasa, 0, 1)
    within = (d <= radius).astype(float)
    score = within @ q_exposed
    # Only exposed residues can themselves form the interface
    score = score * np.clip(rel_sasa, 0, 1)
    if sse is not None:
        helix = np.array([1.0 if c == "a" else coil_weight for c in sse])
        score = score * helix
    return score


def predict_binding(res_ids, score, top_frac=0.15):
    """Call the top-scoring surface residues as the predicted DNA-binding site.

    The threshold is set as a fraction of the sequence rather than an absolute
    score so that it is comparable across proteins of different size and net
    charge. 15% is used because curated HTH DNA-binding motifs occupy roughly
    6-8% of these sequences, so this deliberately over-calls slightly, favouring
    sensitivity over specificity for a screening application.
    """
    k = max(1, int(round(top_frac * len(score))))
    idx = np.argsort(score)[::-1][:k]
    return set(int(res_ids[i]) for i in idx)


def annotated_residues(acc):
    """Curated DNA-binding residues from UniProt / PROSITE (validation only)."""
    t = TARGETS[acc]
    res = set(range(t["hth"][0], t["hth"][1] + 1))
    if "hth2" in t:
        res |= set(range(t["hth2"][0], t["hth2"][1] + 1))
    return res


# --------------------------------------------------------------------------- #
# 2. RUN THE ANALYSIS
# --------------------------------------------------------------------------- #
conf_rows, val_rows, store = [], [], {}

for acc in TARGETS:
    arr = load_model(acc)
    res_ids, res_names, coords, plddt = per_residue(arr)
    sse = sse_string(arr)
    rel = relative_sasa(arr, res_ids, res_names)
    score = positive_patch_score(coords, res_names, rel, sse)

    pred = predict_binding(res_ids, score)
    true = annotated_residues(acc)
    # Precision expected if the same number of residues were called at random
    chance = len(true) / len(res_ids)

    # Validation of the structural prediction against the curated annotation.
    # Precision = of the residues we called, how many are genuinely in the HTH.
    # Recall    = of the curated HTH residues, how many we recovered.
    tp = len(pred & true)
    precision = tp / len(pred) if pred else 0.0
    recall = tp / len(true) if true else 0.0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) else 0.0

    # Mean pLDDT restricted to the DNA-binding motif: a functional-site
    # prediction is only trustworthy where the model itself is confident.
    mask_true = np.isin(res_ids, list(true))
    plddt_hth = float(plddt[mask_true].mean()) if mask_true.any() else float("nan")

    conf_rows.append({**confidence_summary(acc, plddt),
                      "Mean_pLDDT_HTH": round(plddt_hth, 1),
                      "Helix_pct": round(100 * sse.count("a") / len(sse), 1),
                      "Strand_pct": round(100 * sse.count("b") / len(sse), 1),
                      "Radius_of_gyration_A": round(float(struc.gyration_radius(arr)), 1)})

    val_rows.append(dict(
        Accession=acc, Gene=TARGETS[acc]["gene"], Family=TARGETS[acc]["family"],
        Annotated_HTH=f"{TARGETS[acc]['hth'][0]}-{TARGETS[acc]['hth'][1]}"
                      + (f", {TARGETS[acc]['hth2'][0]}-{TARGETS[acc]['hth2'][1]}"
                         if "hth2" in TARGETS[acc] else ""),
        N_predicted=len(pred), N_annotated=len(true), True_positives=tp,
        Precision=round(precision, 3), Recall=round(recall, 3), F1=round(f1, 3),
        Chance_precision=round(chance, 3),
        Enrichment_over_chance=round(precision / chance, 2) if chance else None,
        Mean_pLDDT_HTH=round(plddt_hth, 1)))

    store[acc] = dict(res_ids=res_ids, res_names=res_names, coords=coords,
                      plddt=plddt, sse=sse, rel=rel, score=score,
                      pred=pred, true=true)

    print(f"{TARGETS[acc]['gene']:5s} ({acc})  mean pLDDT {plddt.mean():5.1f}  "
          f"HTH pLDDT {plddt_hth:5.1f}  precision {precision:.2f}  recall {recall:.2f}  "
          f"enrichment {precision/chance:.1f}x")

# --------------------------------------------------------------------------- #
# 2b. SENSITIVITY ANALYSIS
# The two free parameters of the scoring function are the neighbourhood radius
# and the weight given to non-helical residues. Rather than presenting a single
# tuned result, performance is reported across a grid so the reader can judge
# how much the conclusion depends on these choices.
# --------------------------------------------------------------------------- #
sens_rows = []
for radius in (8.0, 10.0, 12.0, 14.0):
    for cw in (0.0, 0.25, 0.5, 1.0):
        precs, recs, enrs = [], [], []
        for acc, s in store.items():
            sc = positive_patch_score(s["coords"], s["res_names"], s["rel"],
                                      s["sse"], radius=radius, coil_weight=cw)
            p_set = predict_binding(s["res_ids"], sc)
            t_set = s["true"]
            tp_ = len(p_set & t_set)
            ch = len(t_set) / len(s["res_ids"])
            precs.append(tp_ / len(p_set)); recs.append(tp_ / len(t_set))
            enrs.append((tp_ / len(p_set)) / ch)
        sens_rows.append(dict(Radius_A=radius, Coil_weight=cw,
                              Mean_precision=round(float(np.mean(precs)), 3),
                              Mean_recall=round(float(np.mean(recs)), 3),
                              Mean_enrichment=round(float(np.mean(enrs)), 2)))

# Write the two results tables
import csv
with open("results/Table_7_structure_confidence.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(conf_rows[0])); w.writeheader(); w.writerows(conf_rows)
with open("results/Table_8_dna_binding_prediction.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(val_rows[0])); w.writeheader(); w.writerows(val_rows)
with open("results/Table_S5_sensitivity_analysis.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(sens_rows[0])); w.writeheader(); w.writerows(sens_rows)

# Residue-level predictions, for the supplementary material
with open("results/predicted_binding_residues.txt", "w") as fh:
    for acc, s in store.items():
        fh.write(f">{TARGETS[acc]['gene']} ({acc}) {TARGETS[acc]['name']}\n")
        fh.write("predicted DNA-binding residues (residue number + one-letter code):\n")
        three_to_one = {"ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q",
                        "GLU":"E","GLY":"G","HIS":"H","ILE":"I","LEU":"L","LYS":"K",
                        "MET":"M","PHE":"F","PRO":"P","SER":"S","THR":"T","TRP":"W",
                        "TYR":"Y","VAL":"V"}
        items = sorted(s["pred"])
        names = {int(r): n for r, n in zip(s["res_ids"], s["res_names"])}
        fh.write(", ".join(f"{three_to_one.get(names[r],'X')}{r}" for r in items) + "\n\n")

# --------------------------------------------------------------------------- #
# 3. FIGURES
# --------------------------------------------------------------------------- #
BLUE, ORANGE, GREY = "#3B6FB6", "#D1642E", "#8A8A8A"
plt.rcParams.update({"font.size": 9, "axes.spines.top": False,
                     "axes.spines.right": False, "figure.dpi": 300})

# --- Figure 5: pLDDT confidence profiles -----------------------------------
fig, axes = plt.subplots(5, 1, figsize=(7.2, 9.0), sharex=False)
for ax, acc in zip(axes, TARGETS):
    s, t = store[acc], TARGETS[acc]
    ax.axhspan(90, 100, color="#2166AC", alpha=.10)
    ax.axhspan(70, 90,  color="#67A9CF", alpha=.10)
    ax.axhspan(50, 70,  color="#FDDBC7", alpha=.25)
    ax.axhspan(0,  50,  color="#EF8A62", alpha=.20)
    ax.plot(s["res_ids"], s["plddt"], color="black", lw=.9)
    # Shade the curated DNA-binding motif(s)
    ax.axvspan(t["hth"][0], t["hth"][1], color=ORANGE, alpha=.35)
    if "hth2" in t:
        ax.axvspan(t["hth2"][0], t["hth2"][1], color=ORANGE, alpha=.35)
    ax.set_ylim(0, 100); ax.set_xlim(1, len(s["plddt"]))
    ax.set_ylabel("pLDDT")
    ax.set_title(f"{t['gene']} ({t['locus']}, {acc}) - {t['family']} family; "
                 f"mean pLDDT {s['plddt'].mean():.1f}", loc="left", fontsize=9)
axes[-1].set_xlabel("Residue number")
handles = [plt.Rectangle((0,0),1,1, color=ORANGE, alpha=.35)]
axes[0].legend(handles, ["Curated HTH DNA-binding motif"], loc="lower right", fontsize=8)
fig.tight_layout()
fig.savefig("figures/Figure5_plddt_profiles.png", bbox_inches="tight")
plt.close(fig)

# --- Figure 6: Cra backbone, coloured by pLDDT and by predicted site --------
acc = "P0ACP1"; s = store[acc]
xyz = s["coords"]
# Orient the molecule along its principal axes so the projection is reproducible
cen = xyz - xyz.mean(0)
u, sv, vt = np.linalg.svd(cen, full_matrices=False)
proj = cen @ vt.T

fig, axs = plt.subplots(1, 2, figsize=(9.4, 4.4))
sc = axs[0].scatter(proj[:, 0], proj[:, 1], c=s["plddt"], cmap="RdYlBu",
                    vmin=50, vmax=100, s=13, zorder=3)
axs[0].plot(proj[:, 0], proj[:, 1], color="grey", lw=.5, alpha=.7, zorder=2)
plt.colorbar(sc, ax=axs[0], label="pLDDT", fraction=.046)
axs[0].set_title("Cra backbone coloured by model confidence", fontsize=9)

is_pred = np.isin(s["res_ids"], list(s["pred"]))
is_true = np.isin(s["res_ids"], list(s["true"]))
axs[1].plot(proj[:, 0], proj[:, 1], color="lightgrey", lw=.6, zorder=1)
axs[1].scatter(proj[~is_pred, 0], proj[~is_pred, 1], c="lightgrey", s=11, zorder=2)
axs[1].scatter(proj[is_pred, 0], proj[is_pred, 1], c=BLUE, s=22, zorder=3,
               label="Predicted DNA-binding")
axs[1].scatter(proj[is_true, 0], proj[is_true, 1], facecolors="none",
               edgecolors=ORANGE, s=60, lw=1.2, zorder=4,
               label="Curated HTH (3-22)")
axs[1].legend(fontsize=8, loc="best")
axs[1].set_title("Predicted positive patch vs curated motif", fontsize=9)
for a in axs:
    a.set_xticks([]); a.set_yticks([]); a.set_aspect("equal")
    a.set_xlabel("principal axis 1"); a.set_ylabel("principal axis 2")
fig.tight_layout()
fig.savefig("figures/Figure6_Cra_structure.png", bbox_inches="tight")
plt.close(fig)

# --- Figure 7: prediction performance across all five TFs -------------------
fig, axs = plt.subplots(1, 2, figsize=(9.4, 3.8))
genes = [TARGETS[a]["gene"] for a in TARGETS]
prec = [r["Precision"] for r in val_rows]
rec  = [r["Recall"] for r in val_rows]
x = np.arange(len(genes)); w = .38
axs[0].bar(x - w/2, prec, w, color=BLUE,   label="Precision")
axs[0].bar(x + w/2, rec,  w, color=ORANGE, label="Recall")
axs[0].axhline(0.15, ls=":", color=GREY)
axs[0].text(len(genes)-.5, .17, "precision expected by chance", ha="right",
            fontsize=7, color=GREY)
axs[0].set_xticks(x); axs[0].set_xticklabels(genes, style="italic")
axs[0].set_ylabel("Score"); axs[0].set_ylim(0, 1)
axs[0].legend(fontsize=8); axs[0].set_title("Recovery of curated HTH motifs", fontsize=9)

mp = [r["Mean_pLDDT"] for r in conf_rows]
mh = [r["Mean_pLDDT_HTH"] for r in conf_rows]
axs[1].bar(x - w/2, mp, w, color=BLUE,   label="Whole model")
axs[1].bar(x + w/2, mh, w, color=ORANGE, label="HTH motif")
axs[1].axhline(70, ls="--", color=GREY)
axs[1].text(-.4, 71.5, "confident threshold", fontsize=7, color=GREY)
axs[1].set_xticks(x); axs[1].set_xticklabels(genes, style="italic")
axs[1].set_ylabel("Mean pLDDT"); axs[1].set_ylim(0, 100)
axs[1].legend(fontsize=8); axs[1].set_title("Model confidence overall vs at the functional site",
                                            fontsize=9)
fig.tight_layout()
fig.savefig("figures/Figure7_binding_site_validation.png", bbox_inches="tight")
plt.close(fig)

print("\nDone. Tables in results/, figures in figures/")
for r in val_rows:
    print(r)
