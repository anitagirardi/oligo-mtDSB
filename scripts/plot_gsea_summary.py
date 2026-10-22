"""Compact GSEA summary figures for the README (dot matrices of a few selected pathways).

Colour = NES, size = -log10 FDR, filled = FDR < 0.05, number = NES.
Input: ../results/cellcharter/enrichment/gsea_prerank_all.csv (notebook 02b).
Run from scripts/:  python plot_gsea_summary.py
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.lines import Line2D

IN = '../results/cellcharter/enrichment/gsea_prerank_all.csv'
OUTDIR = '../results/cellcharter/enrichment'
FDR = 0.05

# reference palette tokens (dataviz skill): diverging blue <-> red with a neutral grey midpoint
BLUE, GREY, RED = '#2a78d6', '#f0efec', '#e34948'
INK, INK2, SURFACE, RULE = '#0b0b0b', '#52514e', '#fcfcfb', '#d9d8d4'
CMAP = LinearSegmentedColormap.from_list('div', [BLUE, GREY, RED])
NORM = TwoSlopeNorm(vmin=-2.5, vcenter=0, vmax=2.5)


def size(q):
    # -log10 FDR, capped at 3 (FDR reported as 0 by the permutation test -> cap)
    return 40 + 70 * np.clip(-np.log10(np.maximum(q, 1e-3)), 0, 3)


def dot_matrix(gsea, group, terms, scopes, title, out, n_up=None, note=''):
    """terms: {term: label}; scopes: {compartment: column label}; columns = scopes x (responder, mtDSB).

    n_up: number of 'up' rows before the 'down' block (None = no down block).
    """
    columns = [(cmp, c) for cmp in ['responder', 'mtDSB'] for c in scopes]
    g = gsea[(gsea['group'] == group) & gsea['Term'].isin(terms)]
    nes = g.pivot_table(index='Term', columns=['comparison', 'compartment'], values='NES')
    fdr = g.pivot_table(index='Term', columns=['comparison', 'compartment'], values='FDR q-val')
    nes, fdr = nes.reindex(index=list(terms), columns=columns), fdr.reindex(index=list(terms), columns=columns)

    k = len(scopes)
    step, gap = 1.55, 0.45
    xs = [j * step + (gap if j >= k else 0) for j in range(2 * k)]
    rows = [i if (n_up is None or i < n_up) else i + 0.6 for i in range(len(terms))]
    fig, ax = plt.subplots(figsize=(3.2 + 1.05 * len(xs), 1.2 + 0.47 * len(terms)), facecolor=SURFACE)
    ax.set_facecolor(SURFACE)
    for j, col in enumerate(columns):
        for i, term in enumerate(terms):
            v, q = nes.loc[term, col], fdr.loc[term, col]
            if np.isnan(v):
                ax.text(xs[j], rows[i], 'n.t.', ha='center', va='center', fontsize=7.5, color=INK2)
                continue
            c = CMAP(NORM(v))
            sig = q < FDR
            ax.scatter(xs[j], rows[i], s=size(q), facecolor=c if sig else SURFACE,
                       edgecolor=SURFACE if sig else c, linewidth=1.5 if sig else 2, zorder=3)
            ax.text(xs[j] + 0.27, rows[i], f'{v:+.1f}', va='center', ha='left', fontsize=7.5, color=INK2)

    ax.set_yticks(rows)
    ax.set_yticklabels(list(terms.values()), fontsize=9.5, color=INK)
    ax.set_xticks(xs)
    ax.set_xticklabels(list(scopes.values()) * 2, fontsize=8.5, color=INK2)
    ax.xaxis.tick_top()
    for x0, lab in zip([np.mean(xs[:k]), np.mean(xs[k:])],
                       [f'Responders vs controls\n(3 vs 6 mice)', f'All mtDSB vs controls\n(6 vs 6 mice)']):
        ax.text(x0, -1.9, lab, ha='center', va='bottom', fontsize=9.5, color=INK, fontweight='bold')
    ax.set_xlim(-0.5, xs[-1] + 0.7)
    ax.set_ylim(rows[-1] + 0.7, -0.7)
    ax.axvline((xs[k - 1] + xs[k]) / 2, color=RULE, linewidth=0.8)
    if n_up is not None:
        ax.axhline(n_up - 0.2, color=RULE, linewidth=0.8)
        for y, lab in [(-0.55, 'up'), (n_up + 0.05, 'down')]:
            ax.text(-0.45, y, lab, fontsize=8, color=INK2, style='italic', va='center')
    for s in ax.spines.values():
        s.set_visible(False)
    ax.tick_params(length=0)

    fig.subplots_adjust(left=0.30, right=0.80, top=0.80, bottom=0.04)
    cax = fig.add_axes([0.85, 0.52, 0.015, 0.25])
    cb = fig.colorbar(plt.cm.ScalarMappable(norm=NORM, cmap=CMAP), cax=cax)
    cb.set_label('NES', fontsize=8.5, color=INK2)
    cb.ax.tick_params(labelsize=7.5, colors=INK2)
    cb.outline.set_visible(False)
    handles = [Line2D([], [], marker='o', linestyle='', markersize=np.sqrt(size(q)), markerfacecolor=INK2,
                      markeredgecolor=SURFACE, label=f'FDR {q:g}') for q in [0.05, 0.01, 0.001]]
    handles.append(Line2D([], [], marker='o', linestyle='', markersize=np.sqrt(size(0.2)), markerfacecolor=SURFACE,
                          markeredgecolor=INK2, markeredgewidth=2, label=f'FDR ≥ {FDR} (hollow)'))
    fig.legend(handles=handles, loc='lower left', bbox_to_anchor=(0.83, 0.06), frameon=False, fontsize=8,
               labelcolor=INK2, labelspacing=1.1)
    fig.suptitle(title, x=0.02, ha='left', fontsize=12, color=INK, y=1.02)
    fig.text(0.02, -0.03, 'Pseudobulk DESeq2 per mouse, preranked GSEA; numbers = NES. '
             'FDRs are optimistic (genes permuted, not mice).' + (f'\n{note}' if note else ''),
             fontsize=7.5, color=INK2, va='top')
    fig.savefig(f'{OUTDIR}/{out}', dpi=200, bbox_inches='tight', facecolor=SURFACE)
    print('saved', out)


gsea = pd.read_csv(IN, dtype={'compartment': str})

dot_matrix(
    gsea, 'mature',
    terms={
        'GOBP_CELLULAR_RESPONSE_TO_IONIZING_RADIATION': 'Response to ionizing radiation (GO)',
        'HALLMARK_P53_PATHWAY': 'p53 pathway',
        'HALLMARK_APOPTOSIS': 'Apoptosis',
        'HALLMARK_UNFOLDED_PROTEIN_RESPONSE': 'Unfolded protein response',
        'HALLMARK_MYC_TARGETS_V1': 'MYC targets',
        'HALLMARK_E2F_TARGETS': 'E2F targets',
        'HALLMARK_HYPOXIA': 'Hypoxia',
        'GOBP_AXON_ENSHEATHMENT_IN_CENTRAL_NERVOUS_SYSTEM': 'Axon ensheathment in CNS (GO)',
        'GOBP_OLIGODENDROCYTE_DIFFERENTIATION': 'Oligodendrocyte differentiation (GO)',
    },
    scopes={'3': 'Pallidum /\nint. capsule (3)', '10': 'Corpus\ncallosum (10)'},
    title='Mature oligodendrocytes: pathway enrichment (GSEA)', out='ol_gsea_summary.png', n_up=7)

dot_matrix(
    gsea, 'Microglia',
    terms={
        'GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_PEPTIDE_ANTIGEN': 'Antigen presentation (GO)',
        'GOBP_VACUOLAR_ACIDIFICATION': 'Lysosomal acidification (GO)',
        'GOBP_CELLULAR_RESPONSE_TO_LIPOPROTEIN_PARTICLE_STIMULUS': 'Lipoprotein / lipid handling (GO)',
        'GOBP_POSITIVE_REGULATION_OF_RESPONSE_TO_CYTOKINE_STIMULUS': 'Response to cytokines (GO)',
        'HALLMARK_INFLAMMATORY_RESPONSE': 'Inflammatory response',
        'HALLMARK_TNFA_SIGNALING_VIA_NFKB': 'TNFα signalling via NF-κB',
        'HALLMARK_CHOLESTEROL_HOMEOSTASIS': 'Cholesterol homeostasis',
    },
    scopes={'all': 'Whole\nsection', '3': 'Pallidum /\nint. caps. (3)', '4': 'Thalamus\n(4)', '10': 'Corpus\ncallosum (10)'},
    title='Microglia: pathway enrichment (GSEA)', out='microglia_gsea_summary.png',
    note='One representative term per theme; down-regulated terms in microglia were spill-over from neighbouring '
         'OLs / neurons and are not shown.')
