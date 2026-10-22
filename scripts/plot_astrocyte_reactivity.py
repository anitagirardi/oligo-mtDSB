"""Astrocyte reactivity in compartment 3 vs the corpus callosum (compartment 10), for the README.

Per mouse (pseudobulk of its astrocytes in the compartment, min 50 cells):
- four reactive-astrocyte genes, log2 CPM; panel note = DESeq2 log2FC and padj, responders vs controls (notebook 02);
- two programme scores (Hallmark IL6/JAK-STAT3; interferon alpha + gamma): mean over the programme's genes of the
  log2 CPM expressed in SD units of the controls (0 = control mean); p = exact permutation, responders vs controls.
Inputs: ../data/adata_5k_cellcharter.h5ad (only the astrocyte rows are read), ../data/genesets/, ../results/cellcharter/.
Run from scripts/:  python plot_astrocyte_reactivity.py
"""
import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from anndata.io import read_elem, sparse_dataset

import compartment_stats as cs

DATA = '../data/adata_5k_cellcharter.h5ad'
HALLMARK = '../data/genesets/mh.all.v2024.1.Mm.symbols.gmt'
DE = '../results/cellcharter/pseudobulk_DE_comp{}_Astrocytes_responder.csv'
OUT = '../results/cellcharter/astrocyte_reactivity_comp3_vs_10.png'
COMPS = {'3': 'Pallidum /\ninternal capsule (3)', '10': 'Corpus\ncallosum (10)'}
GENES = ['Gfap', 'Serpina3n', 'Osmr', 'Stat3']
PROGRAMMES = {'IL6/JAK-STAT3': ['HALLMARK_IL6_JAK_STAT3_SIGNALING'],
              'Interferon (α + γ)': ['HALLMARK_INTERFERON_ALPHA_RESPONSE', 'HALLMARK_INTERFERON_GAMMA_RESPONSE']}
GROUPS = ['control', 'non-responder', 'responder']
XLABELS = ['Control', 'Non-\nresp.', 'Resp.']
MIN_CELLS, MIN_CPM = 50, 10

# reference palette tokens (dataviz skill): emphasis form = responders in red, other groups grey
RED, GREY, INK, INK2, SURFACE, RULE = '#e34948', '#9a9893', '#0b0b0b', '#52514e', '#fcfcfb', '#d9d8d4'
COLOURS = {'control': GREY, 'non-responder': GREY, 'responder': RED}

# --- per-mouse pseudobulk of astrocytes in compartments 3 and 10 (reads only those rows)
with h5py.File(DATA, 'r') as h:
    obs = pd.DataFrame({c: np.asarray(read_elem(h['obs'][c])).astype(str)
                        for c in ['sample_id', 'responder_status', 'compartment', 'cell_class_work']})
    genes = pd.Index(read_elem(h['var']).index)
    idx = np.flatnonzero((obs['cell_class_work'] == 'Astrocytes') & obs['compartment'].isin(list(COMPS)))
    counts = sparse_dataset(h['layers/counts'])[idx]
cells = obs.iloc[idx]
status = obs.drop_duplicates('sample_id').set_index('sample_id')['responder_status']

hall = {l.split('\t')[0]: l.rstrip('\n').split('\t')[2:] for l in open(HALLMARK)}
data = {}
for comp in COMPS:
    m = (cells['compartment'] == comp).to_numpy()
    df = pd.DataFrame(counts[m].toarray(), columns=genes).groupby(cells['sample_id'].to_numpy()[m])
    n = df.size()
    pb = df.sum().loc[n[n >= MIN_CELLS].index]
    cpm = pb.div(pb.sum(1), axis=0) * 1e6
    logcpm = np.log2(cpm + 1)
    ctrl = status.reindex(logcpm.index) == 'control'
    vals = logcpm[GENES].copy()
    for name, sets in PROGRAMMES.items():
        g = sorted({x for s in sets for x in hall[s]} & set(genes))
        g = [x for x in g if cpm[x].mean() >= MIN_CPM]                 # expressed genes only
        z = (logcpm[g] - logcpm.loc[ctrl, g].mean()) / (logcpm.loc[ctrl, g].std() + 1e-6)
        vals[name] = z.mean(1)
    data[comp] = vals

# --- figure: rows = compartments, columns = genes + programmes
cols = GENES + list(PROGRAMMES)
rng = np.random.default_rng(0)
fig, axes = plt.subplots(len(COMPS), len(cols), figsize=(14, 6.6), sharey='col', facecolor=SURFACE)
for r, comp in enumerate(COMPS):
    de = pd.read_csv(DE.format(comp), index_col=0)
    vals = data[comp]
    grp = status.reindex(vals.index)
    for c, col in enumerate(cols):
        ax = axes[r, c]
        ax.set_facecolor(SURFACE)
        for x, g in enumerate(GROUPS):
            v = vals.loc[grp == g, col]
            ax.scatter(x + rng.uniform(-0.12, 0.12, len(v)), v, s=40, color=COLOURS[g], edgecolor=SURFACE,
                       linewidth=1.1, zorder=3)
            ax.plot([x - 0.25, x + 0.25], [v.mean()] * 2, color=INK, linewidth=2, zorder=4, solid_capstyle='round')
        if col in GENES:
            fc, q = de.loc[col, 'log2FoldChange'], de.loc[col, 'padj']
            note = f'log2FC {fc:+.1f}, ' + ('padj < 0.001' if q < 0.001 else f'padj {q:.2f}')
        else:
            keep = grp.isin(['control', 'responder'])
            t = cs.stratified_permutation_test(vals.loc[keep, col], grp[keep], pd.Series('all', index=grp[keep].index),
                                               treat='responder', ref='control')
            note = f"p = {t['p_perm']:.3f}"
        ax.text(0, 1.02, note, transform=ax.transAxes, fontsize=7.5, color=INK2, va='bottom')
        if r == 0:
            ax.set_title(col, fontsize=11, color=INK, loc='left', pad=18)
        ax.set_xticks(range(len(GROUPS)))
        ax.set_xticklabels(XLABELS if r == len(COMPS) - 1 else [], fontsize=8, color=INK2)
        ax.set_xlim(-0.5, len(GROUPS) - 0.5)
        ax.tick_params(axis='y', labelsize=7.5, colors=INK2, length=0)
        ax.tick_params(axis='x', length=0)
        ax.grid(axis='y', color=RULE, linewidth=0.6)
        ax.set_axisbelow(True)
        for s in ['top', 'right', 'left']:
            ax.spines[s].set_visible(False)
        ax.spines['bottom'].set_color(RULE)
    axes[r, 0].set_ylabel(f'{COMPS[comp]}\n\nlog2 CPM', fontsize=9, color=INK)
    axes[r, len(GENES)].set_ylabel('score (SD from controls)', fontsize=8.5, color=INK2)

fig.suptitle('Astrocytes: strongly reactive in the pallidum / internal capsule, weaker in the corpus callosum',
             x=0.01, ha='left', fontsize=12, color=INK, y=1.02)
fig.text(0.01, -0.03, 'One dot per mouse (pseudobulk of its astrocytes in the compartment); line = group mean. Genes: '
         'DESeq2 responders vs controls. Programme scores: mean of the\nprogramme\'s expressed genes in SD units of the '
         'controls; p = exact permutation responders vs controls (smallest possible 0.012).',
         fontsize=7.5, color=INK2, va='top')
fig.tight_layout(w_pad=1.5, h_pad=3.0)
fig.savefig(OUT, dpi=200, bbox_inches='tight', facecolor=SURFACE)
print('saved', OUT)
