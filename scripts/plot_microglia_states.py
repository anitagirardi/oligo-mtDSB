"""Summary figure of the microglial state shift in responders, for the README.

Validated states only (notebook 05): Surveillance and DAM-like (Phagocytosis + Neuroprotection), from the classifier
(primary method). Rows: compartment 3 (globus pallidus + internal capsule) and 4 (thalamus). One dot per mouse = mean
probability over its microglia in the compartment; line = group mean. p = exact permutation responders vs controls,
for the classifier and for Tangram.
Inputs: ../data/microglia_classifier_probs.csv (notebook 05), ../data/microglia_tangram.h5ad (notebook 03).
Run from scripts/:  python plot_microglia_states.py
"""
import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import compartment_stats as cs

MG = '../data/microglia_tangram.h5ad'
CLS = '../data/microglia_classifier_probs.csv'
OUT = '../results/cellcharter/tangram/microglia_states_summary.png'
COMPS = {'3': 'Pallidum /\ninternal capsule (3)', '4': 'Thalamus (4)'}
STATES = {'Surveillance': ['Surveillance'], 'DAM-like': ['Phagocytosis', 'Neuroprotection']}
GROUPS = ['control', 'non-responder', 'responder']
XLABELS = ['Control', 'Non-\nresponder', 'Responder']
MIN_MG = 50

# reference palette tokens (dataviz skill): emphasis form = responders in red, other groups grey
RED, GREY, INK, INK2, SURFACE, RULE = '#e34948', '#9a9893', '#0b0b0b', '#52514e', '#fcfcfb', '#d9d8d4'
COLOURS = {'control': GREY, 'non-responder': GREY, 'responder': RED}

q = ad.read_h5ad(MG)
P_tg = pd.DataFrame(q.obsm['state_prob'], index=q.obs_names, columns=list(q.uns['state_prob_columns']))
P_cls = pd.read_csv(CLS, index_col=0).loc[q.obs_names]
V = {m: pd.DataFrame({k: P[v].sum(1) for k, v in STATES.items()}) for m, P in [('classifier', P_cls), ('Tangram', P_tg)]}
status = q.obs[['sample_id', 'responder_status']].drop_duplicates().set_index('sample_id')['responder_status'].astype(str)


def mouse_means(vals, comp):
    m = (q.obs['compartment'].astype(str) == comp).values
    g = vals[m].groupby(q.obs['sample_id'][m].astype(str).values)
    return g.mean()[g.size() >= MIN_MG]


def p_responders(means, state):
    grp = status.reindex(means.index)
    keep = grp.isin(['control', 'responder'])
    return cs.stratified_permutation_test(means.loc[keep, state], grp[keep], pd.Series('all', index=grp[keep].index),
                                          treat='responder', ref='control')


rng = np.random.default_rng(0)
fig, axes = plt.subplots(len(COMPS), len(STATES), figsize=(7.4, 6.6), facecolor=SURFACE)
for r, comp in enumerate(COMPS):
    means = {m: mouse_means(v, comp) for m, v in V.items()}
    for c, st in enumerate(STATES):
        ax = axes[r, c]
        ax.set_facecolor(SURFACE)
        vals = means['classifier'][st] * 100
        grp = status.reindex(vals.index)
        for x, g in enumerate(GROUPS):
            v = vals[grp == g]
            ax.scatter(x + rng.uniform(-0.12, 0.12, len(v)), v, s=46, color=COLOURS[g], edgecolor=SURFACE,
                       linewidth=1.2, zorder=3)
            ax.plot([x - 0.25, x + 0.25], [v.mean()] * 2, color=INK, linewidth=2, zorder=4, solid_capstyle='round')
        t_cls, t_tg = p_responders(means['classifier'], st), p_responders(means['Tangram'], st)
        arrow = '↑' if t_cls['effect'] > 0 else '↓'
        if np.isclose(t_cls['p_perm'], t_tg['p_perm']):
            note = f"p = {t_cls['p_perm']:.3f}, classifier and Tangram"
        else:
            note = f"p = {t_cls['p_perm']:.3f} (classifier), {t_tg['p_perm']:.3f} (Tangram)"
        ax.set_title(f'{st} {arrow}', fontsize=11, color=INK, loc='left', pad=18)
        ax.text(0, 1.02, note, transform=ax.transAxes, fontsize=7.5, color=INK2, va='bottom')
        ax.set_xticks(range(len(GROUPS)))
        ax.set_xticklabels(XLABELS if r == len(COMPS) - 1 else [], fontsize=8.5, color=INK2)
        ax.set_xlim(-0.5, len(GROUPS) - 0.5)
        ax.tick_params(axis='y', labelsize=8, colors=INK2, length=0)
        ax.tick_params(axis='x', length=0)
        ax.grid(axis='y', color=RULE, linewidth=0.6)
        ax.set_axisbelow(True)
        for s in ['top', 'right', 'left']:
            ax.spines[s].set_visible(False)
        ax.spines['bottom'].set_color(RULE)
    axes[r, 0].set_ylabel(f'{COMPS[comp]}\n\nmean probability (%)', fontsize=9, color=INK)

fig.suptitle('Microglia of responders leave surveillance for a DAM-like state', x=0.01, ha='left',
             fontsize=12, color=INK, y=1.02)
fig.text(0.01, -0.02, 'One dot per mouse (mean classifier probability over its microglia in the compartment); line = '
         'group mean.\nStates from the Prinz myeloid atlas, validated on Xenium (notebook 05); DAM-like = Phagocytosis '
         '+ Neuroprotection.\nProbabilities are soft, so only differences between groups are meaningful. p = responders vs '
         'controls,\nexact permutation (smallest possible 0.012).', fontsize=7.5, color=INK2, va='top')
fig.tight_layout(w_pad=2.5, h_pad=3.0)
fig.savefig(OUT, dpi=200, bbox_inches='tight', facecolor=SURFACE)
print('saved', OUT)
