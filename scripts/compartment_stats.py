"""Helpers for per-section analyses of the CellCharter compartments.

The mouse (= section) is the replicate in every test: all functions first summarise cells into one value
per section, then compare control vs mtDSB across sections.
"""
from itertools import combinations, permutations, product

import anndata as ad
import cellcharter as cc
import numpy as np
import pandas as pd
import scipy.sparse as sp
import statsmodels.api as sm
from scipy import stats
from statsmodels.stats.multitest import multipletests


# ---------------------------------------------------------------- per-section tables
def sample_table(obs, sample_key, cols):
    """One row per sample with sample-level columns (checks they are constant within a sample)."""
    tab = obs[[sample_key] + cols].drop_duplicates()
    if tab[sample_key].duplicated().any():
        raise ValueError(f'{cols} are not constant within {sample_key}')
    return tab.set_index(sample_key)


def fraction_table(obs, sample_key, group_key, label_key, labels, min_cells=100):
    """Cells of each `label` among the cells of each (sample, group).

    Returns a long table with k (cells of the label), n (cells in the group) and frac.
    frac is NaN when n < min_cells (absent or too small, not a real zero).
    """
    n = obs.groupby([sample_key, group_key], observed=True).size()
    rows = []
    for lab in labels:
        k = (obs[obs[label_key] == lab]
             .groupby([sample_key, group_key], observed=True).size()
             .reindex(n.index, fill_value=0))
        df = pd.DataFrame({'k': k, 'n': n}).reset_index()
        df['label'] = lab
        rows.append(df)
    out = pd.concat(rows, ignore_index=True)
    out['frac'] = (out['k'] / out['n']).where(out['n'] >= min_cells)
    return out


def logit(k, n, eps=0.5):
    """Empirical logit with a pseudocount, so that zero counts stay finite."""
    k, n = np.asarray(k, float), np.asarray(n, float)
    return np.log((k + eps) / (n - k + eps))


# ---------------------------------------------------------------- tests
def stratified_permutation_test(y, cond, strata, treat='mtDSB', ref='control'):
    """Exact permutation test, labels permuted only within strata (e.g. DOX timepoint).

    Statistic: mean over strata of (mean treat - mean ref), i.e. the blocked difference.
    y, cond, strata: Series indexed by sample. Samples with NaN y are dropped.
    With 3 vs 3 in each of 2 strata there are 20 * 20 = 400 relabelings (min two-sided p = 2/400).
    """
    ok = y.notna()
    y, cond, strata = y[ok], cond[ok], strata[ok]

    def stat(is_treat):
        diffs = [y[is_treat & (strata == s)].mean() - y[~is_treat & (strata == s)].mean()
                 for s in strata.unique()
                 if (is_treat & (strata == s)).any() and (~is_treat & (strata == s)).any()]
        return np.mean(diffs) if diffs else np.nan

    observed = stat(cond == treat)
    blocks = []
    for s in strata.unique():
        idx = strata.index[strata == s]
        n_treat = int((cond[idx] == treat).sum())
        blocks.append([set(c) for c in combinations(idx, n_treat)])
    null = np.array([stat(pd.Series(y.index.isin(set().union(*combo)), index=y.index))
                     for combo in product(*blocks)])
    p = np.mean(np.abs(null) >= abs(observed) - 1e-12)
    return {'effect': observed, 'p_perm': p, 'n_perm': len(null),
            'n_treat': int((cond == treat).sum()), 'n_ref': int((cond == ref).sum())}


def quasibinomial_test(k, n, cond, covariate=None, treat='mtDSB'):
    """Quasi-binomial GLM, one row per sample: (k, n-k) ~ condition [+ covariate].

    Returns the log odds ratio for treat and a p-value from the t distribution on the residual df
    (dispersion estimated by Pearson chi2, so between-mouse overdispersion is accounted for).
    """
    df = pd.DataFrame({'k': k, 'n': n, 'treat': (cond == treat).astype(float)})
    if covariate is not None:
        df = df.join(pd.get_dummies(covariate, prefix='cov', drop_first=True, dtype=float))
    df = df.dropna()
    X = sm.add_constant(df.drop(columns=['k', 'n']))
    # proportions + var_weights=n: with a 2-column (k, n-k) endog, scale='X2' gives a wrong dispersion in statsmodels
    res = sm.GLM(df['k'] / df['n'], X, family=sm.families.Binomial(),
                 var_weights=df['n']).fit(scale='X2')
    t = res.params['treat'] / res.bse['treat']
    return {'logOR': res.params['treat'], 'p_qb': 2 * stats.t.sf(abs(t), res.df_resid)}


def test_fractions(long, sample_info, sample_key, by, cond_key='condition', strata_key='timepoint',
                   treat='mtDSB', ref='control', min_samples=3):
    """Run both tests for every group in `by` (columns of `long`, e.g. ['label', 'compartment']).

    long: output of fraction_table (k, n, frac per sample). The permutation test works on logit fractions;
    the quasi-binomial GLM on counts (pairs with n < min_cells are dropped, same as frac NaN).
    Samples whose cond_key is neither treat nor ref (e.g. non-responders) are left out.
    """
    rows = []
    for key, sub in long.groupby(by, observed=True):
        sub = sub.set_index(sample_key)
        sub = sub[sub['frac'].notna() & sample_info.loc[sub.index, cond_key].isin([treat, ref]).values]
        info = sample_info.loc[sub.index]
        if ((info[cond_key] == treat).sum() < min_samples) or ((info[cond_key] == ref).sum() < min_samples):
            continue
        y = pd.Series(logit(sub['k'], sub['n']), index=sub.index)
        perm = stratified_permutation_test(y, info[cond_key], info[strata_key], treat, ref)
        qb = quasibinomial_test(sub['k'], sub['n'], info[cond_key], info[strata_key], treat)
        rec = dict(zip(by, key if isinstance(key, tuple) else (key,)))
        rec.update({f'mean_{ref}': sub.loc[info[cond_key] == ref, 'frac'].mean(),
                    f'mean_{treat}': sub.loc[info[cond_key] == treat, 'frac'].mean()})
        rec.update(perm)
        rec.update(qb)
        rows.append(rec)
    return pd.DataFrame(rows)


def add_fdr(res, cols=('p_perm', 'p_qb'), family=None):
    """BH-FDR per p-value column, optionally within families (e.g. per cell type)."""
    res = res.copy()
    groups = res.groupby(family) if family else [(None, res)]
    for col in cols:
        res[col.replace('p_', 'padj_')] = np.nan
        for _, sub in groups:
            ok = sub[col].notna()
            if ok.any():
                res.loc[sub.index[ok], col.replace('p_', 'padj_')] = multipletests(sub.loc[ok, col], method='fdr_bh')[1]
    return res


# ---------------------------------------------------------------- neighborhood enrichment per section
def per_section_nhood(obs, connectivities, sample_key, cluster_key, **nhood_kwargs):
    """cc.gr.nhood_enrichment run separately in each section (sections are disconnected in the graph).

    Returns {sample: enrichment DataFrame}. Categories absent from a section give NaN rows/columns.
    """
    light = ad.AnnData(obs=obs[[sample_key, cluster_key]].copy(),
                       obsp={'cc_connectivities': connectivities})
    out = {}
    for s in light.obs[sample_key].unique():
        sub = light[light.obs[sample_key] == s].copy()
        enr = cc.gr.nhood_enrichment(sub, cluster_key=cluster_key, connectivity_key='cc',
                                     copy=True, **nhood_kwargs)['enrichment']
        present = sub.obs[cluster_key].value_counts() > 0
        absent = present.index[~present.values]
        enr.loc[absent, :] = np.nan
        enr.loc[:, absent] = np.nan
        out[s] = enr
    return out


# ---------------------------------------------------------------- pseudobulk
def pseudobulk(adata, groupby, layer='counts', min_cells=50):
    """Sum raw counts over cells for each combination of `groupby` columns.

    Returns (counts DataFrame groups x genes, metadata DataFrame with groupby columns and n_cells).
    Groups with fewer than min_cells cells are dropped.
    """
    keys = adata.obs[groupby].astype(str).agg('|'.join, axis=1)
    codes, uniques = pd.factorize(keys)
    G = sp.csr_matrix((np.ones(len(codes)), (codes, np.arange(len(codes)))), shape=(len(uniques), len(codes)))
    X = adata.layers[layer] if layer else adata.X
    summed = np.asarray((G @ X).todense()) if sp.issparse(X) else G @ X
    meta = pd.DataFrame([u.split('|') for u in uniques], columns=groupby, index=uniques)
    meta['n_cells'] = np.bincount(codes)
    counts = pd.DataFrame(np.rint(summed).astype(int), index=uniques, columns=adata.var_names)
    keep = meta['n_cells'] >= min_cells
    return counts[keep], meta[keep]


def run_pydeseq2(counts, meta, design='~timepoint + condition', contrast=('condition', 'mtDSB', 'control'),
                 min_count=10, min_samples=3):
    """PyDESeq2 on a pseudobulk (one row per mouse). Genes need >= min_count counts in >= min_samples samples."""
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats

    counts = counts.loc[:, (counts >= min_count).sum(0) >= min_samples]
    dds = DeseqDataSet(counts=counts, metadata=meta, design=design, quiet=True)
    dds.deseq2()
    res = DeseqStats(dds, contrast=list(contrast), quiet=True)
    res.summary()
    return res.results_df.sort_values('pvalue')


def clr(X, eps=1e-6):
    """Centred log-ratio transform of compositions (rows sum to 1, e.g. mean state probabilities)."""
    L = np.log(np.asarray(X, float) + eps)
    return L - L.mean(1, keepdims=True)


def permanova(X, labels, strata):
    """PERMANOVA (Euclidean distance) with an exact permutation test: labels are permuted only within strata.

    The multivariate analogue of a one-way ANOVA with permutation p-values, so it needs no normality assumption.
    X: DataFrame samples x features (e.g. CLR of mean state probabilities); labels, strata: Series indexed like X.
    Works for 2 or more groups; all distinct relabelings within strata are enumerated.
    """
    X = X.loc[labels.index].to_numpy(float)
    lab = labels.astype(str).to_numpy()
    strata = strata.loc[labels.index].astype(str).to_numpy()
    n, groups = len(lab), np.unique(lab)
    if len(groups) < 2 or min((lab == g).sum() for g in groups) < 2:
        return {'pseudo_F': np.nan, 'p_perm': np.nan, 'n_perm': 0}      # too few mice per group

    def pseudo_f(lb):
        ss_w = sum(((X[lb == g] - X[lb == g].mean(0)) ** 2).sum() for g in groups)
        ss_t = ((X - X.mean(0)) ** 2).sum()
        return ((ss_t - ss_w) / (len(groups) - 1)) / (ss_w / (n - len(groups)))

    observed = pseudo_f(lab)
    blocks = []
    for s in np.unique(strata):
        idx = np.where(strata == s)[0]
        blocks.append([(idx, np.array(p)) for p in set(permutations(lab[idx]))])
    null = []
    for combo in product(*blocks):
        lb = lab.copy()
        for idx, p in combo:
            lb[idx] = p
        null.append(pseudo_f(lb))
    null = np.array(null)
    return {'pseudo_F': observed, 'p_perm': np.mean(null >= observed - 1e-12), 'n_perm': len(null)}


# ---------------------------------------------------------------- plotting
GROUP_COLOURS = {'control': '#3182bd', 'non-responder': '#fdae6b', 'responder': '#d62728'}
TP_MARKERS = {'P21': 'o', 'P60': '^'}


def strip_by_mouse(d, x, y, group_key='responder_status', tp_key='timepoint', col=None, order=None,
                   sharey=True, height=3.2, aspect=1.2, title=None, fname=None):
    """One dot per mouse: colour = control / non-responder / responder, marker = DOX timepoint (P21 circle, P60 triangle).

    d: long table with one row per mouse (and per x / col category) holding x, y, group_key and tp_key.
    """
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    d = d.copy()
    d[x] = d[x].astype(str)
    xs = [str(v) for v in order] if order is not None else list(dict.fromkeys(d[x]))
    panels = list(dict.fromkeys(d[col])) if col else [None]
    fig, axes = plt.subplots(1, len(panels), figsize=(height * aspect * len(panels), height),
                             sharey=sharey, squeeze=False)
    rng = np.random.default_rng(0)
    offsets = dict(zip(GROUP_COLOURS, [-0.25, 0.0, 0.25]))
    for ax, p in zip(axes[0], panels):
        sub = d if p is None else d[d[col] == p]
        for i, xv in enumerate(xs):
            for g, colour in GROUP_COLOURS.items():
                for tp, marker in TP_MARKERS.items():
                    v = sub.loc[(sub[x] == xv) & (sub[group_key] == g) & (sub[tp_key] == tp), y].dropna()
                    ax.scatter(i + offsets[g] + rng.uniform(-0.06, 0.06, len(v)), v, color=colour, marker=marker,
                               s=36, edgecolor='white', linewidth=0.5, zorder=3)
        ax.set_xticks(range(len(xs)))
        ax.set_xticklabels(xs)
        ax.set_xlabel(x)
        if p is not None:
            ax.set_title(f'{col} = {p}', fontsize=10)
    axes[0, 0].set_ylabel(y)
    handles = [Line2D([], [], marker='o', linestyle='', color=c, label=g) for g, c in GROUP_COLOURS.items()]
    handles += [Line2D([], [], marker=m, linestyle='', color='grey', label=t) for t, m in TP_MARKERS.items()]
    fig.legend(handles=handles, loc='center left', bbox_to_anchor=(1.0, 0.5), frameon=False)
    if title:
        fig.suptitle(title, y=1.03)
    fig.tight_layout()
    if fname:
        fig.savefig(fname, dpi=150, bbox_inches='tight')
    return fig
