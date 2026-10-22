"""Train scVI on the full Xenium 5k panel and save the latent space (robustness check for the CellCharter compartments).

Run in a SEPARATE conda env: in karospace311 the `scvi` module is the unrelated scvi 0.2.2, and scvi-tools
installs a module with the same name.
    conda create -n scvi-env python=3.12 -y && conda activate scvi-env
    pip install scvi-tools anndata h5py
    python scvi_full_panel.py

Output: ../data/X_scVI_5k.npy (cells in the same order as the full object) and the model in ../data/scvi_5k_model/.
"""
import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scvi
from anndata.io import read_elem

FULL = '../data/mtDNA_DSB_5k_clustered_annotation_with_rbd_2_cytetype_brain_novae4.companion.ready.h5ad'

# read only raw counts + sample_id (no log-normalized X, saves ~4 GB of RAM)
with h5py.File(FULL, 'r') as h:
    counts = read_elem(h['layers/counts'])
    obs = pd.DataFrame({'sample_id': read_elem(h['obs/sample_id'])}, index=read_elem(h['obs/_index']))
    var = pd.DataFrame(index=read_elem(h['var']).index)

adata = ad.AnnData(X=counts, obs=obs, var=var)
print(adata)

# same choices as the 250-HVG run where known: batch = sample, 10 latent dimensions
scvi.model.SCVI.setup_anndata(adata, batch_key='sample_id')
model = scvi.model.SCVI(adata, n_latent=10)
model.train(max_epochs=50, batch_size=1024, early_stopping=True, accelerator='gpu')

np.save('../data/X_scVI_5k.npy', model.get_latent_representation())
model.save('../data/scvi_5k_model', overwrite=True)
