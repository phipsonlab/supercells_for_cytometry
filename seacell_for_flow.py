# TODO still work in progress. Use at your own risk!

# coding: utf-8
import pandas as pd
import scanpy as sc
import anndata as ad
import SEACells
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# For covid healthy samples
dat = pd.read_csv("covid_dataset/flow_cell_dat.csv")
dat_healthy = dat[dat['Group'] == 'Healthy']
dat_healthy_asinh = dat_healthy.iloc[:, 17:30]

adata = ad.AnnData(dat_healthy_asinh)

# Assign cell id
adata.obs_names = [f"Cell_{i:d}" for i in range(adata.n_obs)]
dat_healthy_asinh['CellID']=[f"Cell_{i:d}" for i in range(adata.n_obs)]

# SEACell require PCA
sc.tl.pca(adata, n_comps=5)

# SEAcell settings
n_SEACells = 90
build_kernel_on = 'X_pca'
n_waypoint_eigs = 10


model = SEACells.core.SEACells(adata, build_kernel_on=build_kernel_on, n_SEACells=n_SEACells, n_waypoint_eigs=n_waypoint_eigs,convergence_epsilon=1e-5)
model.construct_kernel_matrix()
M = model.kernel_matrix
model.initialize_archetypes()
model.fit(min_iter=10, max_iter=50)

populations_l1 = dat_healthy['Population L1']
adata.obs['CellType_L1'] = populations_l1
sc.tl.umap(adata)


populations_l1 = dat_healthy['Population L1']
adata.obs['celltype'] = populations_l1
adata.obs['celltype'] = pd.DataFrame(populations_l1, index=dat_healthy['CellID'])

dat_healthy['CellID'] = dat_healthy_asinh['CellID']

populations_l1 = pd.DataFrame({'CellID': dat_healthy['CellID'], 'celltype':dat_healthy['Population L1']})
populations_l1 = populations_l1.set_index('CellID')

adata.obs['celltype'] = populations_l1['celltype']

# Just for testing
SEACells.plot.plot_2D(adata, key='X_umap', colour_metacells=True, save_as='SEACELL.pdf')

# Purity plot
SEACell_purity = SEACells.evaluate.compute_celltype_purity(adata,'celltype')
plt.figure(figsize=(4,4))
sns.boxplot(data=SEACell_purity, y='celltype_purity')
plt.title('Celltype Purity')
sns.despine()
plt.savefig("purity.pdf")


sc.pl.umap(adata, color=['celltype','SEACell'], save='umap.pdf')
seacells = adata.obs['SEACell']

dat_seacell = dat_healthy[dat_healthy['CellID'].isin(np.unique(seacells))]
adata_seacell = ad.AnnData(dat_seacell.iloc[:, 17:30])
adata_seacell.obs['celltype'] = dat_seacell['Population L1']
sc.tl.pca(adata_seacell, n_comps=5)
sc.pp.neighbors(adata_seacell)
sc.tl.umap(adata_seacell)
sc.pl.umap(adata_seacell, color=['celltype'], save='_seacell_only.pdf')

adata_seacell.obs['celltype'] = np.array(dat_seacell['Population L1'])
sc.pl.umap(adata_seacell, color=['celltype'], save='_seacell_only.pdf')
