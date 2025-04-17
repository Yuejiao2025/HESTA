import scanpy as sc
import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import math
from scipy.sparse import isspmatrix

dest_dir = "/data/work/human_embryo/"

sample_table = pd.read_csv(dest_dir + "sample_info.tab", sep="\t", header=0)
chip_dict = pd.Series(sample_table['chip'].values, index=sample_table['sample'].values).to_dict()
bin50_dict = pd.Series(sample_table['bin50'].values, index=sample_table['sample'].values).to_dict()
cellbin_dict = pd.Series(sample_table['cellbin'].values, index=sample_table['sample'].values).to_dict()
snRNA_dict = pd.Series(sample_table['snRNA'].values, index=sample_table['sample'].values).to_dict()


def initial_path(sample, dest_dir, cut):
    chip = chip_dict[sample]
    sample_name = sample + "_" + chip
    bin50_folder = dest_dir + "/bin50/"
    cut_folder = dest_dir + cut + "/"
    return chip, sample_name, bin50_folder, cut_folder


def bin_cluster_mean(adata_sc, domain_name):
    # adata_sc.X = adata_sc.layers["norm"].copy()
    adata_sc.obs[domain_name] = adata_sc.obs[domain_name].astype('category')
    adata_sc.varm['means_per_cluster'] = np.empty((len(adata_sc.var_names), len(adata_sc.obs[domain_name].cat.categories)))
    means_per_cluster = pd.DataFrame(columns=adata_sc.var_names, index=adata_sc.obs[domain_name].cat.categories)
    for clust in adata_sc.obs[domain_name].cat.categories:
        means_per_cluster.loc[clust] = adata_sc[adata_sc.obs[domain_name].isin([clust]), :].X.mean(0)
    adata_sc.varm['means_per_cluster'] = means_per_cluster.T.to_numpy()
    inf_aver = means_per_cluster.T
    return inf_aver


def deconv_assign_cell_label(bin_name, deconv):
    adata_cellbin = bin_name.copy()
    cell_types = []
    scores_list = []
    if deconv in "cell2location":
        adata_cellbin.obs[adata_cellbin.uns['mod']['factor_names']] = adata_cellbin.obsm[
            'q05_cell_abundance_w_sf']
        pred_label = adata_cellbin.obsm['q05_cell_abundance_w_sf'].idxmax(axis=1).values  # + "_cluster"
        scores_list = adata_cellbin.obsm['q05_cell_abundance_w_sf'].max(axis=1).values  # np.max(tmp, axis=1)
        adata_cellbin.obs['cell_label'] = [str(i).removeprefix("q05cell_abundance_w_sf_") for i in pred_label]
        cell_types = adata_cellbin.uns['mod']['factor_names']  # np.unique(adata_cellbin.obs['cell_label'])  # adata_cellbin.uns['mod']['factor_names']
        # adata_cellbin.obs[cell_types.astype(str)] = adata_cellbin.obs[cell_types]
    return adata_cellbin, cell_types, scores_list


def plot_scanpy_spatial_cluster(data, bin_path, domain_key_list, domain_file_name, sample):
    if "cellbin" in bin_path:
        sc.pl.spatial(data, img_key=None, spot_size=20, library_id=None, color=domain_key_list, scale_factor=1, frameon=False, ncols=3, show=False)  # cmap='tab20',
    else:
        sc.pl.spatial(data, img_key=None, spot_size=50, library_id=None, color=domain_key_list, scale_factor=1, frameon=False, ncols=3, show=False)  # cmap='tab20',
    plt.savefig(bin_path + "." + domain_file_name + ".spatial_cluster.png", dpi=300, bbox_inches='tight', transparent=True)   
        

def plot_scanpy_umap_cluster(bin_name, bin_path, domain_key_list, domain_file_name):
    data = bin_name
    sc.pl.umap(data, color=domain_key_list, ncols=5, wspace=0.4, frameon=False, show=False)  #  cmap='tab20',
    plt.savefig(bin_path + "." + domain_file_name + ".umap.pdf", dpi=300, bbox_inches='tight', transparent=True)


def plot_scanpy_spatial_cluster_palette(data, bin_path, domain_key, domain_file_name, sample, color_tab_path):
    data.obs[domain_key] = data.obs[domain_key].astype('category')
    domain_key_list = [domain_key]
    cluster_table = pd.read_csv(color_tab_path, sep="\t", header=0)
    palette_dict = pd.Series(cluster_table['colors'].values, index=cluster_table['clusters'].values.astype(str)).to_dict()
    # fig, axs = plt.subplots(1, 1, figsize=(6, 6))
    if "cellbin" in bin_path:
        sc.pl.spatial(data, img_key=None, spot_size=20, library_id=None, color=domain_key_list, scale_factor=1, frameon=False, ncols=3, show=False, palette=palette_dict)  # cmap='tab20',
    else:
        sc.pl.spatial(data, img_key=None, spot_size=50, library_id=None, color=domain_key_list, scale_factor=1, frameon=False, ncols=3, show=False, palette=palette_dict)  # cmap='tab20',
    plt.savefig(bin_path + "." + domain_file_name + ".spatial_cluster.pdf", dpi=600, bbox_inches='tight', transparent=True)
    
        
def bin_louvain_cluster(bin_name, bin_path, sample, res_array):
    import dynamo as dyn
    data = bin_name.copy()
    res_list = []
    for res in res_array:
        # Vanilla louvain clustering
        dyn.tl.louvain(data, resolution=res, result_key="louvain_" + str(res) + "_clusters")
        res_list.append("louvain_" + str(res) + "_clusters")
    plot_scanpy_umap_cluster(data, bin_path, res_list, "louvain_clusters")
    if sample != "all":
        plot_scanpy_spatial_cluster(data, bin_path, res_list, "louvain_clusters", sample)
    data.write(bin_path + ".louvain_clusters.h5ad")
    return data


def scanpy_bin_highly_variable(bin_name):
    data = bin_name.copy()
    sc.pp.normalize_total(data, inplace=True, target_sum=1e4)
    sc.pp.log1p(data, base=None, layer=None)
    data.layers["norm"] = data.X
    sc.pp.highly_variable_genes(data, flavor="seurat", n_top_genes=5000, layer=None, subset=False, batch_key=None, replace=False)
    sc.pp.pca(data, n_comps=50, use_highly_variable=None, svd_solver='arpack')
    return data
    
    
def bin_merge_by_domain_cluster_denovo(sample_list, final_bin_path, domain_name, domain_key_list, domain_file_name, plot_cluster_list, dest_dir, cut):
    # import scanorama
    import scanpy as sc
    # import dynamo as dyn
    import anndata as ad
    adata_sts = []
    for sample in sample_list:
        print(sample + " running!")
        chip, sample_name, bin50_folder, cut_folder, register_folder = initial_path(sample, dest_dir, cut)
        # adata_bin50_path = bin50_folder + "h5ad/" + sample + ".bin50"
        adata_bin50_path = "/data/input/Files/human_embryo_data/human_embryo_h5ad/" + sample + ".bin50"
        if "cellbin" in final_bin_path:
            adata_bin50_path = cut_folder + sample + ".cellbin"
        adata_bin50 = sc.read(adata_bin50_path + ".h5ad")
        adata_bin50.X = adata_bin50.layers['counts'].copy()
        if sample not in adata_bin50.obs_names[1]:
            adata_bin50.obs_names = sample + "_" + adata_bin50.obs_names
        adata_bin50_domain = bin_subset_domain(adata_bin50, domain_name, domain_key_list)
        if len(adata_bin50_domain.obs_names) > 0:
            adata_bin50_domain = scanpy_bin_highly_variable(adata_bin50_domain)
            adata_bin50_domain.X = adata_bin50_domain.layers['log']
            adata_sts.append(adata_bin50_domain)
    final_domain_bin_path = final_bin_path + "." + domain_file_name
    # before correction
    adata_bin50_domain_merge = sc.concat(adata_sts, join='outer', fill_value=0, merge="same") # uns_merge="same"
    # sc.pp.neighbors(adata_bin50_domain_merge)
    # sc.tl.umap(adata_bin50_domain_merge)
    adata_bin50_domain_merge.X = adata_bin50_domain_merge.layers["counts"]
    adata_bin50_domain_merge.write(final_domain_bin_path + ".before_correct.h5ad")
    plot_scanpy_umap_cluster(adata_bin50_domain_merge, final_domain_bin_path + ".before_correct", plot_cluster_list, "before_correct")
    return adata_bin50_domain_merge


def scvi_bin_merge_tutorial(adata_scvi, batch_key, final_bin_path, plot_cluster_list):
    # https://www.sc-best-practices.org/cellular_structure/integration.html#id219
    import scvi
    import scanpy as sc
    # https://nbisweden.github.io/workshop-scRNAseq/labs/scanpy/scanpy_03_integration.html
    scvi.model.SCVI.setup_anndata(adata_scvi, layer="counts", batch_key=batch_key)
    model_scvi = scvi.model.SCVI(adata_scvi)  # , n_layers=2, n_latent=30, gene_likelihood="zinb", n_hidden = 60, dropout_rate = 0.1, dispersion = 'gene-batch',
    # model_scvi.view_anndata_setup()
    max_epochs_scvi = 5 # np.min([round((20000 / adata_scvi.n_obs) * 400), 400])
    # https://discourse.scverse.org/t/error-when-using-scvi-model-train/2423
    kwargs = {'lr': 0.0023}
    model_scvi.train(max_epochs = int(max_epochs_scvi), train_size = 0.9, validation_size = 0.1, #accelerator='gpu', 
                     check_val_every_n_epoch=1, early_stopping=True, early_stopping_patience=10, early_stopping_monitor="elbo_validation", plan_kwargs = kwargs)
    # Extracting the embedding
    adata_scvi.obsm["X_scVI"] = model_scvi.get_latent_representation()
    # Calculate a batch-corrected UMAP
    sc.pp.neighbors(adata_scvi, use_rep="X_scVI")
    sc.tl.umap(adata_scvi)
    sc.tl.tsne(adata_scvi)
    adata_scvi.write_h5ad(final_bin_path + ".h5ad")
    plot_scanpy_umap_cluster(adata_scvi, final_bin_path, plot_cluster_list, "plot_clusters")
    # plot_scanpy_tsne_cluster(adata_scvi, final_bin_path, plot_cluster_list, "plot_clusters")
    return adata_scvi


def save_bin_clusters_tab(data, bin_path, clusters_list, domain_file_name):
    domain_table = pd.DataFrame(data.obs[clusters_list])
    domain_table.index = data.obs.index
    domain_table.index.name = "cell"
    domain_table.to_csv(bin_path + "." + domain_file_name + "_anno.tab", sep='\t', index=True, header=True)
    return domain_table


def bin_rank_genes(adata, bin_path, domain_key):
    if 'log1p' in list(adata.uns.keys()):
        del adata.uns['log1p']
    adata.obs[domain_key] = adata.obs[domain_key].astype('category')
    adata.obs["clusters"] = adata.obs[domain_key].cat.as_ordered()
    sc.tl.dendrogram(adata, 'clusters')
    sc.pl.dendrogram(adata, 'clusters', show=False)
    plt.savefig(bin_path + ".dendrogram.pdf", dpi=300, bbox_inches='tight', transparent=True)
    sc.pl.correlation_matrix(adata, 'clusters', figsize=(6, 6), show=False)
    plt.savefig(bin_path + ".correlation_matrix.pdf", dpi=300, bbox_inches='tight', transparent=True)
    sc.tl.rank_genes_groups(adata, "clusters", method="t-test")  # , key_added = "t-test"
    sc.pl.rank_genes_groups(adata, n_genes=20, sharey=False, show=False)
    plt.savefig(bin_path + ".rank_genes.pdf", dpi=300, bbox_inches='tight', transparent=True)
    sc.pl.rank_genes_groups_stacked_violin(adata, n_genes=10, cmap='bwr', show=False)  # viridis_r # key="wilcoxon", groupby="clusters"
    plt.savefig(bin_path + ".rank_genes.stacked_violin.pdf", dpi=100, bbox_inches='tight', transparent=True)
    df = get_markers(adata, "clusters", 0.01)
    df.to_csv(bin_path + ".rank_genes.sig.tab", sep='\t', index=True, header=True)
    result = adata.uns['rank_genes_groups']
    groups = result['names'].dtype.names
    df_gene = pd.DataFrame({"cluster_" + group: result[key][group] for group in groups for key in ['names']}).head(100)  # 'pvals'
    df_gene.to_csv(bin_path + ".rank_genes.top100_genes.tab", sep='\t', index=True, header=True)
    adata.write(bin_path + ".rank.h5ad")
    return adata

