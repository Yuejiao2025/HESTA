# HESTA (Human Embryogenesis Spatiotemporal Transcriptomic Atlas)
The scripts utilized for dataset analysis in the paper “Spatiotemporal transcriptome atlas of human embryos after gastrulation”.

The comprehensive spatiotemporal atlas of gene expression during early human embryonic development is critical for insights into embryogenesis, organogenesis, and disease origins. We used Stereo-seq to generate a detailed atlas of gene expression across 77 sagittal sections of 13 whole human embryos ranging from Carnegie stage 12 to 23. Combined with single-nucleus RNA-seq, this spatial transcriptomic approach has elucidated gene expression patterns within defined cellular contexts, revealing the cellular heterogeneity that drives organ-specific differentiation. Our study has established a regulatory profile for the development of 50 organs and 198 substructures, while also identifying potential tissue-identity regulators. Notably, it uncovered previously uncharacterized gene functions in cardiac and brain development. The atlas not only substantiates and refines current understanding of human organ development but also highlights key organs/cell types susceptible to viral infections and genetic disorders. Furthermore, we characterized the dynamics of allelic gene expression within specific organs at different developmental stages. This work presents a groundbreaking compilation of genome-wide gene expression profiles for each spatially defined cell population, which can be visualized as a spatial display of the embryonic transcriptional landscape. These results offer the most thorough delineation of the spatiotemporal transcriptomic dynamics of human organogenesis.

# Instructions for use
The script Figure4_cell_count.py identifies cell types based on predefined markers and calculates their proportions across various brain substructures. 

The python scrips py_he_function_zy.py provides custom functions specifically designed to support analysis workflow in the jupyter notebook py_he_scripts_zy.ipynb.

The jupyter notebook py_he_scripts_zy.ipynb implements two core python analysis pipelines:
1. Spatial transcriptomics mapping of single-nucleus RNA-seq data using cell2location.
2. Implements a bin aggregation pipeline to identify spatially resolved developmental substructures.

The jupyter notebook R_he_scripts.zy.ipynb implements two core R analysis pipelines:
1. Utilizes Monocle3 for neural lineage-specific single-cell trajectory analysis.
2. Conducts rank-based gene set enrichment analysis using irGSEA at bin50 spatial resolution.

# Demo datasets
The data used in these scripts can be found at https://db.cngb.org/stomics/hesta/data/
