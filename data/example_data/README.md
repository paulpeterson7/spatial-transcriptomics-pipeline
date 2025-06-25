# 📥 Breast Cancer Spatial Dataset (10x Genomics)

This example uses:
- **Filtered expression matrix (.h5)**
- **Spatial image and metadata (.tar.gz)**

Download from:
- https://cf.10xgenomics.com/samples/spatial-exp/1.3.0/Visium_Human_Breast_Cancer/Visium_Human_Breast_Cancer_filtered_feature_bc_matrix.h5
- https://cf.10xgenomics.com/samples/spatial-exp/1.3.0/Visium_Human_Breast_Cancer/Visium_Human_Breast_Cancer_spatial.tar.gz

Extract both to `data/example_data/` and load via:

```python
import scanpy as sc
adata = sc.read_visium("data/example_data/", count_file='filtered_feature_bc_matrix.h5')
