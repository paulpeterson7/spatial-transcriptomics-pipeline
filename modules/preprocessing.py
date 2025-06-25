"""
Preprocessing Module

This module handles normalization, feature selection, and dimensionality reduction
for spatial transcriptomics data following scanpy best practices.

Author: [Your Name]
Date: 2025
Project: Spatial Transcriptomics Analysis Pipeline
"""

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import anndata as ad
from typing import Optional, Dict, List, Tuple, Any, Union
import logging

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class SpatialPreprocessor:
    """
    Comprehensive preprocessing pipeline for spatial transcriptomics data.
    
    This class implements normalization, feature selection, and dimensionality
    reduction following established best practices.
    """
    
    def __init__(self, adata: ad.AnnData):
        """
        Initialize preprocessor.
        
        Parameters
        ----------
        adata : AnnData
            Quality-controlled spatial transcriptomics data
        """
        self.adata = adata.copy()
        self.preprocessing_log = {}
        
        # Store raw data if not already stored
        if self.adata.raw is None:
            self.adata.raw = self.adata
    
    def normalize_total(self, 
                       target_sum: float = 1e4,
                       exclude_highly_expressed: bool = True) -> ad.AnnData:
        """
        Normalize total counts per spot.
        
        Parameters
        ----------
        target_sum : float
            Target sum of counts per spot after normalization
        exclude_highly_expressed : bool
            Whether to exclude highly expressed genes from normalization
            
        Returns
        -------
        adata : AnnData
            Normalized data
        """
        logger.info(f"Normalizing total counts to {target_sum}")
        
        # Store pre-normalization statistics
        pre_norm_counts = np.array(self.adata.X.sum(axis=1)).flatten()
        self.preprocessing_log['pre_norm_total_counts'] = {
            'mean': pre_norm_counts.mean(),
            'std': pre_norm_counts.std(),
            'min': pre_norm_counts.min(),
            'max': pre_norm_counts.max()
        }
        
        # Normalize
        sc.pp.normalize_total(
            self.adata, 
            target_sum=target_sum,
            exclude_highly_expressed=exclude_highly_expressed,
            inplace=True
        )
        
        # Store post-normalization statistics
        post_norm_counts = np.array(self.adata.X.sum(axis=1)).flatten()
        self.preprocessing_log['post_norm_total_counts'] = {
            'mean': post_norm_counts.mean(),
            'std': post_norm_counts.std(),
            'target': target_sum
        }
        
        logger.info("Total count normalization complete")
        return self.adata
    
    def log_transform(self, base: Optional[float] = None) -> ad.AnnData:
        """
        Apply log(x+1) transformation.
        
        Parameters
        ----------
        base : float, optional
            Base of logarithm. If None, uses natural log
            
        Returns
        -------
        adata : AnnData
            Log-transformed data
        """
        logger.info("Applying log(x+1) transformation")
        
        # Store pre-log statistics
        pre_log_values = self.adata.X.data if hasattr(self.adata.X, 'data') else self.adata.X.flatten()
        nonzero_values = pre_log_values[pre_log_values > 0]
        
        self.preprocessing_log['pre_log_stats'] = {
            'min_nonzero': nonzero_values.min() if len(nonzero_values) > 0 else 0,
            'max': pre_log_values.max(),
            'mean_nonzero': nonzero_values.mean() if len(nonzero_values) > 0 else 0
        }
        
        # Apply log transformation
        if base is None:
            sc.pp.log1p(self.adata)
        else:
            sc.pp.log1p(self.adata, base=base)
        
        # Store post-log statistics
        post_log_values = self.adata.X.data if hasattr(self.adata.X, 'data') else self.adata.X.flatten()
        self.preprocessing_log['post_log_stats'] = {
            'min': post_log_values.min(),
            'max': post_log_values.max(),
            'mean': post_log_values.mean()
        }
        
        logger.info("Log transformation complete")
        return self.adata
    
    def find_highly_variable_genes(self,
                                  method: str = "seurat_v3",
                                  n_top_genes: int = 2000,
                                  min_mean: float = 0.0125,
                                  max_mean: float = 3,
                                  min_disp: float = 0.5,
                                  n_bins: int = 20) -> ad.AnnData:
        """
        Identify highly variable genes for downstream analysis.
        
        Parameters
        ----------
        method : str
            Method for HVG selection ('seurat', 'seurat_v3', 'cell_ranger')
        n_top_genes : int
            Number of top variable genes to select
        min_mean, max_mean : float
            Mean expression thresholds for 'seurat' method
        min_disp : float
            Minimum dispersion threshold for 'seurat' method
        n_bins : int
            Number of expression bins for 'seurat' method
            
        Returns
        -------
        adata : AnnData
            Data with HVG annotations added
        """
        logger.info(f"Finding highly variable genes using {method} method")
        
        if method == "seurat":
            sc.pp.highly_variable_genes(
                self.adata,
                flavor="seurat",
                min_mean=min_mean,
                max_mean=max_mean,
                min_disp=min_disp,
                n_bins=n_bins
            )
        elif method == "seurat_v3":
            sc.pp.highly_variable_genes(
                self.adata,
                flavor="seurat_v3",
                n_top_genes=n_top_genes
            )
        elif method == "cell_ranger":
            sc.pp.highly_variable_genes(
                self.adata,
                flavor="cell_ranger",
                n_top_genes=n_top_genes
            )
        else:
            raise ValueError(f"Unknown method: {method}")
        
        # Store HVG statistics
        hvg_count = self.adata.var['highly_variable'].sum()
        hvg_fraction = hvg_count / self.adata.n_vars
        
        self.preprocessing_log['hvg_selection'] = {
            'method': method,
            'n_hvg': hvg_count,
            'hvg_fraction': hvg_fraction,
            'total_genes': self.adata.n_vars
        }
        
        logger.info(f"Found {hvg_count} highly variable genes ({hvg_fraction:.1%} of total)")
        return self.adata
    
    def plot_hvg(self, 
                save_path: Optional[str] = None,
                figsize: Tuple[int, int] = (12, 8)) -> plt.Figure:
        """
        Plot highly variable gene selection results.
        
        Parameters
        ----------
        save_path : str, optional
            Path to save the plot
        figsize : tuple
            Figure size
            
        Returns
        -------
        fig : matplotlib.Figure
            The created figure
        """
        if 'highly_variable' not in self.adata.var.columns:
            raise ValueError("No HVG selection performed. Run find_highly_variable_genes() first.")
        
        logger.info("Plotting highly variable genes")
        
        # Create HVG plot
        fig = sc.pl.highly_variable_genes(self.adata, show=False)
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            logger.info(f"HVG plot saved to {save_path}")
        
        return fig
    
    def scale_data(self, 
                  max_value: Optional[float] = 10,
                  zero_center: bool = True,
                  use_hvg_only: bool = True) -> ad.AnnData:
        """
        Scale data to unit variance and optionally clip values.
        
        Parameters
        ----------
        max_value : float, optional
            Maximum value after scaling (clips extreme values)
        zero_center : bool
            Whether to center data around zero
        use_hvg_only : bool
            Whether to scale only highly variable genes
            
        Returns
        -------
        adata : AnnData
            Scaled data
        """
        logger.info("Scaling data")
        
        # Subset to HVG if requested
        if use_hvg_only and 'highly_variable' in self.adata.var.columns:
            adata_hvg = self.adata[:, self.adata.var.highly_variable].copy()
            logger.info(f"Scaling {adata_hvg.n_vars} highly variable genes")
        else:
            adata_hvg = self.adata.copy()
            logger.info(f"Scaling all {adata_hvg.n_vars} genes")
        
        # Store pre-scaling statistics
        pre_scale_data = adata_hvg.X
        if hasattr(pre_scale_data, 'toarray'):
            pre_scale_data = pre_scale_data.toarray()
        
        self.preprocessing_log['pre_scaling'] = {
            'mean': np.mean(pre_scale_data),
            'std': np.std(pre_scale_data),
            'min': np.min(pre_scale_data),
            'max': np.max(pre_scale_data)
        }
        
        # Scale data
        sc.pp.scale(adata_hvg, max_value=max_value, zero_center=zero_center)
        
        # Store post-scaling statistics
        post_scale_data = adata_hvg.X
        if hasattr(post_scale_data, 'toarray'):
            post_scale_data = post_scale_data.toarray()
        
        self.preprocessing_log['post_scaling'] = {
            'mean': np.mean(post_scale_data),
            'std': np.std(post_scale_data),
            'min': np.min(post_scale_data),
            'max': np.max(post_scale_data),
            'max_clip_value': max_value
        }
        
        # Update the main object
        if use_hvg_only:
            # Store scaled HVG data in a separate layer
            self.adata.layers['scaled_hvg'] = np.zeros_like(self.adata.X)
            hvg_mask = self.adata.var.highly_variable
            self.adata.layers['scaled_hvg'][:, hvg_mask] = adata_hvg.X
        else:
            self.adata.X = adata_hvg.X
        
        logger.info("Data scaling complete")
        return self.adata
    
    def principal_component_analysis(self,
                                   n_comps: int = 50,
                                   use_highly_variable: bool = True,
                                   svd_solver: str = 'arpack') -> ad.AnnData:
        """
        Perform principal component analysis.
        
        Parameters
        ----------
        n_comps : int
            Number of principal components to compute
        use_highly_variable : bool
            Whether to use only highly variable genes
        svd_solver : str
            SVD solver to use ('arpack', 'randomized', 'auto')
            
        Returns
        -------
        adata : AnnData
            Data with PCA results added
        """
        logger.info(f"Computing PCA with {n_comps} components")
        
        # Perform PCA
        sc.tl.pca(
            self.adata,
            n_comps=n_comps,
            use_highly_variable=use_highly_variable,
            svd_solver=svd_solver
        )
        
        # Store PCA statistics
        pca_var_ratio = self.adata.uns['pca']['variance_ratio']
        
        self.preprocessing_log['pca'] = {
            'n_components': n_comps,
            'variance_explained_top10': pca_var_ratio[:10].sum(),
            'variance_explained_top50': pca_var_ratio[:min(50, len(pca_var_ratio))].sum(),
            'use_hvg': use_highly_variable
        }
        
        logger.info(f"PCA complete. Top 10 PCs explain {pca_var_ratio[:10].sum():.1%} variance")
        return self.adata
    
    def plot_pca_variance(self,
                         n_pcs_show: int = 50,
                         save_path: Optional[str] = None,
                         figsize: Tuple[int, int] = (10, 6)) -> plt.Figure:
        """
        Plot PCA variance explained.
        
        Parameters
        ----------
        n_pcs_show : int
            Number of PCs to show in plot
        save_path : str, optional
            Path to save plot
        figsize : tuple
            Figure size
            
        Returns
        -------
        fig : matplotlib.Figure
            The created figure
        """
        if 'pca' not in self.adata.uns:
            raise ValueError("No PCA performed. Run principal_component_analysis() first.")
        
        logger.info("Plotting PCA variance explained")
        
        variance_ratio = self.adata.uns['pca']['variance_ratio']
        n_show = min(n_pcs_show, len(variance_ratio))
        
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)
        
        # Individual variance explained
        ax1.bar(range(1, n_show + 1), variance_ratio[:n_show])
        ax1.set_xlabel('Principal Component')
        ax1.set_ylabel('Variance Explained')
        ax1.set_title('Variance Explained by PC')
        ax1.grid(True, alpha=0.3)
        
        # Cumulative variance explained
        cumvar = np.cumsum(variance_ratio[:n_show])
        ax2.plot(range(1, n_show + 1), cumvar, 'o-')
        ax2.axhline(y=0.8, color='r', linestyle='--', alpha=0.7, label='80%')
        ax2.axhline(y=0.9, color='r', linestyle='--', alpha=0.7, label='90%')
        ax2.set_xlabel('Principal Component')
        ax2.set_ylabel('Cumulative Variance Explained')
        ax2.set_title('Cumulative Variance Explained')
        ax2.legend()
        ax2.grid(True, alpha=0.3)
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            logger.info(f"PCA plot saved to {save_path}")
        
        return fig
    
    def get_preprocessing_summary(self) -> str:
        """
        Get a summary of preprocessing steps performed.
        
        Returns
        -------
        summary : str
            Formatted preprocessing summary
        """
        if not self.preprocessing_log:
            return "No preprocessing performed yet."
        
        summary = []
        summary.append("=== PREPROCESSING SUMMARY ===")
        summary.append(f"Final data shape: {self.adata.n_obs} spots × {self.adata.n_vars} genes")
        
        # Normalization
        if 'post_norm_total_counts' in self.preprocessing_log:
            norm_info = self.preprocessing_log['post_norm_total_counts']
            summary.append(f"\n📏 NORMALIZATION:")
            summary.append(f"  Target sum: {norm_info['target']}")
            summary.append(f"  Mean total counts: {norm_info['mean']:.1f}")
        
        # Log transformation
        if 'post_log_stats' in self.preprocessing_log:
            log_info = self.preprocessing_log['post_log_stats']
            summary.append(f"\n📊 LOG TRANSFORMATION:")
            summary.append(f"  Value range: {log_info['min']:.3f} - {log_info['max']:.3f}")
        
        # HVG selection
        if 'hvg_selection' in self.preprocessing_log:
            hvg_info = self.preprocessing_log['hvg_selection']
            summary.append(f"\n🧬 HIGHLY VARIABLE GENES:")
            summary.append(f"  Method: {hvg_info['method']}")
            summary.append(f"  Selected: {hvg_info['n_hvg']} ({hvg_info['hvg_fraction']:.1%})")
        
        # Scaling
        if 'post_scaling' in self.preprocessing_log:
            scale_info = self.preprocessing_log['post_scaling']
            summary.append(f"\n⚖️ SCALING:")
            summary.append(f"  Mean: {scale_info['mean']:.3f}, Std: {scale_info['std']:.3f}")
            summary.append(f"  Clipped at: {scale_info['max_clip_value']}")
        
        # PCA
        if 'pca' in self.preprocessing_log:
            pca_info = self.preprocessing_log['pca']
            summary.append(f"\n🔄 PCA:")
            summary.append(f"  Components: {pca_info['n_components']}")
            summary.append(f"  Variance explained (top 10): {pca_info['variance_explained_top10']:.1%}")
        
        return "\n".join(summary)


def run_standard_preprocessing(adata: ad.AnnData,
                              target_sum: float = 1e4,
                              n_top_genes: int = 2000,
                              hvg_method: str = "seurat_v3",
                              max_value: float = 10,
                              n_pcs: int = 50,
                              plot_results: bool = True,
                              save_dir: Optional[str] = None) -> Tuple[ad.AnnData, SpatialPreprocessor]:
    """
    Run standard preprocessing pipeline.
    
    Parameters
    ----------
    adata : AnnData
        Input data (should be quality-controlled)
    target_sum : float
        Target sum for normalization
    n_top_genes : int
        Number of HVGs to select
    hvg_method : str
        HVG selection method
    max_value : float
        Maximum value for scaling
    n_pcs : int
        Number of principal components
    plot_results : bool
        Whether to generate plots
    save_dir : str, optional
        Directory to save plots
        
    Returns
    -------
    adata_processed : AnnData
        Preprocessed data
    preprocessor : SpatialPreprocessor
        Preprocessor object with logs
    """
    logger.info("Running standard preprocessing pipeline")
    
    # Initialize preprocessor
    preprocessor = SpatialPreprocessor(adata)
    
    # Step 1: Normalize total counts
    adata_processed = preprocessor.normalize_total(target_sum=target_sum)
    
    # Step 2: Log transform
    adata_processed = preprocessor.log_transform()
    
    # Step 3: Find highly variable genes
    adata_processed = preprocessor.find_highly_variable_genes(
        method=hvg_method,
        n_top_genes=n_top_genes
    )
    
    # Generate HVG plot
    if plot_results and save_dir:
        hvg_plot_path = f"{save_dir}/highly_variable_genes.png"
        preprocessor.plot_hvg(save_path=hvg_plot_path)
    
    # Step 4: Scale data
    adata_processed = preprocessor.scale_data(max_value=max_value)
    
    # Step 5: PCA
    adata_processed = preprocessor.principal_component_analysis(n_comps=n_pcs)
    
    # Generate PCA plot
    if plot_results and save_dir:
        pca_plot_path = f"{save_dir}/pca_variance.png"
        preprocessor.plot_pca_variance(save_path=pca_plot_path)
    
    logger.info("Preprocessing pipeline complete")
    return adata_processed, preprocessor


if __name__ == "__main__":
    # Example usage
    import sys
    sys.path.append('..')
    from modules.data_loader import load_and_validate_visium
    from modules.quality_control import run_comprehensive_qc
    
    try:
        # Load and QC data
        adata, _ = load_and_validate_visium("../data/example_data/")
        adata_qc, _ = run_comprehensive_qc(adata, plot_results=False)
        
        # Run preprocessing
        adata_processed, preprocessor = run_standard_preprocessing(
            adata_qc,
            plot_results=True,
            save_dir="../results/figures/"
        )
        
        # Print summary
        print(preprocessor.get_preprocessing_summary())
        
    except Exception as e:
        print(f"Error in preprocessing: {e}")
