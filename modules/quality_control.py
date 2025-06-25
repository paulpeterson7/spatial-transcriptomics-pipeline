"""
Quality Control Module

This module provides comprehensive quality control analysis for spatial 
transcriptomics data, following best practices from the scanpy tutorial.

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
from typing import Optional, Dict, List, Tuple, Any
import logging

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class SpatialQualityControl:
    """
    Comprehensive quality control analysis for spatial transcriptomics data.
    
    This class implements QC metrics calculation, filtering, and visualization
    following established best practices.
    """
    
    def __init__(self, adata: ad.AnnData):
        """
        Initialize QC analysis.
        
        Parameters
        ----------
        adata : AnnData
            Spatial transcriptomics data
        """
        self.adata = adata.copy()
        self.qc_metrics = {}
        self.filtering_stats = {}
        
    def calculate_qc_metrics(self, 
                           qc_vars: Optional[List[str]] = None) -> ad.AnnData:
        """
        Calculate standard QC metrics for spatial transcriptomics data.
        
        Parameters
        ----------
        qc_vars : list of str, optional
            Variables to calculate QC metrics for (e.g., ['mt'] for mitochondrial genes)
            
        Returns
        -------
        adata : AnnData
            Data with QC metrics added to .obs
        """
        logger.info("Calculating QC metrics...")
        
        # Default QC variables
        if qc_vars is None:
            qc_vars = ['mt']
        
        # Identify mitochondrial genes
        if 'mt' in qc_vars:
            self.adata.var["mt"] = self.adata.var_names.str.startswith("MT-")
            logger.info(f"Found {self.adata.var['mt'].sum()} mitochondrial genes")
        
        # Identify ribosomal genes
        if 'ribo' in qc_vars:
            self.adata.var["ribo"] = (
                self.adata.var_names.str.startswith("RPS") | 
                self.adata.var_names.str.startswith("RPL")
            )
            logger.info(f"Found {self.adata.var['ribo'].sum()} ribosomal genes")
        
        # Identify hemoglobin genes
        if 'hb' in qc_vars:
            self.adata.var["hb"] = self.adata.var_names.str.contains("^HB[^(P)]")
            logger.info(f"Found {self.adata.var['hb'].sum()} hemoglobin genes")
        
        # Calculate QC metrics using scanpy
        sc.pp.calculate_qc_metrics(
            self.adata, 
            qc_vars=qc_vars, 
            inplace=True,
            percent_top=None,
            log1p=False
        )
        
        # Additional spatial-specific metrics
        self._calculate_spatial_metrics()
        
        # Store QC summary
        self._summarize_qc_metrics()
        
        logger.info("QC metrics calculation complete")
        return self.adata
    
    def _calculate_spatial_metrics(self):
        """Calculate additional spatial-specific QC metrics."""
        
        # Spatial neighborhood metrics (if spatial coordinates available)
        if 'spatial' in self.adata.obsm:
            coords = self.adata.obsm['spatial']
            
            # Calculate distance to tissue center
            center_x = np.mean(coords[:, 0])
            center_y = np.mean(coords[:, 1])
            distances = np.sqrt(
                (coords[:, 0] - center_x)**2 + 
                (coords[:, 1] - center_y)**2
            )
            self.adata.obs['distance_to_center'] = distances
            
            # Calculate local density (spots within radius)
            # This is a simplified version - for production use proper spatial tools
            radius = np.percentile(distances, 10)  # Use 10th percentile as radius
            local_density = []
            for i, coord in enumerate(coords):
                dists_to_others = np.sqrt(
                    np.sum((coords - coord)**2, axis=1)
                )
                density = np.sum(dists_to_others < radius) - 1  # Exclude self
                local_density.append(density)
            
            self.adata.obs['local_spot_density'] = local_density
    
    def _summarize_qc_metrics(self):
        """Summarize calculated QC metrics."""
        
        obs_cols = self.adata.obs.columns
        
        # Basic metrics
        basic_metrics = [
            'total_counts', 'n_genes_by_counts', 
            'total_counts_mt', 'pct_counts_mt'
        ]
        
        # Additional metrics that might be present
        additional_metrics = [
            col for col in obs_cols 
            if any(keyword in col for keyword in ['ribo', 'hb', 'top_'])
        ]
        
        all_metrics = basic_metrics + additional_metrics
        available_metrics = [m for m in all_metrics if m in obs_cols]
        
        summary = {}
        for metric in available_metrics:
            values = self.adata.obs[metric]
            summary[metric] = {
                'mean': values.mean(),
                'median': values.median(),
                'std': values.std(),
                'min': values.min(),
                'max': values.max(),
                'q25': values.quantile(0.25),
                'q75': values.quantile(0.75)
            }
        
        self.qc_metrics = summary
    
    def plot_qc_distributions(self, 
                            figsize: Tuple[int, int] = (15, 10),
                            save_path: Optional[str] = None) -> plt.Figure:
        """
        Create comprehensive QC distribution plots.
        
        Parameters
        ----------
        figsize : tuple
            Figure size (width, height)
        save_path : str, optional
            Path to save the figure
            
        Returns
        -------
        fig : matplotlib.Figure
            The created figure
        """
        logger.info("Creating QC distribution plots...")
        
        # Determine available metrics
        obs_cols = self.adata.obs.columns
        
        # Core metrics to plot
        plot_metrics = []
        metric_configs = {
            'total_counts': {'title': 'Total UMI Counts', 'log': True},
            'n_genes_by_counts': {'title': 'Genes Detected', 'log': False},
            'pct_counts_mt': {'title': 'Mitochondrial %', 'log': False},
            'pct_counts_ribo': {'title': 'Ribosomal %', 'log': False},
            'distance_to_center': {'title': 'Distance to Center', 'log': False},
            'local_spot_density': {'title': 'Local Spot Density', 'log': False}
        }
        
        for metric, config in metric_configs.items():
            if metric in obs_cols:
                plot_metrics.append((metric, config))
        
        # Create subplots
        n_metrics = len(plot_metrics)
        cols = 3
        rows = (n_metrics + cols - 1) // cols
        
        fig, axes = plt.subplots(rows, cols, figsize=figsize)
        if rows == 1:
            axes = axes.reshape(1, -1)
        axes = axes.flatten()
        
        for i, (metric, config) in enumerate(plot_metrics):
            ax = axes[i]
            
            data = self.adata.obs[metric]
            
            # Create histogram
            if config['log'] and data.min() > 0:
                bins = np.logspace(np.log10(data.min()), np.log10(data.max()), 50)
                ax.hist(data, bins=bins, alpha=0.7, edgecolor='black', linewidth=0.5)
                ax.set_xscale('log')
            else:
                ax.hist(data, bins=50, alpha=0.7, edgecolor='black', linewidth=0.5)
            
            # Add statistics text
            stats_text = (
                f"Mean: {data.mean():.1f}\n"
                f"Median: {data.median():.1f}\n"
                f"Std: {data.std():.1f}"
            )
            ax.text(0.7, 0.7, stats_text, transform=ax.transAxes, 
                   verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
            
            ax.set_title(config['title'])
            ax.set_xlabel(metric.replace('_', ' ').title())
            ax.set_ylabel('Frequency')
            ax.grid(True, alpha=0.3)
        
        # Hide unused subplots
        for i in range(n_metrics, len(axes)):
            axes[i].set_visible(False)
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            logger.info(f"QC plots saved to {save_path}")
        
        return fig
    
    def identify_outliers(self, 
                         metrics: Optional[List[str]] = None,
                         n_mad: float = 3.0) -> Dict[str, np.ndarray]:
        """
        Identify outlier spots based on QC metrics using median absolute deviation.
        
        Parameters
        ----------
        metrics : list of str, optional
            Metrics to use for outlier detection
        n_mad : float
            Number of median absolute deviations to use as threshold
            
        Returns
        -------
        outliers : dict
            Dictionary mapping metric names to boolean arrays indicating outliers
        """
        logger.info("Identifying outlier spots...")
        
        if metrics is None:
            metrics = ['total_counts', 'n_genes_by_counts', 'pct_counts_mt']
        
        # Filter to available metrics
        available_metrics = [m for m in metrics if m in self.adata.obs.columns]
        
        outliers = {}
        
        for metric in available_metrics:
            values = self.adata.obs[metric]
            
            # Calculate median absolute deviation
            median = values.median()
            mad = np.median(np.abs(values - median))
            
            # Define thresholds
            if metric == 'pct_counts_mt':
                # For mitochondrial %, only flag high values
                threshold_high = median + n_mad * mad
                outlier_mask = values > threshold_high
            else:
                # For other metrics, flag both high and low values
                threshold_low = median - n_mad * mad
                threshold_high = median + n_mad * mad
                outlier_mask = (values < threshold_low) | (values > threshold_high)
            
            outliers[metric] = outlier_mask
            
            logger.info(f"{metric}: {outlier_mask.sum()} outliers ({outlier_mask.mean():.1%})")
        
        return outliers
    
    def filter_spots(self,
                    min_counts: Optional[int] = None,
                    max_counts: Optional[int] = None,
                    min_genes: Optional[int] = None,
                    max_genes: Optional[int] = None,
                    max_pct_mt: Optional[float] = None,
                    tissue_only: bool = True) -> ad.AnnData:
        """
        Filter spots based on QC criteria.
        
        Parameters
        ----------
        min_counts : int, optional
            Minimum total UMI counts per spot
        max_counts : int, optional  
            Maximum total UMI counts per spot
        min_genes : int, optional
            Minimum number of genes per spot
        max_genes : int, optional
            Maximum number of genes per spot
        max_pct_mt : float, optional
            Maximum mitochondrial gene percentage
        tissue_only : bool
            Whether to keep only spots marked as 'in tissue'
            
        Returns
        -------
        adata_filtered : AnnData
            Filtered data
        """
        logger.info("Filtering spots based on QC criteria...")
        
        adata_filtered = self.adata.copy()
        initial_spots = adata_filtered.n_obs
        
        self.filtering_stats = {'initial_spots': initial_spots}
        
        # Filter by tissue annotation
        if tissue_only and 'in_tissue' in adata_filtered.obs:
            tissue_spots = adata_filtered.obs['in_tissue'].sum()
            adata_filtered = adata_filtered[adata_filtered.obs['in_tissue']].copy()
            logger.info(f"Tissue filtering: {tissue_spots}/{initial_spots} spots retained")
            self.filtering_stats['after_tissue'] = adata_filtered.n_obs
        
        # Filter by total counts
        if min_counts is not None:
            before = adata_filtered.n_obs
            adata_filtered = adata_filtered[adata_filtered.obs['total_counts'] >= min_counts].copy()
            logger.info(f"Min counts filter: {adata_filtered.n_obs}/{before} spots retained")
            self.filtering_stats['after_min_counts'] = adata_filtered.n_obs
            
        if max_counts is not None:
            before = adata_filtered.n_obs
            adata_filtered = adata_filtered[adata_filtered.obs['total_counts'] <= max_counts].copy()
            logger.info(f"Max counts filter: {adata_filtered.n_obs}/{before} spots retained")
            self.filtering_stats['after_max_counts'] = adata_filtered.n_obs
        
        # Filter by number of genes
        if min_genes is not None:
            before = adata_filtered.n_obs
            adata_filtered = adata_filtered[adata_filtered.obs['n_genes_by_counts'] >= min_genes].copy()
            logger.info(f"Min genes filter: {adata_filtered.n_obs}/{before} spots retained")
            self.filtering_stats['after_min_genes'] = adata_filtered.n_obs
            
        if max_genes is not None:
            before = adata_filtered.n_obs
            adata_filtered = adata_filtered[adata_filtered.obs['n_genes_by_counts'] <= max_genes].copy()
            logger.info(f"Max genes filter: {adata_filtered.n_obs}/{before} spots retained")
            self.filtering_stats['after_max_genes'] = adata_filtered.n_obs
        
        # Filter by mitochondrial percentage
        if max_pct_mt is not None and 'pct_counts_mt' in adata_filtered.obs:
            before = adata_filtered.n_obs
            adata_filtered = adata_filtered[adata_filtered.obs['pct_counts_mt'] <= max_pct_mt].copy()
            logger.info(f"MT% filter: {adata_filtered.n_obs}/{before} spots retained")
            self.filtering_stats['after_mt_filter'] = adata_filtered.n_obs
        
        self.filtering_stats['final_spots'] = adata_filtered.n_obs
        self.filtering_stats['spots_removed'] = initial_spots - adata_filtered.n_obs
        self.filtering_stats['retention_rate'] = adata_filtered.n_obs / initial_spots
        
        logger.info(f"Filtering complete: {adata_filtered.n_obs}/{initial_spots} spots retained ({adata_filtered.n_obs/initial_spots:.1%})")
        
        return adata_filtered
    
    def filter_genes(self,
                    min_cells: int = 10,
                    min_counts: int = 1) -> ad.AnnData:
        """
        Filter genes based on expression criteria.
        
        Parameters
        ----------
        min_cells : int
            Minimum number of cells expressing the gene
        min_counts : int
            Minimum total counts for the gene
            
        Returns
        -------
        adata_filtered : AnnData
            Data with filtered genes
        """
        logger.info("Filtering genes...")
        
        initial_genes = self.adata.n_vars
        
        # Filter by minimum cells
        sc.pp.filter_genes(self.adata, min_cells=min_cells)
        
        logger.info(f"Gene filtering: {self.adata.n_vars}/{initial_genes} genes retained")
        
        return self.adata
    
    def plot_spatial_qc(self,
                       metrics: Optional[List[str]] = None,
                       figsize: Tuple[int, int] = (15, 10),
                       save_path: Optional[str] = None) -> plt.Figure:
        """
        Plot QC metrics in spatial coordinates.
        
        Parameters
        ----------
        metrics : list of str, optional
            QC metrics to plot spatially
        figsize : tuple
            Figure size
        save_path : str, optional
            Path to save figure
            
        Returns
        -------
        fig : matplotlib.Figure
            The created figure
        """
        logger.info("Creating spatial QC plots...")
        
        if metrics is None:
            metrics = ['total_counts', 'n_genes_by_counts']
        
        # Filter to available metrics
        available_metrics = [m for m in metrics if m in self.adata.obs.columns]
        
        if not available_metrics:
            raise ValueError("No valid metrics found for spatial plotting")
        
        # Create spatial plots using scanpy
        sc.pl.spatial(
            self.adata,
            color=available_metrics,
            img_key="hires" if "hires" in self.adata.uns.get('spatial', {}).get(list(self.adata.uns.get('spatial', {}).keys())[0] if self.adata.uns.get('spatial') else '', {}).get('images', {}) else None,
            size=1.5,
            alpha=0.7,
            ncols=2,
            save=save_path.replace('.png', '_spatial_qc.png') if save_path else None
        )
        
        return plt.gcf()
    
    def get_qc_summary(self) -> str:
        """
        Get a comprehensive QC summary report.
        
        Returns
        -------
        summary : str
            Formatted QC summary
        """
        if not self.qc_metrics:
            return "No QC metrics calculated. Run calculate_qc_metrics() first."
        
        summary = []
        summary.append("=== QUALITY CONTROL SUMMARY ===")
        summary.append(f"Data shape: {self.adata.n_obs} spots × {self.adata.n_vars} genes")
        
        # QC metrics summary
        summary.append(f"\n📊 QC METRICS:")
        for metric, stats in self.qc_metrics.items():
            summary.append(f"  {metric}:")
            summary.append(f"    Mean: {stats['mean']:.1f}, Median: {stats['median']:.1f}")
            summary.append(f"    Range: {stats['min']:.1f} - {stats['max']:.1f}")
        
        # Filtering summary
        if self.filtering_stats:
            summary.append(f"\n🔧 FILTERING RESULTS:")
            for step, count in self.filtering_stats.items():
                if step != 'retention_rate':
                    summary.append(f"  {step}: {count}")
            summary.append(f"  Overall retention: {self.filtering_stats.get('retention_rate', 0):.1%}")
        
        return "\n".join(summary)


def run_comprehensive_qc(adata: ad.AnnData,
                        min_counts: int = 1000,
                        max_counts: int = 50000,
                        min_genes: int = 500,
                        max_pct_mt: float = 20.0,
                        plot_results: bool = True,
                        save_dir: Optional[str] = None) -> Tuple[ad.AnnData, SpatialQualityControl]:
    """
    Run comprehensive quality control analysis.
    
    Parameters
    ----------
    adata : AnnData
        Input spatial transcriptomics data
    min_counts, max_counts : int
        Count filtering thresholds
    min_genes : int
        Minimum genes per spot
    max_pct_mt : float
        Maximum mitochondrial percentage
    plot_results : bool
        Whether to generate QC plots
    save_dir : str, optional
        Directory to save plots
        
    Returns
    -------
    adata_qc : AnnData
        Quality-controlled data
    qc_analyzer : SpatialQualityControl
        QC analyzer object with results
    """
    # Initialize QC analyzer
    qc_analyzer = SpatialQualityControl(adata)
    
    # Calculate QC metrics
    adata_qc = qc_analyzer.calculate_qc_metrics()
    
    # Generate plots
    if plot_results:
        # QC distribution plots
        qc_plot_path = f"{save_dir}/qc_distributions.png" if save_dir else None
        qc_analyzer.plot_qc_distributions(save_path=qc_plot_path)
        
        # Spatial QC plots
        spatial_plot_path = f"{save_dir}/spatial_qc.png" if save_dir else None
        qc_analyzer.plot_spatial_qc(save_path=spatial_plot_path)
    
    # Apply filtering
    adata_filtered = qc_analyzer.filter_spots(
        min_counts=min_counts,
        max_counts=max_counts,
        min_genes=min_genes,
        max_pct_mt=max_pct_mt
    )
    
    # Filter genes
    adata_final = qc_analyzer.filter_genes(min_cells=10)
    
    return adata_final, qc_analyzer


if __name__ == "__main__":
    # Example usage
    import sys
    sys.path.append('..')
    from modules.data_loader import load_and_validate_visium
    
    try:
        # Load data
        adata, _ = load_and_validate_visium("../data/example_data/")
        
        # Run QC
        adata_qc, qc_analyzer = run_comprehensive_qc(
            adata,
            plot_results=True,
            save_dir="../results/figures/"
        )
        
        # Print summary
        print(qc_analyzer.get_qc_summary())
        
    except Exception as e:
        print(f"Error in QC analysis: {e}")