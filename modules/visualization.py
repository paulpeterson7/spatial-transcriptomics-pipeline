"""
Visualization Module

This module provides comprehensive visualization functions for spatial
transcriptomics data, including spatial plots, heatmaps, and summary figures.

Author: [Your Name]
Date: 2025
Project: Spatial Transcriptomics Analysis Pipeline
"""

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import seaborn as sns
import anndata as ad
from typing import Optional, Dict, List, Tuple, Any, Union
import logging
from pathlib import Path
import warnings

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class SpatialVisualizer:
    """
    Comprehensive visualization toolkit for spatial transcriptomics data.
    
    This class provides methods for creating publication-ready figures
    including spatial plots, heatmaps, and integrated analysis visualizations.
    """
    
    def __init__(self, adata: ad.AnnData, save_dir: Optional[str] = None):
        """
        Initialize visualizer.
        
        Parameters
        ----------
        adata : AnnData
            Spatial transcriptomics data with analysis results
        save_dir : str, optional
            Directory to save figures
        """
        self.adata = adata.copy()
        self.save_dir = Path(save_dir) if save_dir else Path("./figures")
        self.save_dir.mkdir(parents=True, exist_ok=True)
        
        # Configure matplotlib for publication-quality figures
        plt.rcParams.update({
            'figure.dpi': 300,
            'savefig.dpi': 300,
            'font.size': 12,
            'axes.titlesize': 14,
            'axes.labelsize': 12,
            'xtick.labelsize': 10,
            'ytick.labelsize': 10,
            'legend.fontsize': 10,
            'figure.titlesize': 16
        })
    
    def plot_spatial_overview(self,
                             figsize: Tuple[int, int] = (20, 12),
                             save_name: str = "spatial_overview.png") -> plt.Figure:
        """
        Create comprehensive spatial overview figure.
        
        Parameters
        ----------
        figsize : tuple
            Figure size (width, height)
        save_name : str
            Filename for saving
            
        Returns
        -------
        fig : matplotlib.Figure
            The created figure
        """
        logger.info("Creating spatial overview figure")
        
        # Check available metadata
        obs_cols = self.adata.obs.columns
        spatial_metrics = ['total_counts', 'n_genes_by_counts']
        cluster_col = None
        
        # Find clustering results
        for col in ['spatial_clusters', 'leiden', 'louvain']:
            if col in obs_cols:
                cluster_col = col
                break
        
        # Determine subplot layout
        n_plots = len([col for col in spatial_metrics if col in obs_cols])
        if cluster_col:
            n_plots += 1
        
        cols = 3
        rows = (n_plots + cols - 1) // cols
        
        fig, axes = plt.subplots(rows, cols, figsize=figsize)
        if rows == 1:
            axes = axes.reshape(1, -1)
        axes = axes.flatten()
        
        plot_idx = 0
        
        # Plot spatial metrics
        for metric in spatial_metrics:
            if metric in obs_cols and plot_idx < len(axes):
                self._plot_spatial_metric(
                    metric=metric,
                    ax=axes[plot_idx],
                    title=metric.replace('_', ' ').title()
                )
                plot_idx += 1
        
        # Plot clusters if available
        if cluster_col and plot_idx < len(axes):
            self._plot_spatial_categorical(
                column=cluster_col,
                ax=axes[plot_idx],
                title="Spatial Clusters"
            )
            plot_idx += 1
        
        # Hide unused subplots
        for i in range(plot_idx, len(axes)):
            axes[i].set_visible(False)
        
        plt.tight_layout()
        
        # Save figure
        save_path = self.save_dir / save_name
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        logger.info(f"Spatial overview saved to {save_path}")
        
        return fig
    
    def plot_gene_expression_panel(self,
                                  genes: List[str],
                                  ncols: int = 3,
                                  spot_size: float = 2.0,
                                  alpha_img: float = 0.7,
                                  cmap: str = 'hot',
                                  figsize_per_plot: Tuple[int, int] = (6, 6),
                                  save_name: str = "gene_expression_panel.png") -> plt.Figure:
        """
        Create multi-gene expression panel.
        
        Parameters
        ----------
        genes : list of str
            Genes to visualize
        ncols : int
            Number of columns in subplot grid
        spot_size : float
            Size of spatial spots
        alpha_img : float
            Transparency of background image
        cmap : str
            Colormap for expression
        figsize_per_plot : tuple
            Size of each subplot
        save_name : str
            Filename for saving
            
        Returns
        -------
        fig : matplotlib.Figure
            The created figure
        """
        logger.info(f"Creating gene expression panel for {len(genes)} genes")
        
        # Filter to available genes
        available_genes = [g for g in genes if g in self.adata.var_names]
        if not available_genes:
            raise ValueError("No valid genes found in dataset")
        
        logger.info(f"Plotting {len(available_genes)}/{len(genes)} available genes")
        
        # Calculate subplot layout
        nrows = (len(available_genes) + ncols - 1) // ncols
        figsize = (figsize_per_plot[0] * ncols, figsize_per_plot[1] * nrows)
        
        fig, axes = plt.subplots(nrows, ncols, figsize=figsize)
        if nrows == 1 and ncols == 1:
            axes = [axes]
        elif nrows == 1:
            axes = axes.reshape(1, -1)
        axes = axes.flatten()
        
        # Plot each gene
        for i, gene in enumerate(available_genes):
            if i < len(axes):
                self._plot_spatial_gene(
                    gene=gene,
                    ax=axes[i],
                    spot_size=spot_size,
                    alpha_img=alpha_img,
                    cmap=cmap
                )
        
        # Hide unused subplots
        for i in range(len(available_genes), len(axes)):
            axes[i].set_visible(False)
        
        plt.tight_layout()
        
        # Save figure
        save_path = self.save_dir / save_name
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        logger.info(f"Gene expression panel saved to {save_path}")
        
        return fig
    
    def plot_gene_panels_summary(self,
                                gene_panels: Dict[str, List[str]],
                                figsize: Tuple[int, int] = (20, 15),
                                save_name: str = "gene_panels_summary.png") -> plt.Figure:
        """
        Create summary visualization of gene panel scores.
        
        Parameters
        ----------
        gene_panels : dict
            Dictionary mapping panel names to gene lists
        figsize : tuple
            Figure size
        save_name : str
            Filename for saving
            
        Returns
        -------
        fig : matplotlib.Figure
            The created figure
        """
        logger.info(f"Creating gene panels summary for {len(gene_panels)} panels")
        
        # Check which panels have scores calculated
        available_panels = {}
        for panel_name, genes in gene_panels.items():
            score_col = f"{panel_name}_score"
            if score_col in self.adata.obs.columns:
                available_panels[panel_name] = score_col
        
        if not available_panels:
            raise ValueError("No gene panel scores found. Run spatial analysis first.")
        
        # Create subplot layout
        n_panels = len(available_panels)
        cols = 3
        rows = (n_panels + cols - 1) // cols
        
        fig, axes = plt.subplots(rows, cols + 1, figsize=figsize)  # +1 for correlation heatmap
        if rows == 1:
            axes = axes.reshape(1, -1)
        
        # Plot spatial distribution of each panel
        for i, (panel_name, score_col) in enumerate(available_panels.items()):
            row = i // cols
            col = i % cols
            
            self._plot_spatial_metric(
                metric=score_col,
                ax=axes[row, col],
                title=f"{panel_name.title()} Panel Score",
                cmap='viridis'
            )
        
        # Create correlation heatmap in the last column
        panel_scores = pd.DataFrame({
            name: self.adata.obs[score_col] 
            for name, score_col in available_panels.items()
        })
        
        # Calculate correlation matrix
        corr_matrix = panel_scores.corr()
        
        # Plot correlation heatmap
        heatmap_ax = axes[0, -1] if rows == 1 else axes[0, -1]
        sns.heatmap(
            corr_matrix,
            annot=True,
            cmap='RdBu_r',
            center=0,
            square=True,
            ax=heatmap_ax,
            cbar_kws={'label': 'Correlation'}
        )
        heatmap_ax.set_title('Panel Score Correlations')
        
        # Hide unused subplot spaces
        for row in range(rows):
            for col in range(cols):
                idx = row * cols + col
                if idx >= n_panels:
                    axes[row, col].set_visible(False)
        
        # Hide remaining correlation plot spaces
        for row in range(1, rows):
            axes[row, -1].set_visible(False)
        
        plt.tight_layout()
        
        # Save figure
        save_path = self.save_dir / save_name
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        logger.info(f"Gene panels summary saved to {save_path}")
        
        return fig
    
    def plot_clustering_summary(self,
                               cluster_col: str = 'spatial_clusters',
                               figsize: Tuple[int, int] = (20, 12),
                               save_name: str = "clustering_summary.png") -> plt.Figure:
        """
        Create clustering analysis summary figure.
        
        Parameters
        ----------
        cluster_col : str
            Column containing cluster assignments
        figsize : tuple
            Figure size
        save_name : str
            Filename for saving
            
        Returns
        -------
        fig : matplotlib.Figure
            The created figure
        """
        logger.info(f"Creating clustering summary for {cluster_col}")
        
        if cluster_col not in self.adata.obs.columns:
            raise ValueError(f"Column {cluster_col} not found")
        
        fig = plt.figure(figsize=figsize)
        gs = fig.add_gridspec(3, 4, height_ratios=[2, 1, 1], width_ratios=[1, 1, 1, 1])
        
        # Main spatial plot
        ax_spatial = fig.add_subplot(gs[0, :2])
        self._plot_spatial_categorical(
            column=cluster_col,
            ax=ax_spatial,
            title="Spatial Clusters",
            legend=True
        )
        
        # UMAP plot (if available)
        if 'X_umap' in self.adata.obsm:
            ax_umap = fig.add_subplot(gs[0, 2])
            sc.pl.umap(
                self.adata,
                color=cluster_col,
                ax=ax_umap,
                show=False,
                frameon=False,
                title='UMAP Clusters'
            )
        
        # Cluster size barplot
        ax_sizes = fig.add_subplot(gs[0, 3])
        cluster_sizes = self.adata.obs[cluster_col].value_counts().sort_index()
        cluster_sizes.plot(kind='bar', ax=ax_sizes, color='skyblue')
        ax_sizes.set_title('Cluster Sizes')
        ax_sizes.set_xlabel('Cluster')
        ax_sizes.set_ylabel('Number of Spots')
        ax_sizes.tick_params(axis='x', rotation=45)
        
        # QC metrics by cluster
        ax_qc1 = fig.add_subplot(gs[1, :2])
        if 'total_counts' in self.adata.obs.columns:
            sns.boxplot(
                data=self.adata.obs,
                x=cluster_col,
                y='total_counts',
                ax=ax_qc1
            )
            ax_qc1.set_title('Total UMI Counts by Cluster')
            ax_qc1.tick_params(axis='x', rotation=45)
        
        ax_qc2 = fig.add_subplot(gs[1, 2:])
        if 'n_genes_by_counts' in self.adata.obs.columns:
            sns.boxplot(
                data=self.adata.obs,
                x=cluster_col,
                y='n_genes_by_counts',
                ax=ax_qc2
            )
            ax_qc2.set_title('Genes Detected by Cluster')
            ax_qc2.tick_params(axis='x', rotation=45)
        
        # Marker genes heatmap (if available)
        if 'rank_genes_groups' in self.adata.uns:
            ax_markers = fig.add_subplot(gs[2, :])
            
            try:
                # Get top marker genes for each cluster
                n_genes = 3
                marker_genes = []
                clusters = self.adata.obs[cluster_col].unique()
                
                for cluster in sorted(clusters):
                    cluster_str = str(cluster)
                    if cluster_str in self.adata.uns['rank_genes_groups']['names'].dtype.names:
                        genes = self.adata.uns['rank_genes_groups']['names'][cluster_str][:n_genes]
                        marker_genes.extend(genes)
                
                # Remove duplicates while preserving order
                unique_markers = []
                for gene in marker_genes:
                    if gene not in unique_markers and gene in self.adata.var_names:
                        unique_markers.append(gene)
                
                if unique_markers:
                    # Create heatmap
                    sc.pl.heatmap(
                        self.adata,
                        unique_markers,
                        groupby=cluster_col,
                        ax=ax_markers,
                        show=False,
                        cmap='RdBu_r',
                        standard_scale='var'
                    )
                    ax_markers.set_title('Top Marker Genes by Cluster')
                else:
                    ax_markers.text(0.5, 0.5, 'No marker genes available', 
                                   ha='center', va='center', transform=ax_markers.transAxes)
                    ax_markers.set_title('Marker Genes')
                    
            except Exception as e:
                logger.warning(f"Could not create marker genes heatmap: {e}")
                ax_markers.text(0.5, 0.5, 'Marker genes heatmap failed', 
                               ha='center', va='center', transform=ax_markers.transAxes)
        
        plt.tight_layout()
        
        # Save figure
        save_path = self.save_dir / save_name
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        logger.info(f"Clustering summary saved to {save_path}")
        
        return fig
    
    def plot_qc_summary(self,
                       figsize: Tuple[int, int] = (16, 12),
                       save_name: str = "qc_summary.png") -> plt.Figure:
        """
        Create quality control summary figure.
        
        Parameters
        ----------
        figsize : tuple
            Figure size
        save_name : str
            Filename for saving
            
        Returns
        -------
        fig : matplotlib.Figure
            The created figure
        """
        logger.info("Creating QC summary figure")
        
        fig, axes = plt.subplots(3, 3, figsize=figsize)
        axes = axes.flatten()
        
        # QC metrics to plot
        qc_metrics = [
            ('total_counts', 'Total UMI Counts'),
            ('n_genes_by_counts', 'Genes Detected'),
            ('pct_counts_mt', 'Mitochondrial %')
        ]
        
        plot_idx = 0
        
        # Distribution plots
        for metric, title in qc_metrics:
            if metric in self.adata.obs.columns and plot_idx < len(axes):
                data = self.adata.obs[metric]
                
                # Histogram
                axes[plot_idx].hist(data, bins=50, alpha=0.7, edgecolor='black')
                axes[plot_idx].axvline(data.mean(), color='red', linestyle='--', 
                                      label=f'Mean: {data.mean():.1f}')
                axes[plot_idx].set_title(f'{title} Distribution')
                axes[plot_idx].set_xlabel(title)
                axes[plot_idx].set_ylabel('Frequency')
                axes[plot_idx].legend()
                plot_idx += 1
        
        # Spatial QC plots
        for metric, title in qc_metrics:
            if metric in self.adata.obs.columns and plot_idx < len(axes):
                self._plot_spatial_metric(
                    metric=metric,
                    ax=axes[plot_idx],
                    title=f'{title} (Spatial)',
                    spot_size=1.5
                )
                plot_idx += 1
        
        # Correlation plots
        if plot_idx < len(axes):
            # UMI vs Genes scatter
            if 'total_counts' in self.adata.obs and 'n_genes_by_counts' in self.adata.obs:
                axes[plot_idx].scatter(
                    self.adata.obs['total_counts'],
                    self.adata.obs['n_genes_by_counts'],
                    alpha=0.6, s=1
                )
                axes[plot_idx].set_xlabel('Total UMI Counts')
                axes[plot_idx].set_ylabel('Genes Detected')
                axes[plot_idx].set_title('UMI vs Genes Correlation')
                plot_idx += 1
        
        # Hide unused subplots
        for i in range(plot_idx, len(axes)):
            axes[i].set_visible(False)
        
        plt.tight_layout()
        
        # Save figure
        save_path = self.save_dir / save_name
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        logger.info(f"QC summary saved to {save_path}")
        
        return fig
    
    def _plot_spatial_metric(self,
                            metric: str,
                            ax: plt.Axes,
                            title: str,
                            cmap: str = 'viridis',
                            spot_size: float = 2.0) -> None:
        """Plot a continuous metric in spatial coordinates."""
        
        if 'spatial' not in self.adata.obsm:
            ax.text(0.5, 0.5, 'No spatial coordinates', 
                   ha='center', va='center', transform=ax.transAxes)
            ax.set_title(title)
            return
        
        coords = self.adata.obsm['spatial']
        values = self.adata.obs[metric]
        
        # Try to plot with tissue image background
        try:
            sc.pl.spatial(
                self.adata,
                color=metric,
                ax=ax,
                show=False,
                title=title,
                spot_size=spot_size,
                img_key='hires',
                alpha_img=0.7,
                cmap=cmap,
                frameon=False
            )
        except:
            # Fallback to simple scatter plot
            scatter = ax.scatter(
                coords[:, 0], coords[:, 1],
                c=values, cmap=cmap, s=spot_size, alpha=0.8
            )
            ax.set_title(title)
            ax.set_aspect('equal')
            ax.axis('off')
            
            # Add colorbar
            plt.colorbar(scatter, ax=ax, label=metric)
    
    def _plot_spatial_categorical(self,
                                 column: str,
                                 ax: plt.Axes,
                                 title: str,
                                 legend: bool = False) -> None:
        """Plot a categorical variable in spatial coordinates."""
        
        if 'spatial' not in self.adata.obsm:
            ax.text(0.5, 0.5, 'No spatial coordinates', 
                   ha='center', va='center', transform=ax.transAxes)
            ax.set_title(title)
            return
        
        try:
            sc.pl.spatial(
                self.adata,
                color=column,
                ax=ax,
                show=False,
                title=title,
                spot_size=2.0,
                img_key='hires',
                alpha_img=0.7,
                frameon=False,
                legend_loc='right margin' if legend else None
            )
        except:
            # Fallback to simple scatter plot
            coords = self.adata.obsm['spatial']
            categories = self.adata.obs[column]
            
            # Create color map for categories
            unique_cats = categories.unique()
            colors = plt.cm.tab20(np.linspace(0, 1, len(unique_cats)))
            color_map = dict(zip(unique_cats, colors))
            
            for cat in unique_cats:
                mask = categories == cat
                ax.scatter(
                    coords[mask, 0], coords[mask, 1],
                    c=[color_map[cat]], label=str(cat), s=2, alpha=0.8
                )
            
            ax.set_title(title)
            ax.set_aspect('equal')
            ax.axis('off')
            
            if legend:
                ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    
    def _plot_spatial_gene(self,
                          gene: str,
                          ax: plt.Axes,
                          spot_size: float = 2.0,
                          alpha_img: float = 0.7,
                          cmap: str = 'hot') -> None:
        """Plot gene expression in spatial coordinates."""
        
        if gene not in self.adata.var_names:
            ax.text(0.5, 0.5, f'Gene {gene} not found', 
                   ha='center', va='center', transform=ax.transAxes)
            ax.set_title(f"{gene} (Not Found)")
            return
        
        try:
            sc.pl.spatial(
                self.adata,
                color=gene,
                ax=ax,
                show=False,
                title=gene,
                spot_size=spot_size,
                img_key='hires',
                alpha_img=alpha_img,
                cmap=cmap,
                frameon=False
            )
        except:
            # Fallback to simple scatter plot
            coords = self.adata.obsm['spatial']
            expression = self.adata[:, gene].X
            if hasattr(expression, 'toarray'):
                expression = expression.toarray().flatten()
            
            scatter = ax.scatter(
                coords[:, 0], coords[:, 1],
                c=expression, cmap=cmap, s=spot_size, alpha=0.8
            )
            ax.set_title(gene)
            ax.set_aspect('equal')
            ax.axis('off')
            
            plt.colorbar(scatter, ax=ax, label='Expression')
    
    def create_publication_figure(self,
                                 figsize: Tuple[int, int] = (20, 16),
                                 save_name: str = "publication_figure.png") -> plt.Figure:
        """
        Create a comprehensive publication-ready figure.
        
        Parameters
        ----------
        figsize : tuple
            Figure size
        save_name : str
            Filename for saving
            
        Returns
        -------
        fig : matplotlib.Figure
            The created figure
        """
        logger.info("Creating publication-ready figure")
        
        fig = plt.figure(figsize=figsize)
        gs = fig.add_gridspec(4, 6, height_ratios=[1, 1, 1, 0.3])
        
        # Panel A: Raw tissue image and QC metrics
        ax_tissue = fig.add_subplot(gs[0, :2])
        if 'total_counts' in self.adata.obs.columns:
            self._plot_spatial_metric(
                metric='total_counts',
                ax=ax_tissue,
                title='A. Total UMI Counts',
                spot_size=1.5
            )
        
        # Panel B: Spatial clusters
        ax_clusters = fig.add_subplot(gs[0, 2:4])
        if 'spatial_clusters' in self.adata.obs.columns:
            self._plot_spatial_categorical(
                column='spatial_clusters',
                ax=ax_clusters,
                title='B. Spatial Clusters'
            )
        
        # Panel C: Key marker gene
        ax_marker = fig.add_subplot(gs[0, 4:])
        # Find a good marker gene to display
        marker_gene = None
        for gene in ['KRT19', 'COL1A1', 'CD68', 'EPCAM', 'VIM']:
            if gene in self.adata.var_names:
                marker_gene = gene
                break
        
        if marker_gene:
            self._plot_spatial_gene(
                gene=marker_gene,
                ax=ax_marker,
                title=f'C. {marker_gene} Expression',
                spot_size=1.5
            )
        
        # Panel D: Gene panel scores (if available)
        panel_cols = [col for col in self.adata.obs.columns if col.endswith('_score')]
        if panel_cols:
            n_panels = min(3, len(panel_cols))
            for i, col in enumerate(panel_cols[:n_panels]):
                ax_panel = fig.add_subplot(gs[1, i*2:(i+1)*2])
                panel_name = col.replace('_score', '').title()
                self._plot_spatial_metric(
                    metric=col,
                    ax=ax_panel,
                    title=f'D{i+1}. {panel_name} Score',
                    spot_size=1.5,
                    cmap='viridis'
                )
        
        # Panel E: Summary statistics
        ax_stats = fig.add_subplot(gs[2, :3])
        self._create_summary_table(ax_stats)
        
        # Panel F: Marker genes heatmap
        if 'rank_genes_groups' in self.adata.uns and 'spatial_clusters' in self.adata.obs.columns:
            ax_heatmap = fig.add_subplot(gs[2, 3:])
            self._create_marker_heatmap(ax_heatmap)
        
        # Add figure labels
        label_props = dict(fontsize=16, fontweight='bold')
        fig.text(0.02, 0.95, 'A', **label_props)
        fig.text(0.35, 0.95, 'B', **label_props)
        fig.text(0.68, 0.95, 'C', **label_props)
        
        plt.tight_layout()
        
        # Save figure
        save_path = self.save_dir / save_name
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        logger.info(f"Publication figure saved to {save_path}")
        
        return fig
    
    def _create_summary_table(self, ax: plt.Axes) -> None:
        """Create a summary statistics table."""
        
        stats_data = []
        stats_data.append(['Metric', 'Value'])
        stats_data.append(['Total Spots', f"{self.adata.n_obs:,}"])
        stats_data.append(['Total Genes', f"{self.adata.n_vars:,}"])
        
        if 'total_counts' in self.adata.obs.columns:
            mean_umi = self.adata.obs['total_counts'].mean()
            stats_data.append(['Mean UMI/spot', f"{mean_umi:,.0f}"])
        
        if 'spatial_clusters' in self.adata.obs.columns:
            n_clusters = len(self.adata.obs['spatial_clusters'].unique())
            stats_data.append(['Spatial Clusters', str(n_clusters)])
        
        # Create table
        table = ax.table(
            cellText=stats_data[1:],
            colLabels=stats_data[0],
            cellLoc='center',
            loc='center',
            colWidths=[0.5, 0.5]
        )
        table.auto_set_font_size(False)
        table.set_fontsize(10)
        table.scale(1, 2)
        
        ax.set_title('E. Dataset Summary')
        ax.axis('off')
    
    def _create_marker_heatmap(self, ax: plt.Axes) -> None:
        """Create marker genes heatmap."""
        
        try:
            # Get top 2 marker genes per cluster
            marker_genes = []
            clusters = self.adata.obs['spatial_clusters'].unique()
            
            for cluster in sorted(clusters):
                cluster_str = str(cluster)
                if cluster_str in self.adata.uns['rank_genes_groups']['names'].dtype.names:
                    genes = self.adata.uns['rank_genes_groups']['names'][cluster_str][:2]
                    marker_genes.extend(genes)
            
            # Remove duplicates
            unique_markers = []
            for gene in marker_genes:
                if gene not in unique_markers and gene in self.adata.var_names:
                    unique_markers.append(gene)
            
            if unique_markers:
                sc.pl.heatmap(
                    self.adata,
                    unique_markers[:10],  # Limit to 10 genes
                    groupby='spatial_clusters',
                    ax=ax,
                    show=False,
                    cmap='RdBu_r',
                    standard_scale='var'
                )
                ax.set_title('F. Cluster Marker Genes')
            else:
                ax.text(0.5, 0.5, 'No marker genes available', 
                       ha='center', va='center', transform=ax.transAxes)
                ax.set_title('F. Marker Genes')
                
        except Exception as e:
            ax.text(0.5, 0.5, f'Heatmap error: {str(e)[:50]}...', 
                   ha='center', va='center', transform=ax.transAxes)
            ax.set_title('F. Marker Genes (Error)')


def create_comprehensive_visualization(adata: ad.AnnData,
                                     gene_panels: Optional[Dict[str, List[str]]] = None,
                                     save_dir: str = "./figures") -> Dict[str, plt.Figure]:
    """
    Create comprehensive visualization suite.
    
    Parameters
    ----------
    adata : AnnData
        Analyzed spatial transcriptomics data
    gene_panels : dict, optional
        Gene panels for visualization
    save_dir : str
        Directory to save figures
        
    Returns
    -------
    figures : dict
        Dictionary of created figures
    """
    logger.info("Creating comprehensive visualization suite")
    
    # Initialize visualizer
    visualizer = SpatialVisualizer(adata, save_dir=save_dir)
    
    figures = {}
    
    # QC summary
    try:
        figures['qc_summary'] = visualizer.plot_qc_summary()
    except Exception as e:
        logger.warning(f"QC summary failed: {e}")
    
    # Spatial overview
    try:
        figures['spatial_overview'] = visualizer.plot_spatial_overview()
    except Exception as e:
        logger.warning(f"Spatial overview failed: {e}")
    
    # Clustering summary
    if 'spatial_clusters' in adata.obs.columns:
        try:
            figures['clustering_summary'] = visualizer.plot_clustering_summary()
        except Exception as e:
            logger.warning(f"Clustering summary failed: {e}")
    
    # Gene panels
    if gene_panels:
        try:
            figures['gene_panels'] = visualizer.plot_gene_panels_summary(gene_panels)
        except Exception as e:
            logger.warning(f"Gene panels visualization failed: {e}")
        
        # Individual gene expression panels
        for panel_name, genes in gene_panels.items():
            try:
                figures[f'genes_{panel_name}'] = visualizer.plot_gene_expression_panel(
                    genes=genes,
                    save_name=f"gene_panel_{panel_name}.png"
                )
            except Exception as e:
                logger.warning(f"Gene panel {panel_name} failed: {e}")
    
    # Publication figure
    try:
        figures['publication'] = visualizer.create_publication_figure()
    except Exception as e:
        logger.warning(f"Publication figure failed: {e}")
    
    logger.info(f"Created {len(figures)} visualization figures")
    return figures


if __name__ == "__main__":
    # Example usage
    import sys
    sys.path.append('..')
    from modules.data_loader import load_and_validate_visium
    from modules.quality_control import run_comprehensive_qc
    from modules.preprocessing import run_standard_preprocessing
    from modules.spatial_analysis import run_spatial_analysis_pipeline
    
    try:
        # Load and process data
        adata, _ = load_and_validate_visium("../data/example_data/")
        adata_qc, _ = run_comprehensive_qc(adata, plot_results=False)
        adata_processed, _ = run_standard_preprocessing(adata_qc, plot_results=False)
        
        # Gene panels
        gene_panels = {
            'epithelial': ['KRT19', 'KRT8', 'EPCAM'],
            'stromal': ['COL1A1', 'VIM', 'ACTA2'],
            'immune': ['CD68', 'CD3E', 'PTPRC']
        }
        
        # Spatial analysis
        adata_spatial, _ = run_spatial_analysis_pipeline(
            adata_processed,
            gene_panels=gene_panels
        )
        
        # Create visualizations
        figures = create_comprehensive_visualization(
            adata_spatial,
            gene_panels=gene_panels,
            save_dir="../results/figures/"
        )
        
        print(f"Created {len(figures)} visualization figures")
        
    except Exception as e:
        print(f"Error in visualization: {e}")

        