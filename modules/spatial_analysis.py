"""
Spatial Analysis Module

This module provides spatial pattern analysis, clustering, and neighborhood
analysis for spatial transcriptomics data.

Author: [Your Name]
Date: 2025
Project: Spatial Transcriptomics Analysis Pipeline
"""

import scanpy as sc
import squidpy as sq
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import anndata as ad
from typing import Optional, Dict, List, Tuple, Any, Union
import logging
from sklearn.metrics import adjusted_rand_score, silhouette_score
from scipy.spatial.distance import pdist, squareform
from scipy import stats

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class SpatialAnalyzer:
    """
    Comprehensive spatial analysis for spatial transcriptomics data.
    
    This class implements spatial clustering, neighborhood analysis,
    marker gene detection, and spatial pattern characterization.
    """
    
    def __init__(self, adata: ad.AnnData):
        """
        Initialize spatial analyzer.
        
        Parameters
        ----------
        adata : AnnData
            Preprocessed spatial transcriptomics data
        """
        self.adata = adata.copy()
        self.spatial_results = {}
        
        # Validate spatial data
        if 'spatial' not in self.adata.obsm:
            raise ValueError("No spatial coordinates found. Ensure data has 'spatial' in .obsm")
    
    def build_spatial_graph(self,
                           n_neighbors: int = 10,
                           n_rings: int = 1,
                           coord_type: str = 'generic',
                           radius: Optional[float] = None) -> ad.AnnData:
        """
        Build spatial neighborhood graph.
        
        Parameters
        ----------
        n_neighbors : int
            Number of spatial neighbors
        n_rings : int
            Number of rings of neighbors to include
        coord_type : str
            Type of coordinates ('visium', 'generic')
        radius : float, optional
            Radius for radius-based neighbors
            
        Returns
        -------
        adata : AnnData
            Data with spatial graph added
        """
        logger.info(f"Building spatial graph with {n_neighbors} neighbors")
        
        # Use squidpy for advanced spatial graph construction
        try:
            if coord_type == 'visium':
                sq.gr.spatial_neighbors(
                    self.adata,
                    coord_type='visium',
                    n_rings=n_rings
                )
            else:
                sq.gr.spatial_neighbors(
                    self.adata,
                    coord_type='generic',
                    n_neighs=n_neighbors,
                    radius=radius
                )
        except ImportError:
            logger.warning("Squidpy not available, using scanpy for graph construction")
            # Fallback to scanpy
            sc.pp.neighbors(
                self.adata,
                n_neighbors=n_neighbors,
                use_rep='spatial',
                metric='euclidean'
            )
        
        # Store graph statistics
        if 'spatial_connectivities' in self.adata.obsp:
            connectivity_matrix = self.adata.obsp['spatial_connectivities']
        elif 'connectivities' in self.adata.obsp:
            connectivity_matrix = self.adata.obsp['connectivities']
        else:
            logger.warning("No connectivity matrix found")
            return self.adata
        
        # Calculate graph statistics
        n_connections = np.array(connectivity_matrix.sum(axis=1)).flatten()
        self.spatial_results['graph_stats'] = {
            'mean_connections': n_connections.mean(),
            'std_connections': n_connections.std(),
            'min_connections': n_connections.min(),
            'max_connections': n_connections.max()
        }
        
        logger.info(f"Spatial graph built: {n_connections.mean():.1f} ± {n_connections.std():.1f} connections per spot")
        return self.adata
    
    def spatial_clustering(self,
                          method: str = 'leiden',
                          resolution: float = 0.5,
                          n_iterations: int = 2,
                          random_state: int = 42) -> ad.AnnData:
        """
        Perform spatial clustering.
        
        Parameters
        ----------
        method : str
            Clustering method ('leiden', 'louvain')
        resolution : float
            Clustering resolution
        n_iterations : int
            Number of clustering iterations
        random_state : int
            Random seed for reproducibility
            
        Returns
        -------
        adata : AnnData
            Data with clustering results added
        """
        logger.info(f"Performing {method} clustering with resolution {resolution}")
        
        # Check if spatial graph exists
        if ('spatial_connectivities' not in self.adata.obsp and 
            'connectivities' not in self.adata.obsp):
            logger.warning("No spatial graph found, building default graph")
            self.build_spatial_graph()
        
        # Perform clustering
        if method == 'leiden':
            sc.tl.leiden(
                self.adata,
                resolution=resolution,
                n_iterations=n_iterations,
                random_state=random_state,
                key_added='spatial_clusters'
            )
        elif method == 'louvain':
            sc.tl.louvain(
                self.adata,
                resolution=resolution,
                random_state=random_state,
                key_added='spatial_clusters'
            )
        else:
            raise ValueError(f"Unknown clustering method: {method}")
        
        # Calculate clustering statistics
        clusters = self.adata.obs['spatial_clusters'].astype(str)
        n_clusters = len(clusters.unique())
        cluster_sizes = clusters.value_counts().sort_index()
        
        # Calculate silhouette score if possible
        try:
            if 'X_pca' in self.adata.obsm:
                sil_score = silhouette_score(
                    self.adata.obsm['X_pca'][:, :10],
                    clusters
                )
            else:
                # Use spatial coordinates
                sil_score = silhouette_score(
                    self.adata.obsm['spatial'],
                    clusters
                )
        except:
            sil_score = np.nan
        
        self.spatial_results['clustering'] = {
            'method': method,
            'resolution': resolution,
            'n_clusters': n_clusters,
            'cluster_sizes': cluster_sizes.to_dict(),
            'silhouette_score': sil_score
        }
        
        logger.info(f"Clustering complete: {n_clusters} clusters found")
        logger.info(f"Cluster sizes: {dict(cluster_sizes)}")
        
        return self.adata
    
    def find_marker_genes(self,
                         groupby: str = 'spatial_clusters',
                         method: str = 'wilcoxon',
                         n_genes: int = 10,
                         min_fold_change: float = 1.5) -> ad.AnnData:
        """
        Find marker genes for spatial clusters.
        
        Parameters
        ----------
        groupby : str
            Column to group by for marker detection
        method : str
            Statistical method ('wilcoxon', 't-test', 'logreg')
        n_genes : int
            Number of top genes per group
        min_fold_change : float
            Minimum fold change threshold
            
        Returns
        -------
        adata : AnnData
            Data with marker gene results added
        """
        logger.info(f"Finding marker genes for {groupby}")
        
        if groupby not in self.adata.obs.columns:
            raise ValueError(f"Column {groupby} not found in .obs")
        
        # Perform marker gene analysis
        sc.tl.rank_genes_groups(
            self.adata,
            groupby=groupby,
            method=method,
            n_genes=n_genes,
            use_raw=True if self.adata.raw is not None else False
        )
        
        # Extract and organize results
        marker_results = {}
        groups = self.adata.obs[groupby].unique()
        
        for group in groups:
            group_str = str(group)
            
            # Get top genes for this group
            gene_names = self.adata.uns['rank_genes_groups']['names'][group_str][:n_genes]
            scores = self.adata.uns['rank_genes_groups']['scores'][group_str][:n_genes]
            pvals = self.adata.uns['rank_genes_groups']['pvals'][group_str][:n_genes]
            logfold = self.adata.uns['rank_genes_groups']['logfoldchanges'][group_str][:n_genes]
            
            # Filter by fold change
            valid_genes = []
            for i in range(len(gene_names)):
                if abs(logfold[i]) >= np.log2(min_fold_change):
                    valid_genes.append({
                        'gene': gene_names[i],
                        'score': scores[i],
                        'pval': pvals[i],
                        'log2fc': logfold[i],
                        'fold_change': 2**logfold[i]
                    })
            
            marker_results[group_str] = valid_genes
        
        self.spatial_results['marker_genes'] = marker_results
        
        # Log summary
        total_markers = sum(len(genes) for genes in marker_results.values())
        logger.info(f"Found {total_markers} marker genes across {len(groups)} groups")
        
        return self.adata
    
    def spatial_autocorrelation(self,
                               genes: Optional[List[str]] = None,
                               n_perms: int = 1000) -> Dict[str, float]:
        """
        Calculate spatial autocorrelation (Moran's I) for genes.
        
        Parameters
        ----------
        genes : list of str, optional
            Genes to analyze. If None, uses highly variable genes
        n_perms : int
            Number of permutations for significance testing
            
        Returns
        -------
        moran_results : dict
            Moran's I statistics for each gene
        """
        logger.info("Calculating spatial autocorrelation")
        
        # Select genes to analyze
        if genes is None:
            if 'highly_variable' in self.adata.var.columns:
                genes = self.adata.var_names[self.adata.var.highly_variable][:100]
            else:
                # Use top variable genes
                gene_vars = np.var(self.adata.X.toarray() if hasattr(self.adata.X, 'toarray') else self.adata.X, axis=0)
                top_var_idx = np.argsort(gene_vars)[-100:]
                genes = self.adata.var_names[top_var_idx]
        
        genes = [g for g in genes if g in self.adata.var_names]
        logger.info(f"Analyzing spatial autocorrelation for {len(genes)} genes")
        
        # Try using squidpy for Moran's I
        try:
            sq.gr.spatial_autocorr(
                self.adata,
                genes=genes,
                mode='moran',
                n_perms=n_perms
            )
            
            moran_results = self.adata.uns['moranI'].copy()
            
        except ImportError:
            logger.warning("Squidpy not available, using manual Moran's I calculation")
            moran_results = self._calculate_moran_manual(genes)
        
        # Store results
        self.spatial_results['spatial_autocorr'] = moran_results
        
        # Log top spatially autocorrelated genes
        if isinstance(moran_results, pd.DataFrame):
            top_genes = moran_results.nlargest(5, 'I')['I']
            logger.info(f"Top spatially autocorrelated genes: {dict(top_genes)}")
        
        return moran_results
    
    def _calculate_moran_manual(self, genes: List[str]) -> pd.DataFrame:
        """Manual calculation of Moran's I when squidpy is not available."""
        
        # Get spatial coordinates
        coords = self.adata.obsm['spatial']
        
        # Calculate distance matrix
        distances = pdist(coords, metric='euclidean')
        dist_matrix = squareform(distances)
        
        # Create spatial weights (inverse distance)
        weights = 1.0 / (dist_matrix + 1e-10)  # Add small value to avoid division by zero
        np.fill_diagonal(weights, 0)  # No self-weight
        
        # Row-normalize weights
        row_sums = weights.sum(axis=1)
        weights = weights / row_sums[:, np.newaxis]
        
        results = []
        for gene in genes:
            if gene in self.adata.var_names:
                # Get expression values
                expr = self.adata[:, gene].X
                if hasattr(expr, 'toarray'):
                    expr = expr.toarray().flatten()
                
                # Calculate Moran's I
                n = len(expr)
                mean_expr = np.mean(expr)
                
                # Numerator: sum of weighted cross-products
                numerator = 0
                for i in range(n):
                    for j in range(n):
                        numerator += weights[i, j] * (expr[i] - mean_expr) * (expr[j] - mean_expr)
                
                # Denominator: sum of squared deviations
                denominator = np.sum((expr - mean_expr) ** 2)
                
                # Moran's I
                if denominator > 0:
                    moran_i = numerator / denominator
                else:
                    moran_i = 0
                
                results.append({'gene': gene, 'I': moran_i})
        
        return pd.DataFrame(results).set_index('gene')
    
    def analyze_gene_panels(self,
                           gene_panels: Dict[str, List[str]]) -> Dict[str, Any]:
        """
        Analyze expression patterns of predefined gene panels.
        
        Parameters
        ----------
        gene_panels : dict
            Dictionary mapping panel names to gene lists
            
        Returns
        -------
        panel_results : dict
            Analysis results for each gene panel
        """
        logger.info(f"Analyzing {len(gene_panels)} gene panels")
        
        panel_results = {}
        
        for panel_name, genes in gene_panels.items():
            # Filter to available genes
            available_genes = [g for g in genes if g in self.adata.var_names]
            
            if not available_genes:
                logger.warning(f"No genes found for panel {panel_name}")
                continue
            
            logger.info(f"Panel {panel_name}: {len(available_genes)}/{len(genes)} genes available")
            
            # Calculate panel score (mean expression of available genes)
            panel_expr = self.adata[:, available_genes].X
            if hasattr(panel_expr, 'toarray'):
                panel_expr = panel_expr.toarray()
            
            panel_score = np.mean(panel_expr, axis=1)
            
            # Add to observations
            self.adata.obs[f'{panel_name}_score'] = panel_score
            
            # Calculate statistics
            panel_stats = {
                'genes_available': available_genes,
                'genes_missing': [g for g in genes if g not in self.adata.var_names],
                'mean_score': panel_score.mean(),
                'std_score': panel_score.std(),
                'min_score': panel_score.min(),
                'max_score': panel_score.max()
            }
            
            # Spatial autocorrelation of panel score
            try:
                coords = self.adata.obsm['spatial']
                # Simple spatial correlation (correlation between score and spatial coordinates)
                spatial_corr_x = stats.pearsonr(coords[:, 0], panel_score)[0]
                spatial_corr_y = stats.pearsonr(coords[:, 1], panel_score)[0]
                panel_stats['spatial_corr'] = {'x': spatial_corr_x, 'y': spatial_corr_y}
            except:
                panel_stats['spatial_corr'] = {'x': np.nan, 'y': np.nan}
            
            panel_results[panel_name] = panel_stats
        
        self.spatial_results['gene_panels'] = panel_results
        return panel_results
    
    def spatial_domains(self,
                       cluster_key: str = 'spatial_clusters') -> Dict[str, Any]:
        """
        Analyze spatial organization of clusters/domains.
        
        Parameters
        ----------
        cluster_key : str
            Column containing cluster assignments
            
        Returns
        -------
        domain_stats : dict
            Spatial domain analysis results
        """
        logger.info("Analyzing spatial domains")
        
        if cluster_key not in self.adata.obs.columns:
            raise ValueError(f"Column {cluster_key} not found")
        
        coords = self.adata.obsm['spatial']
        clusters = self.adata.obs[cluster_key].astype(str)
        
        domain_stats = {}
        
        for cluster in clusters.unique():
            cluster_mask = clusters == cluster
            cluster_coords = coords[cluster_mask]
            
            if len(cluster_coords) < 2:
                continue
            
            # Calculate cluster spatial properties
            centroid = np.mean(cluster_coords, axis=0)
            
            # Compactness (mean distance to centroid)
            distances_to_centroid = np.sqrt(np.sum((cluster_coords - centroid)**2, axis=1))
            compactness = np.mean(distances_to_centroid)
            
            # Convex hull area (if possible)
            try:
                from scipy.spatial import ConvexHull
                if len(cluster_coords) >= 3:
                    hull = ConvexHull(cluster_coords)
                    convex_area = hull.volume  # In 2D, volume is area
                else:
                    convex_area = 0
            except:
                convex_area = np.nan
            
            # Nearest neighbor analysis within cluster
            if len(cluster_coords) > 1:
                cluster_distances = pdist(cluster_coords, metric='euclidean')
                mean_nn_dist = np.mean(cluster_distances)
            else:
                mean_nn_dist = np.nan
            
            domain_stats[cluster] = {
                'n_spots': int(cluster_mask.sum()),
                'centroid': centroid.tolist(),
                'compactness': compactness,
                'convex_area': convex_area,
                'mean_nn_distance': mean_nn_dist
            }
        
        self.spatial_results['spatial_domains'] = domain_stats
        return domain_stats
    
    def get_analysis_summary(self) -> str:
        """
        Get a summary of spatial analysis results.
        
        Returns
        -------
        summary : str
            Formatted analysis summary
        """
        if not self.spatial_results:
            return "No spatial analysis performed yet."
        
        summary = []
        summary.append("=== SPATIAL ANALYSIS SUMMARY ===")
        summary.append(f"Data shape: {self.adata.n_obs} spots × {self.adata.n_vars} genes")
        
        # Graph statistics
        if 'graph_stats' in self.spatial_results:
            graph = self.spatial_results['graph_stats']
            summary.append(f"\n🕸️ SPATIAL GRAPH:")
            summary.append(f"   Mean connections: {graph['mean_connections']:.1f} ± {graph['std_connections']:.1f}")
        
        # Clustering results
        if 'clustering' in self.spatial_results:
            clust = self.spatial_results['clustering']
            summary.append(f"\n🎯 CLUSTERING:")
            summary.append(f"   Method: {clust['method']} (resolution: {clust['resolution']})")
            summary.append(f"   Clusters found: {clust['n_clusters']}")
            if not np.isnan(clust['silhouette_score']):
                summary.append(f"   Silhouette score: {clust['silhouette_score']:.3f}")
        
        # Marker genes
        if 'marker_genes' in self.spatial_results:
            markers = self.spatial_results['marker_genes']
            total_markers = sum(len(genes) for genes in markers.values())
            summary.append(f"\n🧬 MARKER GENES:")
            summary.append(f"   Total markers found: {total_markers}")
            summary.append(f"   Clusters with markers: {len(markers)}")
        
        # Gene panels
        if 'gene_panels' in self.spatial_results:
            panels = self.spatial_results['gene_panels']
            summary.append(f"\n📋 GENE PANELS:")
            summary.append(f"   Panels analyzed: {len(panels)}")
            for name, stats in panels.items():
                summary.append(f"   {name}: {len(stats['genes_available'])} genes")
        
        # Spatial domains
        if 'spatial_domains' in self.spatial_results:
            domains = self.spatial_results['spatial_domains']
            summary.append(f"\n🗺️ SPATIAL DOMAINS:")
            summary.append(f"   Domains analyzed: {len(domains)}")
            sizes = [stats['n_spots'] for stats in domains.values()]
            summary.append(f"   Size range: {min(sizes)} - {max(sizes)} spots")
        
        return "\n".join(summary)


def run_spatial_analysis_pipeline(adata: ad.AnnData,
                                 resolution: float = 0.5,
                                 n_neighbors: int = 10,
                                 gene_panels: Optional[Dict[str, List[str]]] = None,
                                 find_markers: bool = True) -> Tuple[ad.AnnData, SpatialAnalyzer]:
    """
    Run complete spatial analysis pipeline.
    
    Parameters
    ----------
    adata : AnnData
        Preprocessed spatial transcriptomics data
    resolution : float
        Clustering resolution
    n_neighbors : int
        Number of spatial neighbors
    gene_panels : dict, optional
        Gene panels to analyze
    find_markers : bool
        Whether to find marker genes
        
    Returns
    -------
    adata_spatial : AnnData
        Data with spatial analysis results
    analyzer : SpatialAnalyzer
        Analyzer object with results
    """
    logger.info("Running spatial analysis pipeline")
    
    # Initialize analyzer
    analyzer = SpatialAnalyzer(adata)
    
    # Build spatial graph
    adata_spatial = analyzer.build_spatial_graph(n_neighbors=n_neighbors)
    
    # Perform clustering
    adata_spatial = analyzer.spatial_clustering(resolution=resolution)
    
    # Find marker genes
    if find_markers:
        adata_spatial = analyzer.find_marker_genes()
    
    # Analyze gene panels if provided
    if gene_panels:
        analyzer.analyze_gene_panels(gene_panels)
    
    # Spatial domain analysis
    analyzer.spatial_domains()
    
    # Calculate spatial autocorrelation for top genes
    try:
        analyzer.spatial_autocorrelation()
    except Exception as e:
        logger.warning(f"Spatial autocorrelation failed: {e}")
    
    logger.info("Spatial analysis pipeline complete")
    return adata_spatial, analyzer


if __name__ == "__main__":
    # Example usage
    import sys
    sys.path.append('..')
    from modules.data_loader import load_and_validate_visium
    from modules.quality_control import run_comprehensive_qc
    from modules.preprocessing import run_standard_preprocessing
    
    try:
        # Load and process data
        adata, _ = load_and_validate_visium("../data/example_data/")
        adata_qc, _ = run_comprehensive_qc(adata, plot_results=False)
        adata_processed, _ = run_standard_preprocessing(adata_qc, plot_results=False)
        
        # Define gene panels
        gene_panels = {
            'epithelial': ['KRT19', 'KRT8', 'EPCAM'],
            'stromal': ['COL1A1', 'VIM', 'ACTA2'],
            'immune': ['CD68', 'CD3E', 'PTPRC']
        }
        
        # Run spatial analysis
        adata_spatial, analyzer = run_spatial_analysis_pipeline(
            adata_processed,
            gene_panels=gene_panels
        )
        
        # Print summary
        print(analyzer.get_analysis_summary())
        
    except Exception as e:
        print(f"Error in spatial analysis: {e}")
