"""
Utility Functions Module

This module provides helper functions, constants, and utilities used
across the spatial transcriptomics analysis pipeline.

Author: [Your Name]
Date: 2025
Project: Spatial Transcriptomics Analysis Pipeline
"""

import pandas as pd
import numpy as np
import anndata as ad
import scanpy as sc
from typing import Dict, List, Tuple, Any, Optional, Union
import logging
import json
import yaml
from pathlib import Path
import pickle
import warnings

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


# =====================================
# CONSTANTS AND CONFIGURATION
# =====================================

# Default gene panels for different biological contexts
DEFAULT_GENE_PANELS = {
    'breast_cancer': {
        'epithelial': ['KRT19', 'KRT8', 'KRT18', 'EPCAM', 'CDH1', 'CLDN3', 'CLDN4'],
        'stromal': ['COL1A1', 'COL1A2', 'VIM', 'ACTA2', 'PDGFRA', 'FAP', 'THY1'],
        'immune': ['CD68', 'CD3E', 'CD8A', 'CD4', 'FOXP3', 'PTPRC', 'CD19', 'CD20'],
        'endothelial': ['PECAM1', 'VWF', 'CD34', 'CLDN5', 'FLT1'],
        'cancer_markers': ['ESR1', 'PGR', 'ERBB2', 'MKI67', 'TP53', 'CCND1'],
        'housekeeping': ['ACTB', 'GAPDH', 'B2M', 'PPIA', 'HPRT1'],
        'mitochondrial': ['MT-CO1', 'MT-ND1', 'MT-ATP6', 'MT-CYB', 'MT-ND2']
    },
    'brain': {
        'neurons': ['RBFOX3', 'SYP', 'MAP2', 'TUBB3', 'ENO2'],
        'astrocytes': ['GFAP', 'S100B', 'AQP4', 'ALDH1L1'],
        'oligodendrocytes': ['MBP', 'PLP1', 'MOG', 'OLIG2'],
        'microglia': ['IBA1', 'CX3CR1', 'TMEM119', 'P2RY12'],
        'endothelial': ['PECAM1', 'VWF', 'CD34', 'CLDN5']
    },
    'immune': {
        't_cells': ['CD3E', 'CD3D', 'CD8A', 'CD4', 'IL7R'],
        'b_cells': ['CD19', 'CD20', 'MS4A1', 'CD79A'],
        'macrophages': ['CD68', 'CD163', 'CSF1R', 'MARCO'],
        'dendritic': ['CD1C', 'CLEC9A', 'FCER1A', 'CD141'],
        'nk_cells': ['KLRD1', 'KLRF1', 'NCR1', 'GNLY']
    }
}

# Quality control thresholds
DEFAULT_QC_THRESHOLDS = {
    'min_counts': 1000,
    'max_counts': 50000,
    'min_genes': 500,
    'max_genes': 8000,
    'max_pct_mt': 20.0,
    'max_pct_ribo': 50.0,
    'min_cells_per_gene': 10
}

# Preprocessing parameters
DEFAULT_PREPROCESSING_PARAMS = {
    'target_sum': 10000,
    'n_top_genes': 2000,
    'hvg_method': 'seurat_v3',
    'max_value': 10,
    'n_pcs': 50
}

# Spatial analysis parameters
DEFAULT_SPATIAL_PARAMS = {
    'n_neighbors': 10,
    'resolution': 0.5,
    'n_iterations': 2,
    'coord_type': 'generic'
}

# Color palettes for visualization
COLOR_PALETTES = {
    'clusters': ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', 
                '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf'],
    'expression': ['hot', 'viridis', 'plasma', 'inferno', 'magma'],
    'diverging': ['RdBu_r', 'RdYlBu_r', 'coolwarm', 'seismic']
}


# =====================================
# DATA HANDLING UTILITIES
# =====================================

def save_adata(adata: ad.AnnData, 
               filepath: Union[str, Path],
               compression: str = 'gzip') -> None:
    """
    Save AnnData object with error handling.
    
    Parameters
    ----------
    adata : AnnData
        Data to save
    filepath : str or Path
        Output file path
    compression : str
        Compression method
    """
    try:
        filepath = Path(filepath)
        filepath.parent.mkdir(parents=True, exist_ok=True)
        adata.write(filepath, compression=compression)
        logger.info(f"Data saved to {filepath}")
    except Exception as e:
        logger.error(f"Failed to save data: {e}")
        raise


def load_adata(filepath: Union[str, Path]) -> ad.AnnData:
    """
    Load AnnData object with error handling.
    
    Parameters
    ----------
    filepath : str or Path
        Input file path
        
    Returns
    -------
    adata : AnnData
        Loaded data
    """
    try:
        adata = sc.read_h5ad(filepath)
        logger.info(f"Data loaded from {filepath}")
        return adata
    except Exception as e:
        logger.error(f"Failed to load data: {e}")
        raise


def backup_adata(adata: ad.AnnData, 
                backup_dir: Union[str, Path] = "./backups",
                prefix: str = "adata_backup") -> Path:
    """
    Create a backup of AnnData object.
    
    Parameters
    ----------
    adata : AnnData
        Data to backup
    backup_dir : str or Path
        Backup directory
    prefix : str
        Filename prefix
        
    Returns
    -------
    backup_path : Path
        Path to backup file
    """
    import datetime
    
    backup_dir = Path(backup_dir)
    backup_dir.mkdir(parents=True, exist_ok=True)
    
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    backup_path = backup_dir / f"{prefix}_{timestamp}.h5ad"
    
    save_adata(adata, backup_path)
    return backup_path


# =====================================
# CONFIGURATION UTILITIES
# =====================================

def load_config(config_path: Union[str, Path]) -> Dict[str, Any]:
    """
    Load configuration from YAML or JSON file.
    
    Parameters
    ----------
    config_path : str or Path
        Path to configuration file
        
    Returns
    -------
    config : dict
        Configuration dictionary
    """
    config_path = Path(config_path)
    
    if not config_path.exists():
        raise FileNotFoundError(f"Configuration file not found: {config_path}")
    
    try:
        with open(config_path, 'r') as f:
            if config_path.suffix.lower() in ['.yaml', '.yml']:
                config = yaml.safe_load(f)
            elif config_path.suffix.lower() == '.json':
                config = json.load(f)
            else:
                raise ValueError(f"Unsupported config format: {config_path.suffix}")
        
        logger.info(f"Configuration loaded from {config_path}")
        return config
        
    except Exception as e:
        logger.error(f"Failed to load configuration: {e}")
        raise


def save_config(config: Dict[str, Any], 
               config_path: Union[str, Path]) -> None:
    """
    Save configuration to YAML or JSON file.
    
    Parameters
    ----------
    config : dict
        Configuration dictionary
    config_path : str or Path
        Output file path
    """
    config_path = Path(config_path)
    config_path.parent.mkdir(parents=True, exist_ok=True)
    
    try:
        with open(config_path, 'w') as f:
            if config_path.suffix.lower() in ['.yaml', '.yml']:
                yaml.dump(config, f, default_flow_style=False, indent=2)
            elif config_path.suffix.lower() == '.json':
                json.dump(config, f, indent=2)
            else:
                raise ValueError(f"Unsupported config format: {config_path.suffix}")
        
        logger.info(f"Configuration saved to {config_path}")
        
    except Exception as e:
        logger.error(f"Failed to save configuration: {e}")
        raise


def merge_configs(*configs: Dict[str, Any]) -> Dict[str, Any]:
    """
    Merge multiple configuration dictionaries.
    
    Parameters
    ----------
    *configs : dict
        Configuration dictionaries to merge
        
    Returns
    -------
    merged_config : dict
        Merged configuration
    """
    merged = {}
    for config in configs:
        if isinstance(config, dict):
            for key, value in config.items():
                if key in merged and isinstance(merged[key], dict) and isinstance(value, dict):
                    merged[key] = merge_configs(merged[key], value)
                else:
                    merged[key] = value
    return merged


# =====================================
# GENE PANEL UTILITIES
# =====================================

def get_gene_panel(panel_name: str, 
                  tissue_type: str = 'breast_cancer') -> List[str]:
    """
    Get predefined gene panel.
    
    Parameters
    ----------
    panel_name : str
        Name of the gene panel
    tissue_type : str
        Tissue type context
        
    Returns
    -------
    genes : list of str
        List of genes in the panel
    """
    if tissue_type not in DEFAULT_GENE_PANELS:
        raise ValueError(f"Unknown tissue type: {tissue_type}")
    
    if panel_name not in DEFAULT_GENE_PANELS[tissue_type]:
        raise ValueError(f"Unknown panel: {panel_name} for tissue {tissue_type}")
    
    return DEFAULT_GENE_PANELS[tissue_type][panel_name].copy()


def filter_gene_panel(genes: List[str], 
                     adata: ad.AnnData,
                     min_expression: float = 0.1,
                     min_cells: int = 10) -> Tuple[List[str], List[str]]:
    """
    Filter gene panel to available and well-expressed genes.
    
    Parameters
    ----------
    genes : list of str
        Input gene list
    adata : AnnData
        Data to check against
    min_expression : float
        Minimum mean expression
    min_cells : int
        Minimum number of expressing cells
        
    Returns
    -------
    available_genes : list of str
        Available and well-expressed genes
    missing_genes : list of str
        Missing or poorly expressed genes
    """
    available_genes = []
    missing_genes = []
    
    for gene in genes:
        if gene not in adata.var_names:
            missing_genes.append(gene)
            continue
        
        # Check expression levels
        expr = adata[:, gene].X
        if hasattr(expr, 'toarray'):
            expr = expr.toarray().flatten()
        
        mean_expr = expr.mean()
        n_expressing = (expr > 0).sum()
        
        if mean_expr >= min_expression and n_expressing >= min_cells:
            available_genes.append(gene)
        else:
            missing_genes.append(gene)
            logger.debug(f"Gene {gene} filtered: mean_expr={mean_expr:.3f}, n_cells={n_expressing}")
    
    return available_genes, missing_genes


def create_custom_gene_panel(adata: ad.AnnData,
                            cluster_key: str,
                            n_genes_per_cluster: int = 5,
                            min_fold_change: float = 2.0) -> Dict[str, List[str]]:
    """
    Create custom gene panels based on cluster markers.
    
    Parameters
    ----------
    adata : AnnData
        Data with cluster annotations
    cluster_key : str
        Column containing cluster assignments
    n_genes_per_cluster : int
        Number of genes per cluster
    min_fold_change : float
        Minimum fold change threshold
        
    Returns
    -------
    gene_panels : dict
        Dictionary mapping cluster names to gene lists
    """
    if 'rank_genes_groups' not in adata.uns:
        raise ValueError("No marker genes found. Run marker gene analysis first.")
    
    gene_panels = {}
    
    # Get unique clusters
    clusters = adata.obs[cluster_key].unique()
    
    for cluster in clusters:
        cluster_str = str(cluster)
        
        if cluster_str in adata.uns['rank_genes_groups']['names'].dtype.names:
            # Get top genes for this cluster
            genes = adata.uns['rank_genes_groups']['names'][cluster_str][:n_genes_per_cluster]
            logfold = adata.uns['rank_genes_groups']['logfoldchanges'][cluster_str][:n_genes_per_cluster]
            
            # Filter by fold change
            filtered_genes = []
            for gene, lfc in zip(genes, logfold):
                if abs(lfc) >= np.log2(min_fold_change):
                    filtered_genes.append(gene)
            
            if filtered_genes:
                gene_panels[f'cluster_{cluster}'] = filtered_genes
    
    return gene_panels


# =====================================
# STATISTICAL UTILITIES
# =====================================

def calculate_basic_stats(data: np.ndarray) -> Dict[str, float]:
    """
    Calculate basic statistical measures.
    
    Parameters
    ----------
    data : array-like
        Input data
        
    Returns
    -------
    stats : dict
        Dictionary of statistical measures
    """
    data = np.asarray(data)
    
    return {
        'mean': np.mean(data),
        'median': np.median(data),
        'std': np.std(data),
        'min': np.min(data),
        'max': np.max(data),
        'q25': np.percentile(data, 25),
        'q75': np.percentile(data, 75),
        'iqr': np.percentile(data, 75) - np.percentile(data, 25),
        'cv': np.std(data) / np.mean(data) if np.mean(data) > 0 else 0
    }


def detect_outliers(data: np.ndarray, 
                   method: str = 'iqr',
                   threshold: float = 1.5) -> np.ndarray:
    """
    Detect outliers in data.
    
    Parameters
    ----------
    data : array-like
        Input data
    method : str
        Outlier detection method ('iqr', 'zscore', 'mad')
    threshold : float
        Threshold for outlier detection
        
    Returns
    -------
    outliers : array of bool
        Boolean mask indicating outliers
    """
    data = np.asarray(data)
    
    if method == 'iqr':
        q25, q75 = np.percentile(data, [25, 75])
        iqr = q75 - q25
        lower = q25 - threshold * iqr
        upper = q75 + threshold * iqr
        outliers = (data < lower) | (data > upper)
        
    elif method == 'zscore':
        z_scores = np.abs((data - np.mean(data)) / np.std(data))
        outliers = z_scores > threshold
        
    elif method == 'mad':
        median = np.median(data)
        mad = np.median(np.abs(data - median))
        modified_z_scores = 0.6745 * (data - median) / mad
        outliers = np.abs(modified_z_scores) > threshold
        
    else:
        raise ValueError(f"Unknown outlier detection method: {method}")
    
    return outliers


def calculate_correlation_matrix(adata: ad.AnnData,
                               genes: Optional[List[str]] = None,
                               method: str = 'pearson') -> pd.DataFrame:
    """
    Calculate gene-gene correlation matrix.
    
    Parameters
    ----------
    adata : AnnData
        Input data
    genes : list of str, optional
        Genes to include. If None, uses highly variable genes
    method : str
        Correlation method ('pearson', 'spearman')
        
    Returns
    -------
    corr_matrix : DataFrame
        Correlation matrix
    """
    if genes is None:
        if 'highly_variable' in adata.var.columns:
            genes = adata.var_names[adata.var.highly_variable]
        else:
            # Use top 100 most variable genes
            gene_vars = np.var(adata.X.toarray() if hasattr(adata.X, 'toarray') else adata.X, axis=0)
            top_var_idx = np.argsort(gene_vars)[-100:]
            genes = adata.var_names[top_var_idx]
    
    # Filter to available genes
    genes = [g for g in genes if g in adata.var_names]
    
    # Extract expression data
    expr_data = adata[:, genes].X
    if hasattr(expr_data, 'toarray'):
        expr_data = expr_data.toarray()
    
    # Calculate correlation
    expr_df = pd.DataFrame(expr_data, columns=genes)
    corr_matrix = expr_df.corr(method=method)
    
    return corr_matrix


# =====================================
# LOGGING AND REPORTING UTILITIES
# =====================================

def setup_logging(log_level: str = 'INFO',
                 log_file: Optional[Union[str, Path]] = None) -> logging.Logger:
    """
    Setup logging configuration.
    
    Parameters
    ----------
    log_level : str
        Logging level ('DEBUG', 'INFO', 'WARNING', 'ERROR')
    log_file : str or Path, optional
        Path to log file
        
    Returns
    -------
    logger : Logger
        Configured logger
    """
    log_format = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    
    # Clear existing handlers
    for handler in logging.root.handlers[:]:
        logging.root.removeHandler(handler)
    
    # Configure logging
    logging.basicConfig(
        level=getattr(logging, log_level.upper()),
        format=log_format
    )
    
    if log_file:
        log_file = Path(log_file)
        log_file.parent.mkdir(parents=True, exist_ok=True)
        file_handler = logging.FileHandler(log_file)
        file_handler.setFormatter(logging.Formatter(log_format))
        logging.getLogger().addHandler(file_handler)
    
    return logging.getLogger(__name__)


def create_analysis_report(adata: ad.AnnData,
                          analysis_results: Dict[str, Any],
                          output_path: Union[str, Path]) -> None:
    """
    Create comprehensive analysis report.
    
    Parameters
    ----------
    adata : AnnData
        Analyzed data
    analysis_results : dict
        Results from analysis modules
    output_path : str or Path
        Output file path
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    report = []
    report.append("# Spatial Transcriptomics Analysis Report")
    report.append("=" * 50)
    report.append("")
    
    # Dataset overview
    report.append("## Dataset Overview")
    report.append(f"- **Spots**: {adata.n_obs:,}")
    report.append(f"- **Genes**: {adata.n_vars:,}")
    report.append("")
    
    # QC summary
    if 'qc_results' in analysis_results:
        qc = analysis_results['qc_results']
        report.append("## Quality Control")
        report.append(f"- **Mean UMI per spot**: {qc.get('mean_umi', 'N/A')}")
        report.append(f"- **Mean genes per spot**: {qc.get('mean_genes', 'N/A')}")
        report.append("")
    
    # Spatial analysis
    if 'spatial_results' in analysis_results:
        spatial = analysis_results['spatial_results']
        report.append("## Spatial Analysis")
        if 'clustering' in spatial:
            clust = spatial['clustering']
            report.append(f"- **Clustering method**: {clust.get('method', 'N/A')}")
            report.append(f"- **Number of clusters**: {clust.get('n_clusters', 'N/A')}")
        report.append("")
    
    # Gene panels
    if 'gene_panels' in analysis_results:
        panels = analysis_results['gene_panels']
        report.append("## Gene Panel Analysis")
        for panel_name, panel_data in panels.items():
            n_genes = len(panel_data.get('genes_available', []))
            report.append(f"- **{panel_name}**: {n_genes} genes analyzed")
        report.append("")
    
    # Save report
    with open(output_path, 'w') as f:
        f.write('\n'.join(report))
    
    logger.info(f"Analysis report saved to {output_path}")


# =====================================
# VALIDATION UTILITIES
# =====================================

def validate_adata_structure(adata: ad.AnnData) -> Dict[str, Any]:
    """
    Validate AnnData object structure.
    
    Parameters
    ----------
    adata : AnnData
        Data to validate
        
    Returns
    -------
    validation_results : dict
        Validation results and issues found
    """
    results = {'issues': [], 'warnings': []}
    
    # Check basic structure
    if adata.n_obs == 0:
        results['issues'].append("No observations (spots) found")
    if adata.n_vars == 0:
        results['issues'].append("No variables (genes) found")
    
    # Check for spatial coordinates
    if 'spatial' not in adata.obsm:
        results['issues'].append("No spatial coordinates found")
    elif adata.obsm['spatial'].shape[1] != 2:
        results['issues'].append("Spatial coordinates should be 2D")
    
    # Check for required metadata
    required_obs = ['in_tissue']
    for col in required_obs:
        if col not in adata.obs.columns:
            results['warnings'].append(f"Missing recommended column: {col}")
    
    # Check expression data
    if hasattr(adata.X, 'toarray'):
        X_test = adata.X[:100, :100].toarray()
    else:
        X_test = adata.X[:100, :100]
    
    if np.any(X_test < 0):
        results['warnings'].append("Negative values found in expression matrix")
    
    # Check for duplicated gene names
    if adata.var_names.duplicated().any():
        results['warnings'].append("Duplicated gene names found")
    
    return results


def check_memory_usage(adata: ad.AnnData) -> Dict[str, str]:
    """
    Check memory usage of AnnData object.
    
    Parameters
    ----------
    adata : AnnData
        Data to check
        
    Returns
    -------
    memory_info : dict
        Memory usage information
    """
    import sys
    
    def format_bytes(bytes_value):
        for unit in ['B', 'KB', 'MB', 'GB']:
            if bytes_value < 1024:
                return f"{bytes_value:.1f} {unit}"
            bytes_value /= 1024
        return f"{bytes_value:.1f} TB"
    
    memory_info = {}
    
    # Main expression matrix
    if hasattr(adata.X, 'data'):
        # Sparse matrix
        memory_info['X_matrix'] = format_bytes(adata.X.data.nbytes)
        memory_info['X_type'] = 'sparse'
    else:
        # Dense matrix
        memory_info['X_matrix'] = format_bytes(adata.X.nbytes)
        memory_info['X_type'] = 'dense'
    
    # Observations
    obs_size = sum(col.memory_usage(deep=True) for col in adata.obs.values())
    memory_info['obs'] = format_bytes(obs_size)
    
    # Variables
    var_size = sum(col.memory_usage(deep=True) for col in adata.var.values())
    memory_info['var'] = format_bytes(var_size)
    
    # Total object size
    total_size = sys.getsizeof(adata)
    memory_info['total_object'] = format_bytes(total_size)
    
    return memory_info


# =====================================
# EXAMPLE USAGE AND TESTING
# =====================================

def run_pipeline_example():
    """Example of using the utility functions."""
    
    # Setup logging
    logger = setup_logging('INFO', 'pipeline.log')
    logger.info("Starting pipeline example")
    
    # Load default configuration
    config = {
        'qc_thresholds': DEFAULT_QC_THRESHOLDS,
        'preprocessing': DEFAULT_PREPROCESSING_PARAMS,
        'spatial': DEFAULT_SPATIAL_PARAMS
    }
    
    # Get gene panels
    epithelial_genes = get_gene_panel('epithelial', 'breast_cancer')
    logger.info(f"Epithelial panel: {len(epithelial_genes)} genes")
    
    # Example analysis results
    analysis_results = {
        'qc_results': {'mean_umi': 5000, 'mean_genes': 2000},
        'spatial_results': {'clustering': {'method': 'leiden', 'n_clusters': 8}},
        'gene_panels': {'epithelial': {'genes_available': epithelial_genes}}
    }
    
    logger.info("Pipeline example completed")


if __name__ == "__main__":
    run_pipeline_example()
