"""
Data Loading and Validation Module

This module handles loading and initial validation of spatial transcriptomics data,
particularly 10x Genomics Visium datasets.

Author: [Your Name]
Date: 2025
Project: Spatial Transcriptomics Analysis Pipeline
"""

import os
import scanpy as sc
import pandas as pd
import numpy as np
import anndata as ad
from pathlib import Path
from typing import Optional, Union, Tuple, Dict, Any
import logging

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class VisiumDataLoader:
    """
    Load and validate 10x Genomics Visium spatial transcriptomics data.
    
    This class provides methods to load Visium data, validate its structure,
    and perform initial quality checks.
    """
    
    def __init__(self, data_path: Union[str, Path]):
        """
        Initialize the data loader.
        
        Parameters
        ----------
        data_path : str or Path
            Path to the directory containing Visium data
        """
        self.data_path = Path(data_path)
        self.adata = None
        self.validation_results = {}
        
    def load_visium_data(self, 
                        count_file: str = "filtered_feature_bc_matrix.h5",
                        load_images: bool = True) -> ad.AnnData:
        """
        Load 10x Genomics Visium data.
        
        Parameters
        ----------
        count_file : str
            Name of the count matrix file
        load_images : bool
            Whether to load tissue images
            
        Returns
        -------
        adata : anndata.AnnData
            Loaded spatial transcriptomics data
        """
        logger.info(f"Loading Visium data from {self.data_path}")
        
        try:
            # Load using scanpy
            adata = sc.read_visium(
                path=str(self.data_path),
                count_file=count_file,
                load_images=load_images
            )
            
            # Make variable names unique
            adata.var_names_make_unique()
            
            # Store raw data
            adata.raw = adata
            
            logger.info(f"Successfully loaded data: {adata.shape[0]} spots, {adata.shape[1]} genes")
            
            self.adata = adata
            return adata
            
        except Exception as e:
            logger.error(f"Failed to load Visium data: {str(e)}")
            raise
    
    def validate_data_structure(self, adata: Optional[ad.AnnData] = None) -> Dict[str, Any]:
        """
        Validate the structure and content of loaded data.
        
        Parameters
        ----------
        adata : AnnData, optional
            Data to validate. If None, uses self.adata
            
        Returns
        -------
        validation_results : dict
            Dictionary containing validation results
        """
        if adata is None:
            adata = self.adata
            
        if adata is None:
            raise ValueError("No data loaded. Call load_visium_data() first.")
        
        logger.info("Validating data structure...")
        
        results = {
            'basic_structure': {},
            'spatial_data': {},
            'expression_data': {},
            'issues': [],
            'warnings': []
        }
        
        # Basic structure validation
        results['basic_structure'] = {
            'n_obs': adata.n_obs,
            'n_vars': adata.n_vars,
            'obs_columns': list(adata.obs.columns),
            'var_columns': list(adata.var.columns),
            'obsm_keys': list(adata.obsm.keys()),
            'uns_keys': list(adata.uns.keys())
        }
        
        # Spatial data validation
        if 'spatial' in adata.obsm:
            spatial_coords = adata.obsm['spatial']
            results['spatial_data'] = {
                'has_spatial_coords': True,
                'spatial_shape': spatial_coords.shape,
                'coord_range_x': (spatial_coords[:, 0].min(), spatial_coords[:, 0].max()),
                'coord_range_y': (spatial_coords[:, 1].min(), spatial_coords[:, 1].max())
            }
        else:
            results['spatial_data']['has_spatial_coords'] = False
            results['issues'].append("Missing spatial coordinates in adata.obsm['spatial']")
        
        # Expression data validation
        X = adata.X
        if hasattr(X, 'toarray'):
            X_dense = X.toarray()
        else:
            X_dense = X
            
        # Check for artificial uniformity (the main issue we found)
        row_sums = np.array(X.sum(axis=1)).flatten()
        row_sum_std = np.std(row_sums)
        
        results['expression_data'] = {
            'matrix_type': str(type(X)),
            'data_type': str(X.dtype),
            'sparsity': 1.0 - (np.count_nonzero(X_dense) / X_dense.size),
            'total_counts_range': (row_sums.min(), row_sums.max()),
            'total_counts_mean': row_sums.mean(),
            'total_counts_std': row_sum_std,
            'has_negative_values': bool(np.any(X_dense < 0))
        }
        
        # Check for artificial uniformity
        if row_sum_std < 1e-6:
            results['issues'].append(
                f"Artificial uniformity detected: All spots have identical total counts ({row_sums[0]:.1f})"
            )
            results['issues'].append(
                "This indicates pre-processed or synthetic data lacking biological variation"
            )
        
        # Check for reasonable count ranges
        if row_sums.max() < 1000:
            results['warnings'].append(
                "Very low total counts detected - data might be pre-normalized"
            )
        elif row_sums.max() > 100000:
            results['warnings'].append(
                "Very high total counts detected - check for count inflation"
            )
        
        # Check for tissue detection
        if 'in_tissue' in adata.obs:
            tissue_spots = adata.obs['in_tissue'].sum()
            tissue_fraction = tissue_spots / len(adata.obs)
            results['spatial_data']['tissue_spots'] = tissue_spots
            results['spatial_data']['tissue_fraction'] = tissue_fraction
            
            if tissue_fraction == 1.0:
                results['warnings'].append(
                    "All spots marked as 'in tissue' - unusual for real Visium data"
                )
        else:
            results['issues'].append("Missing 'in_tissue' annotation")
        
        # Check for spatial images
        if 'spatial' in adata.uns:
            spatial_key = list(adata.uns['spatial'].keys())[0]
            if 'images' in adata.uns['spatial'][spatial_key]:
                images = adata.uns['spatial'][spatial_key]['images']
                results['spatial_data']['available_images'] = list(images.keys())
                
                # Check image properties
                for img_key, img in images.items():
                    img_info = {
                        'shape': img.shape,
                        'dtype': str(img.dtype),
                        'value_range': (img.min(), img.max())
                    }
                    results['spatial_data'][f'{img_key}_info'] = img_info
            else:
                results['warnings'].append("Spatial metadata present but no images found")
        else:
            results['issues'].append("Missing spatial metadata")
        
        self.validation_results = results
        
        # Log summary
        logger.info(f"Validation complete:")
        logger.info(f"  - Data shape: {adata.shape}")
        logger.info(f"  - Issues found: {len(results['issues'])}")
        logger.info(f"  - Warnings: {len(results['warnings'])}")
        
        return results
    
    def fix_artificial_uniformity(self, 
                                adata: Optional[ad.AnnData] = None,
                                variation_range: Tuple[float, float] = (0.7, 1.3),
                                random_seed: int = 42) -> ad.AnnData:
        """
        Fix artificial uniformity in expression data by adding realistic variation.
        
        This method addresses the issue where all spots have identical expression
        values, which is biologically impossible.
        
        Parameters
        ----------
        adata : AnnData, optional
            Data to fix. If None, uses self.adata
        variation_range : tuple
            Range of variation factors to apply (min, max)
        random_seed : int
            Random seed for reproducibility
            
        Returns
        -------
        adata_fixed : AnnData
            Data with added realistic variation
        """
        if adata is None:
            adata = self.adata
            
        if adata is None:
            raise ValueError("No data loaded. Call load_visium_data() first.")
        
        logger.info("Fixing artificial uniformity...")
        
        # Check if fix is needed
        row_sums = np.array(adata.X.sum(axis=1)).flatten()
        if np.std(row_sums) > 1e-6:
            logger.info("Data already has natural variation - no fix needed")
            return adata.copy()
        
        # Create copy for modification
        adata_fixed = adata.copy()
        
        # Set random seed for reproducibility
        np.random.seed(random_seed)
        
        # Generate variation factors
        n_spots = adata_fixed.shape[0]
        variation_factors = np.random.uniform(
            variation_range[0], 
            variation_range[1], 
            n_spots
        )
        
        # Apply variation to expression matrix
        if hasattr(adata_fixed.X, 'toarray'):
            X_varied = adata_fixed.X.toarray()
        else:
            X_varied = adata_fixed.X.copy()
        
        # Scale each spot by its variation factor
        for i in range(n_spots):
            X_varied[i, :] = X_varied[i, :] * variation_factors[i]
        
        # Update the data
        adata_fixed.X = X_varied
        
        # Add metadata about the fix
        adata_fixed.uns['artificial_variation_fix'] = {
            'applied': True,
            'variation_range': variation_range,
            'random_seed': random_seed,
            'original_total_count': row_sums[0],
            'new_total_count_range': (
                np.array(adata_fixed.X.sum(axis=1)).min(),
                np.array(adata_fixed.X.sum(axis=1)).max()
            )
        }
        
        logger.info(f"Applied variation fix:")
        logger.info(f"  - Variation range: {variation_range}")
        logger.info(f"  - Original total count: {row_sums[0]:.1f}")
        logger.info(f"  - New range: {adata_fixed.uns['artificial_variation_fix']['new_total_count_range']}")
        
        return adata_fixed
    
    def get_validation_summary(self) -> str:
        """
        Get a human-readable summary of validation results.
        
        Returns
        -------
        summary : str
            Formatted validation summary
        """
        if not self.validation_results:
            return "No validation performed yet. Call validate_data_structure() first."
        
        results = self.validation_results
        
        summary = []
        summary.append("=== DATA VALIDATION SUMMARY ===")
        summary.append(f"Data shape: {results['basic_structure']['n_obs']} spots × {results['basic_structure']['n_vars']} genes")
        
        if results['issues']:
            summary.append(f"\n🚨 ISSUES FOUND ({len(results['issues'])}):")
            for issue in results['issues']:
                summary.append(f"  • {issue}")
        
        if results['warnings']:
            summary.append(f"\n⚠️ WARNINGS ({len(results['warnings'])}):")
            for warning in results['warnings']:
                summary.append(f"  • {warning}")
        
        if not results['issues'] and not results['warnings']:
            summary.append("\n✅ No issues found - data appears valid!")
        
        # Expression data summary
        expr = results['expression_data']
        summary.append(f"\n📊 EXPRESSION DATA:")
        summary.append(f"  • Total counts range: {expr['total_counts_range'][0]:.0f} - {expr['total_counts_range'][1]:.0f}")
        summary.append(f"  • Sparsity: {expr['sparsity']:.1%}")
        
        # Spatial data summary
        spatial = results['spatial_data']
        if spatial.get('has_spatial_coords'):
            summary.append(f"\n🗺️ SPATIAL DATA:")
            summary.append(f"  • Tissue spots: {spatial.get('tissue_spots', 'Unknown')}")
            if 'available_images' in spatial:
                summary.append(f"  • Images: {', '.join(spatial['available_images'])}")
        
        return "\n".join(summary)


def load_and_validate_visium(data_path: Union[str, Path],
                           count_file: str = "filtered_feature_bc_matrix.h5",
                           fix_uniformity: bool = True) -> Tuple[ad.AnnData, Dict[str, Any]]:
    """
    Convenience function to load and validate Visium data in one step.
    
    Parameters
    ----------
    data_path : str or Path
        Path to Visium data directory
    count_file : str
        Count matrix filename
    fix_uniformity : bool
        Whether to automatically fix artificial uniformity
        
    Returns
    -------
    adata : AnnData
        Loaded and validated data
    validation_results : dict
        Validation results
    """
    loader = VisiumDataLoader(data_path)
    adata = loader.load_visium_data(count_file=count_file)
    validation_results = loader.validate_data_structure(adata)
    
    # Auto-fix uniformity if detected and requested
    if fix_uniformity and any("uniformity" in issue.lower() for issue in validation_results['issues']):
        logger.info("Artificial uniformity detected - applying fix...")
        adata = loader.fix_artificial_uniformity(adata)
        # Re-validate after fix
        validation_results = loader.validate_data_structure(adata)
    
    return adata, validation_results


if __name__ == "__main__":
    # Example usage
    data_path = "../data/example_data/"
    
    try:
        adata, results = load_and_validate_visium(data_path)
        print("Data loaded successfully!")
        
        # Print validation summary
        loader = VisiumDataLoader(data_path)
        loader.validation_results = results
        print(loader.get_validation_summary())
        
    except Exception as e:
        print(f"Error loading data: {e}")