# Spatial Transcriptomics Analysis Pipeline

A comprehensive, production-ready pipeline for spatial transcriptomics data analysis, particularly optimized for 10x Genomics Visium datasets. This pipeline implements best practices for quality control, preprocessing, spatial analysis, and visualization.

Inspired by the regulatory logic and spatial enhancer analysis framework presented in:
- Regner et al. (2025). *Defining the regulatory logic of breast cancer using single-cell epigenetic and transcriptome profiling*. *Genome Research*. [PubMed](https://pubmed.ncbi.nlm.nih.gov/39914387/)

## Overview

This pipeline was developed as part of a bioinformatics portfolio to demonstrate:
- **Spatial transcriptomics expertise** with real-world Visium data
- **Workflow engineering** using Nextflow for scalability and reproducibility
- **Containerization** with Docker for cross-platform compatibility
- **Modular design** with clean Python modules and comprehensive testing
- **Data quality assessment** including detection and correction of common issues

## Quick Start

### Method 1: Using Docker (Recommended)

```bash
# Clone the repository
git clone https://github.com/[your-username]/spatial-transcriptomics-pipeline.git
cd spatial-transcriptomics-pipeline

# Build and setup Docker container
./setup_docker.sh

# Run the pipeline
./run_pipeline.sh
```

**Alternative manual Docker setup:**
```bash
# Build Docker container manually
docker build -t spatial-transcriptomics:latest .

# Run pipeline manually
docker run --rm -v $(pwd):/workspace -w /workspace \
    spatial-transcriptomics:latest \
    nextflow run workflows/main.nf \
    --input_dir /workspace/data/example_data \
    --outdir /workspace/results
```

### Method 2: Native Execution (Inside Container)

If you want to run multiple commands or debug inside the container:

```bash
# Build container first
./setup_docker.sh

# Enter container interactively
docker run -it --rm -v $(pwd):/workspace -w /workspace spatial-transcriptomics:latest

# Then run native pipeline inside container
./run_pipeline_native.sh
```

### Method 3: Using Conda/Mamba (Local Installation)

```bash
# Create conda environment
conda env create -f environment.yml
conda activate spatial-transcriptomics

# Install Nextflow
curl -s https://get.nextflow.io | bash
sudo mv nextflow /usr/local/bin/

# Run pipeline natively
./run_pipeline_native.sh
```

## Script Overview

- **`setup_docker.sh`**: Builds the Docker container and handles Dockerfile setup
- **`run_pipeline.sh`**: Runs the complete pipeline using Docker (recommended for most users)
- **`run_pipeline_native.sh`**: Runs pipeline directly with Nextflow (for development/debugging)

## Pipeline Features

### Data Loading & Validation
- Automatic data structure validation
- Detection of artificial uniformity (common in demo datasets)
- Spatial coordinate and image verification
- Comprehensive data quality reporting

### Quality Control
- Standard QC metrics (UMI counts, genes detected, mitochondrial %)
- Outlier detection using robust statistical methods
- Interactive QC visualizations in spatial coordinates
- Gene and spot filtering with customizable thresholds

### Preprocessing
- Normalization (total count, log transformation)
- Highly variable gene selection (Seurat v3, Cell Ranger methods)
- Data scaling and centering
- Principal component analysis with variance explained plots

### Spatial Analysis
- Spatial neighborhood graphs
- Clustering (Leiden algorithm)
- Marker gene identification
- Spatial pattern analysis
- Multi-gene panel analysis

### Visualization
- High-quality spatial plots with H&E tissue background
- Multi-gene expression panels
- Quality control dashboards
- Clustering and marker gene visualizations
- Publication-ready figures (300 DPI, vector formats)

## Architecture

```
├── modules/                    # Core Python modules
│   ├── data_loader.py         # Data loading and validation
│   ├── quality_control.py     # QC metrics and filtering
│   ├── preprocessing.py       # Normalization and feature selection
│   ├── spatial_analysis.py    # Spatial pattern analysis
│   ├── visualization.py       # Plot generation
│   └── utils.py              # Helper functions
├── workflows/                 # Nextflow workflow definitions
│   ├── main.nf               # Main pipeline workflow
│   ├── nextflow.config       # Pipeline configuration
│   └── modules/              # Nextflow module definitions
├── notebooks/                # Interactive analysis notebooks
│   ├── 01_exploratory_analysis.ipynb
│   ├── 02_quality_control.ipynb
│   └── 03_spatial_patterns.ipynb
├── tests/                    # Unit and integration tests
├── docker/                   # Container definitions
│   ├── Dockerfile
│   └── requirements.txt
├── data/                     # Example and reference data
└── docs/                     # Documentation
```

## Requirements

### System Requirements
- **Memory**: 8+ GB RAM (16+ GB recommended)
- **Storage**: 10+ GB free space
- **OS**: Linux, macOS, or Windows (with WSL2)

### Software Dependencies
- **Python**: 3.11+
- **Nextflow**: 22.10.0+
- **Docker**: 20.10+ (optional but recommended)
- **Java**: 11+ (for Nextflow)

### Key Python Packages
- `scanpy >= 1.10.0` - Core spatial transcriptomics analysis
- `squidpy >= 1.3.0` - Spatial analysis extensions
- `pandas >= 2.0.0` - Data manipulation
- `matplotlib >= 3.7.0` - Visualization
- `seaborn >= 0.12.0` - Statistical visualization

## Usage Examples

### Basic Analysis
```bash
nextflow run workflows/main.nf \
    --input_dir data/visium_sample \
    --outdir results/basic_analysis
```

### Custom Quality Control
```bash
nextflow run workflows/main.nf \
    --input_dir data/visium_sample \
    --min_counts 2000 \
    --max_counts 30000 \
    --max_pct_mt 15.0 \
    --outdir results/strict_qc
```

### Advanced Spatial Analysis
```bash
nextflow run workflows/main.nf \
    --input_dir data/visium_sample \
    --gene_panels "immune,stromal,epithelial,cancer" \
    --resolution 0.3 \
    --n_neighbors 15 \
    --outdir results/spatial_analysis
```

### Skip Steps for Debugging
```bash
nextflow run workflows/main.nf \
    --input_dir data/visium_sample \
    --skip_preprocessing \
    --skip_spatial_analysis \
    --outdir results/qc_only
```

## Output Structure

```
results/
├── figures/                   # All generated plots
│   ├── qc_distributions.png
│   ├── spatial_qc_metrics.png
│   ├── highly_variable_genes.png
│   ├── pca_variance.png
│   ├── spatial_clusters.png
│   └── gene_expression_panels.png
├── tables/                    # Summary statistics
│   ├── qc_metrics.csv
│   ├── hvg_list.csv
│   └── cluster_markers.csv
├── reports/                   # Analysis reports
│   ├── pipeline_report.html
│   └── quality_control_summary.txt
├── processed_data/            # Intermediate data files
│   ├── adata_qc.h5ad
│   ├── adata_processed.h5ad
│   └── adata_final.h5ad
└── pipeline_info/             # Execution metadata
    ├── execution_report.html
    ├── execution_timeline.html
    └── pipeline_dag.svg
```

## Configuration Options

### Quality Control Parameters
```groovy
params {
    min_counts = 1000        // Minimum UMI counts per spot
    max_counts = 50000       // Maximum UMI counts per spot
    min_genes = 500          // Minimum genes per spot
    max_pct_mt = 20.0        // Maximum mitochondrial gene %
}
```

### Preprocessing Parameters
```groovy
params {
    target_sum = 10000       // Normalization target sum
    n_top_genes = 2000       // Number of HVGs to select
    hvg_method = "seurat_v3" // HVG selection method
    n_pcs = 50              // Number of principal components
}
```

### Spatial Analysis Parameters
```groovy
params {
    resolution = 0.5         // Clustering resolution
    n_neighbors = 10         // Number of neighbors for graph
    gene_panels = "breast_cancer,immune,stromal"
}
```

## Testing

### Run Unit Tests
```bash
pytest tests/ -v --cov=modules
```

### Test with Example Data
```bash
nextflow run workflows/main.nf -profile test
```

### Continuous Integration Test
```bash
nextflow run workflows/main.nf -profile ci
```

## Documentation

- **[Methods Documentation](docs/METHODS.md)** - Detailed methodology
- **[API Documentation](docs/API.md)** - Code documentation
- **[Troubleshooting Guide](docs/TROUBLESHOOTING.md)** - Common issues
- **[Contributing Guidelines](docs/CONTRIBUTING.md)** - Development guide

## Key Innovations

### Data Quality Assessment
This pipeline includes sophisticated data validation that can:
- Detect artificial uniformity in demo/synthetic datasets
- Automatically fix unrealistic data distributions
- Validate spatial coordinates and image alignment
- Generate comprehensive data quality reports

### Flexible Workflow Design
- **Modular architecture** - Each step can be run independently
- **Skip functionality** - Bypass steps for debugging or custom workflows
- **Multiple execution modes** - Local, cluster, cloud deployment
- **Comprehensive logging** - Detailed execution tracking and error reporting

### Production-Ready Features
- **Container support** - Docker and Singularity compatibility
- **Resource management** - Automatic scaling and retry logic
- **Error handling** - Robust failure recovery
- **Reproducibility** - Version-controlled dependencies and parameters

## Scientific Applications

This pipeline is designed for various spatial transcriptomics research applications:

### Cancer Research
- Tumor microenvironment analysis
- Cancer-stroma interactions
- Immune infiltration patterns
- Biomarker spatial distribution

### Developmental Biology
- Tissue morphogenesis studies
- Cell fate mapping
- Spatial gene regulatory networks
- Organ development analysis

### Neuroscience
- Brain region characterization
- Neuronal cell type distribution
- Disease progression mapping
- Drug response analysis

## Educational Value

This pipeline serves as a comprehensive example of:
- Modern bioinformatics workflow development
- Spatial data analysis best practices
- Reproducible research methodologies
- Software engineering in computational biology

## Troubleshooting

### Common Issues

#### Memory Errors
```bash
# Increase memory allocation
nextflow run workflows/main.nf --max_memory 32.GB
```

#### Data Loading Failures
```bash
# Check data structure
python -c "
from modules.data_loader import load_and_validate_visium
adata, results = load_and_validate_visium('data/your_data/')
print(results)
"
```

#### Uniform Expression Detection
```bash
# Enable automatic fix
nextflow run workflows/main.nf --fix_uniformity true
```

#### Missing Dependencies
```bash
# Rebuild container
docker build --no-cache -t spatial-transcriptomics:latest docker/
```

### Getting Help

1. **Check the logs**: `results/pipeline_info/execution_report.html`
2. **Validate input data**: Use the data validation module
3. **Run with debug profile**: `nextflow run -profile debug`
4. **Check resource usage**: Monitor memory and CPU utilization

## Contributing

We welcome contributions! Please see [CONTRIBUTING.md](docs/CONTRIBUTING.md) for guidelines.

### Development Setup
```bash
# Clone repository
git clone https://github.com/[your-username]/spatial-transcriptomics-pipeline.git
cd spatial-transcriptomics-pipeline

# Install development dependencies
pip install -e ".[dev]"

# Run tests
pytest tests/ -v

# Format code
black modules/ tests/
flake8 modules/ tests/
```

### Adding New Features
1. **Create feature branch**: `git checkout -b feature/new-analysis`
2. **Add tests**: Include unit tests for new functionality
3. **Update documentation**: Add to relevant docs
4. **Submit PR**: Include description and test results

## Performance Benchmarks

| Dataset Size | Processing Time | Memory Usage | Output Size |
|-------------|----------------|-------------|-------------|
| Small (1K spots) | 5-10 minutes | 2-4 GB | 100-200 MB |
| Medium (5K spots) | 15-30 minutes | 4-8 GB | 500 MB - 1 GB |
| Large (10K+ spots) | 30-60 minutes | 8-16 GB | 1-2 GB |

## Roadmap

### Version 1.1 (Planned)
- Multi-sample analysis support
- Integration with scRNA-seq reference datasets
- Advanced spatial statistics (Moran's I, spatial autocorrelation)
- Interactive visualizations with Plotly/Bokeh

### Version 1.2 (Planned)
- GPU acceleration for large datasets
- Cloud deployment templates (AWS, GCP, Azure)
- Real-time analysis capabilities
- Machine learning integration for pattern recognition

### Version 2.0 (Future)
- Multi-modal analysis (spatial + ATAC-seq)
- Temporal analysis support
- 3D spatial reconstruction
- AI-powered tissue annotation

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

- Scanpy development team for the core spatial analysis framework
- 10x Genomics for Visium technology and example datasets
- Nextflow community for workflow management tools
- Open source community for foundational tools and libraries


## Related Projects

- **[Scanpy](https://scanpy.readthedocs.io/)** - Single-cell analysis toolkit
- **[Squidpy](https://squidpy.readthedocs.io/)** - Spatial molecular data analysis
- **[Nextflow](https://www.nextflow.io/)** - Workflow management system
- **[10x Genomics](https://www.10xgenomics.com/products/spatial-gene-expression)** - Visium spatial transcriptomics
