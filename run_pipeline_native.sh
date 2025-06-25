#!/bin/bash
# Native pipeline runner (runs inside container without Docker)

echo "🧬 Starting Spatial Transcriptomics Pipeline (Native)..."

# Check if we're inside the container
if [ ! -f "/usr/local/bin/nextflow" ]; then
    echo "❌ Nextflow not found! This script should run inside the container"
    echo "💡 Run: docker run -it --rm -v \$(pwd):/workspace -w /workspace spatial-transcriptomics:latest"
    exit 1
fi

echo "✅ Running inside container"

# Check if data file exists
DATA_FILE="data/example_data/Visium_Human_Breast_Cancer_filtered_feature_bc_matrix.h5"
if [ ! -f "$DATA_FILE" ]; then
    echo "❌ Data file not found: $DATA_FILE"
    echo "💡 Current working directory: $(pwd)"
    echo "💡 Available files:"
    find . -name "*.h5" 2>/dev/null || echo "   No .h5 files found"
    exit 1
fi

echo "✅ Data file found: $DATA_FILE"

# Create output directories
mkdir -p results/figures results/processed_data results/reports

# Check if workflow exists
if [ ! -f "workflows/main.nf" ]; then
    echo "❌ Nextflow workflow not found: workflows/main.nf"
    echo "💡 Available workflows:"
    find . -name "*.nf" 2>/dev/null || echo "   No .nf files found"
    exit 1
fi

echo "✅ Workflow found: workflows/main.nf"

echo "🚀 Running Nextflow pipeline..."

# Run Nextflow directly (no Docker wrapper)
nextflow run workflows/main.nf \
    --input_dir ./data/example_data \
    --outdir ./results

# Check if pipeline succeeded
if [ $? -eq 0 ]; then
    echo ""
    echo "🎉 Pipeline completed successfully!"
    echo ""
    echo "📁 Results available in:"
    echo "   • results/figures/ - Generated plots"
    echo "   • results/processed_data/ - Processed datasets"
    echo ""
    echo "🖼️ Generated files:"
    find results/ -name "*.png" -o -name "*.h5ad" -o -name "*.txt" 2>/dev/null | while read file; do
        echo "   📄 $file"
    done
else
    echo ""
    echo "❌ Pipeline failed!"
    echo "💡 Check the error messages above"
fi