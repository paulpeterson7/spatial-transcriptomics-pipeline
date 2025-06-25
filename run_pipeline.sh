#!/bin/bash
# Spatial Transcriptomics Pipeline Runner with Proper Volume Mounting

echo "🧬 Starting Spatial Transcriptomics Pipeline..."

# Check if Docker container exists
if ! docker image inspect spatial-transcriptomics:latest >/dev/null 2>&1; then
    echo "❌ Docker container not found!"
    echo "💡 Run ./setup_docker.sh first to build the container"
    exit 1
fi

echo "✅ Docker container found"

# Check if data file exists
DATA_FILE="data/example_data/Visium_Human_Breast_Cancer_filtered_feature_bc_matrix.h5"
if [ ! -f "$DATA_FILE" ]; then
    echo "❌ Data file not found: $DATA_FILE"
    echo "💡 Please ensure your data file is in the correct location"
    echo "   Expected: $(pwd)/$DATA_FILE"
    exit 1
fi

echo "✅ Data file found: $DATA_FILE"

# Create output directories
mkdir -p results/figures results/processed_data results/reports

echo "🚀 Running Nextflow pipeline in Docker container..."

# Mount the entire project directory and run Nextflow
docker run --rm \
    -v "$(pwd):/workspace" \
    -w /workspace \
    spatial-transcriptomics:latest \
    nextflow run workflows/main.nf \
    --input_dir /workspace/data/example_data \
    --outdir /workspace/results

# Check if pipeline succeeded
if [ $? -eq 0 ]; then
    echo ""
    echo "🎉 Pipeline completed successfully!"
    echo ""
    echo "📁 Results available in:"
    echo "   • results/figures/ - Generated plots"
    echo "   • results/processed_data/ - Processed datasets"
    echo "   • results/reports/ - Analysis summaries"
    echo ""
    echo "🖼️ Key files to check:"
    find results/figures/ -name "*.png" 2>/dev/null | head -5 | while read file; do
        echo "   📈 $file"
    done
else
    echo ""
    echo "❌ Pipeline failed!"
    echo "💡 Check the error messages above for details"
    echo "   • Verify data file exists and is readable"
    echo "   • Check Docker container has necessary permissions"
    exit 1
fi