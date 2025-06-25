#!/bin/bash
# Setup Docker for Spatial Transcriptomics Pipeline
# Builds container from root directory (where Dockerfile should be)

echo "🔧 Setting up Docker configuration..."

# Check if Dockerfile exists in root
if [ ! -f "Dockerfile" ]; then
    echo "❌ Dockerfile not found in root directory!"
    echo "💡 Please ensure Dockerfile is in the project root"
    echo "   Current directory: $(pwd)"
    echo "   Expected: $(pwd)/Dockerfile"
    exit 1
fi

echo "✅ Dockerfile found in root"

# Check if requirements.txt exists in docker/
if [ ! -f "docker/requirements.txt" ]; then
    echo "❌ requirements.txt not found!"
    echo "💡 Expected: docker/requirements.txt"
    exit 1
fi

echo "✅ Requirements file found"

# Check Docker is running
if ! docker info >/dev/null 2>&1; then
    echo "❌ Docker is not running!"
    echo "💡 Please start Docker Desktop and try again"
    exit 1
fi

echo "✅ Docker is running"

# Build Docker container from root directory
echo "🐳 Building Docker container from root..."
echo "   Command: docker build -t spatial-transcriptomics:latest ."

docker build -t spatial-transcriptomics:latest .

if [ $? -eq 0 ]; then
    echo ""
    echo "🎉 Docker container built successfully!"
    echo ""
    echo "🧪 Testing container..."

    # Test that Nextflow is installed
    if docker run --rm spatial-transcriptomics:latest nextflow -version >/dev/null 2>&1; then
        echo "✅ Nextflow is working"
    else
        echo "⚠️ Nextflow test failed"
    fi

    # Test that Python packages are installed
    if docker run --rm spatial-transcriptomics:latest python -c "import scanpy; print('✅ Scanpy available:', scanpy.__version__)" 2>/dev/null; then
        echo "✅ Python packages working"
    else
        echo "⚠️ Python package test failed"
    fi

    echo ""
    echo "🚀 Ready to run pipeline!"
    echo ""
    echo "📋 Next steps:"
    echo "   • Run pipeline: ./run_pipeline.sh"
    echo "   • Or interactive: docker run --rm -v \"\$(pwd):/workspace\" -w /workspace spatial-transcriptomics:latest ./run_pipeline_native.sh"

else
    echo ""
    echo "❌ Docker build failed!"
    echo ""
    echo "🔧 Troubleshooting tips:"
    echo "   • Check that Dockerfile is in root directory"
    echo "   • Ensure docker/requirements.txt exists"
    echo "   • Verify Docker Desktop is running"
    echo "   • Try: docker system prune (clean build)"
    exit 1
fi