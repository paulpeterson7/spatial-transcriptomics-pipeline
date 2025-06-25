#!/bin/bash
# Setup Docker for Spatial Transcriptomics Pipeline
# This script moves Dockerfile to root and builds the container

echo "🔧 Setting up Docker configuration..."

# Move Dockerfile to root directory
echo "📁 Moving Dockerfile to root..."
mv docker/Dockerfile ./Dockerfile

echo "✅ Dockerfile moved to root"
echo "⚠️  Remember to update COPY paths in Dockerfile!"

# Build Docker container from root
echo "🐳 Building Docker container..."
docker build -t spatial-transcriptomics:latest .

if [ $? -eq 0 ]; then
    echo "✅ Docker container built successfully!"
else
    echo "❌ Docker build failed!"
    exit 1
fi

echo "🚀 Ready to run pipeline!"