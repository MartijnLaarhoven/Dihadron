#!/bin/bash
# Clean up and recreate directory structure for full analysis

echo "Removing old output directories..."
rm -rf ProcessOutput
rm -rf TemplateFit

echo "Creating ProcessOutput structure..."
mkdir -p ProcessOutput/EtaDiff
mkdir -p ProcessOutput/PeripheralSubtraction

echo "Creating TemplateFit structure..."
mkdir -p TemplateFit/PDFs
mkdir -p TemplateFit/EtaDiff/PDFs
mkdir -p TemplateFit/PeripheralSubtracted/EtaDiff/PDFs
mkdir -p TemplateFit/Uncorrected_FourierFit/EtaDiff
mkdir -p TemplateFit/Comparisons/PDFs

echo "Directory structure created successfully!"
echo ""
echo "Ready to run: ./run_full_analysis.sh"

