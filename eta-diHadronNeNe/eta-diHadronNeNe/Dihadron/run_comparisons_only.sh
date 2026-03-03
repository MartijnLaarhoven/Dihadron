#!/bin/bash

# Quick script to regenerate uncorrected v2 extraction and comparison plots
# Run after Process_Bootstrap.cxx has generated bootstrap samples

set -e

echo "=========================================="
echo "Running Uncorrected & Comparison Analysis"
echo "=========================================="

echo ""
echo "Step 1: Extract uncorrected v2 (raw correlations)..."
root -l -b -q Process_Uncorrected_FourierFit.cxx

echo ""
echo "Step 2: Generate comparison plots (all 4 methods)..."
root -l -b -q Plot_AllMethods_Comparison.cxx

echo ""
echo "=========================================="
echo "✓ Uncorrected & Comparison Analysis Complete"
echo "=========================================="
echo ""
echo "Output files:"
echo "  Uncorrected v2: ./TemplateFit/Uncorrected_FourierFit/EtaDiff/Vn_Combined_*_Cent.root"
echo "  Comparison PDFs: ./TemplateFit/Comparisons/PDFs/"
echo ""
