#!/bin/bash
# Full analysis chain for dihadron correlation
# Run this to process all methods with existing datasets
# Execution order: Data processing -> Bootstrap -> Fits -> Diagnostics -> Plots

set -e  # Exit on error

echo "============================================"
echo "Starting Full Dihadron Analysis Pipeline"
echo "============================================"

# Step 1: Process dPhi vs dEta correlations (read from raw data, output SE/ME histograms)
echo ""
echo "Step 1/8: Processing dPhi vs dEta correlations..."
root -l -b -q Process_dPhidEta.cxx

# Step 2: Create Bootstrap samples for standard Template Fit method
echo ""
echo "Step 2/8: Creating Bootstrap samples for Template Fit..."
root -l -b -q Process_CreateBootstrapSample.cxx

# Step 3: Run Template Fit method
echo ""
echo "Step 3/8: Running Template Fit method..."
root -l -b -q Process_TemplateFit.cxx

# Step 4: Compute ZYAM baseline for comparison
echo ""
echo "Step 4/8: Computing ZYAM baseline corrections..."
root -l -b -q Process_PeripheralSubtraction.cxx

# Step 5: Create Bootstrap samples for ZYAM peripheral subtraction
echo ""
echo "Step 5/8: Creating Bootstrap samples for peripheral subtraction..."
root -l -b -q Process_CreateBootstrapSample_PeripheralSubtracted.cxx

# Step 6: Run Template Fit on ZYAM-subtracted data
echo ""
echo "Step 6/8: Running Template Fit on peripheral-subtracted data..."
root -l -b -q Process_TemplateFit_PeripheralSubtracted.cxx

# Step 7: Compute uncorrected Fourier fit (baseline - no corrections)
echo ""
echo "Step 7/8: Computing uncorrected Fourier fit (baseline for correction comparison)..."
root -l -b -q Process_Uncorrected_FourierFit.cxx 2>&1 | tail -20

# Step 8: Generate comparison plots (remaining methods)
echo ""
echo "Step 8/8: Generating comparison plots..."
root -l -b -q Plot_AllMethods_Comparison.cxx

echo ""
echo "============================================"
echo "Analysis pipeline completed successfully!"
echo "============================================"
echo "Output locations:"
echo "  - Uncorrected baseline: TemplateFit/Uncorrected_FourierFit/EtaDiff/"
echo "  - Template Fit results: TemplateFit/EtaDiff/"
echo "  - ZYAM results: TemplateFit/PeripheralSubtracted/EtaDiff/"
echo "  - Comparison plots: TemplateFit/Comparisons/PDFs/"
echo ""
