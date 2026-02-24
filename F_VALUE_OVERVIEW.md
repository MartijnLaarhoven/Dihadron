# F-Value Overview from Template Fitting

## Overview

The F-parameter is now included in the template fitting output. This document explains what F represents and how to interpret the F-values across different eta regions.

## What is F?

**F is the background/noise amplitude parameter** extracted from the template fitting procedure in the dihadron correlation analysis.

In the template fitting framework:
- The fitting function models the 2D dPhi-dEta correlation with:
  - Flow harmonics: v1, v2, v3, v4 (physics signals)
  - Background components: F (amplitude), G (shape parameter)
  
- **F (par[3])**: Amplitude of the background/baseline shape
  - Represents the overall magnitude of uncorrelated background
  - Large F values indicate significant background contribution relative to the signal
  - Lower F values suggest cleaner correlation measurements

## Where to Find F-Values

When running the template fitting analysis:

```bash
root -l -b -q Process_TemplateFit.cxx
```

The console output now includes F-values for each eta region:

```
print result: LHC25af_pass2_604820_Cent_0_20_TrigEta_[-0.8, -0.7]_AssoEta_[-0.8, -0.75]
v1: [value] +/- [error]
v2: [value] +/- [error]
v3: [value] +/- [error]
v4: [value] +/- [error]
F:  [value] +/- [error]
```

## F-Value Interpretation

### Typical Range
- F values typically range from 2-20 depending on:
  - Trigger pT range
  - Associated pT range
  - Eta region (center vs boundary)
  - Collision system (Pb-Pb, Ne-Ne, pp)
  - Centrality class

### Quality Indicators
- **Stable F values across eta**: Indicates consistent background subtraction
- **Large F error (F_err)**: May indicate:
  - Insufficient statistics in that eta bin
  - Poor fit quality
  - Boundary/edge effects
- **Rapid F variation**: Potential data quality issues or fitting instabilities

### Boundary Eta Regions
F-values may show anomalies at the edges of eta ranges (e.g., eta 0.7-0.8, -0.8 to -0.7):
- Could indicate reduced detector efficiency
- May reflect edge effects in the correlation histograms
- Requires investigation if consistently problematic

## Implementation Details

### Structure Changes (Process_TemplateFit.cxx)

```cpp
struct VnUnit {
    // ... existing v1, v2, v3, v4 and their errors ...
    Double_t F;         // F-parameter value
    Double_t F_err;     // F-parameter error
};
```

### Output Changes

F-values are now:
1. Extracted from the template fitting: `F = par[3]`
2. Stored in the `VnUnit` structure alongside v1, v2, v3, v4
3. Printed to console output immediately after extraction
4. Available for subsequent analysis/plotting

## Usage in Analysis

### Creating F-Value Plots

To compare F-values across eta regions:

```python
# Pseudo-code example for ROOT/Python
import ROOT

# Read template fit results
fResults = ROOT.TFile("path/to/TemplateFit/output.root")

# Extract F-values for each eta bin
F_values = []
F_errors = []
eta_bins = []

for eta_bin in range(n_eta_regions):
    # Access v1-v4 and F from results
    v1 = result_graphs[eta_bin].GetPoint(...)
    F = result_F_values[eta_bin]  # From console output or stored data
    F_errors.append(result_F_errors[eta_bin])
    eta_bins.append(eta_center)
    
# Create comparison plot
g_F = ROOT.TGraphErrors(len(F_values), eta_bins, F_values, ...)
```

### Statistical Checks

For each dataset/centrality combination:
1. Compare F-values across all eta regions
2. Check that F errors remain reasonable (typically < 20% of F value)
3. Identify outliers or anomalies
4. Correlate F variations with flow coefficient (v2) stability

## Example Output Format

```
Analysis: LHC25af_pass2_604820 (Ne-Ne collision)
Centrality: 0-20%

Eta Region Analysis:
Trigger Eta: [-0.8, -0.7]
  - AssocEta [-0.8, -0.75]: F = 8.23 ± 0.45  v2 = 0.00512 ± 0.00089
  - AssocEta [-0.75, -0.7]: F = 7.91 ± 0.43  v2 = 0.00498 ± 0.00091

Trigger Eta: [-0.7, -0.6]
  - AssocEta [-0.8, -0.75]: F = 8.15 ± 0.42  v2 = 0.00520 ± 0.00085
  - AssocEta [-0.75, -0.7]: F = 7.95 ± 0.44  v2 = 0.00505 ± 0.00088

[... continues for all eta regions ...]
```

## Technical Notes

### Reconstruction Details
- F is extracted from the fitting parameter vector at index 3: `F = par[3]`
- F error is from the error vector: `F_err = parerr[3]`
- Parameter order in fitting: [v2, v3, v4, **F**, G, v1]

### Backward Compatibility
- VnUnit structure now includes F with default value 0.0
- Existing code that creates VnUnit with old 8-parameter constructor still works
- New code should pass F and F_err (9th and 10th parameters)

### Debug Information
If issues arise with F-values:
1. Check that template fitting converges (no "invalid fit parameters" warnings)
2. Verify histogram statistics are sufficient (>0.1 integral for eta-diff)
3. Look for "Fit returned all-zero parameters" errors
4. Check for boundary eta issues (empty histograms)

## Related Parameters

For complete background understanding, also review:
- **G**: Background shape parameter (par[4])
- **v1/v2/v3/v4**: Signal components after background subtraction

The combination of F and these other parameters gives a complete picture of the correlation structure.

## Git History

This feature was added with commit:
```
5f7cbed: Add F-value reporting to template fit output
```

Changes focus on:
- Extended VnUnit structure
- Template fit parameter extraction
- Console output reporting

All changes maintain backward compatibility with existing analysis chains.
