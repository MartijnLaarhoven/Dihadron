# Automatic dEta Range Detection Implementation

## Overview
The dPhidEta plots now automatically detect and display only the relevant dEta range where data exists, eliminating empty zero-padded regions that made correlation patterns unreadable.

## Implementation Details

### Feature Location
[Process_dPhidEta.cxx](Process_dPhidEta.cxx#L820-L845) - Lines 820-845 in `Read_dPhidEta_givenRange_EtaDiff()` function

### How It Works

1. **Data Scanning** (Lines 820-831):
   - Iterates through all dEta bins in the 2D histogram
   - For each Y-axis (dEta) bin, sums all content across the X-axis (dPhi)
   - Tracks the first and last bins that contain non-zero content

2. **Range Calculation** (Lines 834-845):
   - If data bins found:
     - `actualDEtaMin` = lower edge of first bin with data
     - `actualDEtaMax` = upper edge of last bin with data
     - Logs: "Actual dEta range found: [min, max]"
   - If no data found (fallback):
     - Uses pre-calculated expected range
     - Logs: "No data found in Y-axis, using calculated range"

3. **Display Application** (Lines 872-880):
   - Creates clones of histograms for drawing
   - Sets Y-axis range to `actualDEtaMin` and `actualDEtaMax` ONLY on display clones
   - **Critical**: Original histograms keep full data range for file writing and downstream processing
   - Both Same/Mixed event and Single/All histograms use the same detected range

### Code Example
```cpp
// Lines 872-880: Apply automatic range to display clones
TH2D* hPhiEtaSMsum_draw = (TH2D*)hPhiEtaSMsum->Clone("hPhiEtaSMsum_draw");
hPhiEtaSMsum_draw->GetYaxis()->SetRangeUser(actualDEtaMin, actualDEtaMax);
hPhiEtaSMsum_draw->Draw("surf1");
```

## Key Design Principles

### ✅ Preserve Original Data
- Original histograms written to files with **full dEta range**
- Bootstrap sampling and template fitting use complete data
- Automatic range detection is display-only

### ✅ Symmetric Eta Handling
- Negative trigger eta [-0.8, -0.7] → Associated eta [0.7, 0.8]
- Positive trigger eta [0.0, 0.8] → Associated eta [-0.8, -0.7]
- Automatic range detects actual correlation patterns for each bin

### ✅ Dynamic Adaptation
- Each trigger eta bin gets its own detected range
- No manual range calculations needed
- Pattern-dependent visualization

## Expected Benefits

1. **Improved Plot Readability**
   - No empty/zero regions dominating the visualization
   - Correlation structures clearly visible
   - Better for presentations and analysis

2. **Automatic Tuning**
   - No need to pre-calculate expected dEta ranges
   - Handles varying data density automatically
   - Robust to different dataset characteristics

3. **Backward Compatible**
   - Full data preserved in output files
   - Downstream analysis unchanged
   - Only display is affected

## Technical Stack

- ROOT 6.x: TH2D histogram manipulation
- C++17: Structured bindings (range declaration)
- THnSparseF: Sparse histogram projections

## Files Modified

- [Process_dPhidEta.cxx](Process_dPhidEta.cxx) - Lines 820-880

## Testing Notes

The automatic range detection requires the full ROOT analysis environment with O2Physics support. Run with:
```bash
cd /home/martijn-laarhoven/Work/dihadronanalysis-master/Dihadron
root -l Process_dPhidEta.cxx
```

Check console output for:
```
Actual dEta range found: [min, max]
```

## Git History

- Commit: Ensure automatic dEta range detection is used for display clones
- Previous implementation: Fixed axis range setting order (project before restricting axes)
- Foundation: Dynamic associated eta parameters (etaAssoMin, etaAssoMax)

