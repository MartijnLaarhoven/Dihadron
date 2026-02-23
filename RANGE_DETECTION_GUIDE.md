# Automatic dEta Range Detection - Quick Reference

## What Changed?

### Before
- Manual `expectedDEtaMin/Max` calculations based on eta bin boundaries
- Plots often showed large empty (zero-valued) regions
- Hard to see correlation patterns due to white space

### After  
- Automatic `actualDEtaMin/max` detected from histogram content
- Only displays bins with actual data
- Correlation structure clearly visible

## Code Location

**File**: [Process_dPhidEta.cxx](Process_dPhidEta.cxx)  
**Function**: `Read_dPhidEta_givenRange_EtaDiff()`  
**Lines**: 820-880

## Algorithm

```cpp
// Step 1: Scan Y-axis (dEta) for data
for (Int_t iy = 1; iy <= nYbins; ++iy) {
    Double_t sumContent = 0;
    for (Int_t ix = 1; ix <= nXbins; ++ix) {
        sumContent += histogram->GetBinContent(ix, iy);
    }
    if (sumContent > 0) {
        if (firstBinWithData == -1) firstBinWithData = iy;
        lastBinWithData = iy;
    }
}

// Step 2: Get actual edges
actualDEtaMin = hist->GetYaxis()->GetBinLowEdge(firstBinWithData);
actualDEtaMax = hist->GetYaxis()->GetBinUpEdge(lastBinWithData);

// Step 3: Apply to display clones only
displayClone->GetYaxis()->SetRangeUser(actualDEtaMin, actualDEtaMax);
```

## Key Points

✅ **Original data preserved**: Full histograms written to files  
✅ **Display only**: Range restriction applied only to canvas clones  
✅ **Dynamic per-bin**: Each eta trigger bin gets its own detected range  
✅ **Robust**: Falls back to calculated range if no data found  

## Expected Console Output

```
Processing trigger eta bin: [0.0, 0.1]
Actual dEta range found: [-1.5, 1.6]
...
Actual dEta range found: [-1.4, 1.5]
...
No data found in Y-axis, using calculated range: [-0.9, 0.9]
```

## Files Affected

- `ProcessOutput/EtaDiff/Mixed_*.root` - Full data with unrestricted axes
- Canvas PDFs in `TemplateFit/EtaDiff/PDFs/` - Display with detected ranges

## Testing

Run the full pipeline to verify:
```bash
root -l Process_dPhidEta.cxx         # Generate dPhidEta histograms
root -l Process_CreateBootstrapSample.cxx  # Bootstrap sampling  
root -l Process_TemplateFit.cxx      # Extract V1/V2 from cleaner data
```

Look for plots that show only the relevant dEta range with no white space.

