# GaussianDiffusion Class Analysis

## Core Components

### 1. GausDesc Structure
- Represents a Gaussian distribution with two key parameters:
  - `center`: The mean/center of the Gaussian distribution
  - `sigma`: The standard deviation of the distribution

### 2. GaussianDiffusion Class
- Main purpose: Models the diffusion of charge deposits in a wire chamber detector
- Key members:
  - `m_deposition`: Pointer to the original charge deposition
  - `m_time_desc`: Gaussian description for time dimension
  - `m_pitch_desc`: Gaussian description for spatial/pitch dimension
  - `m_patch`: 2D array storing the diffused charge distribution
  - `m_qweights`: Vector storing weights for charge interpolation
  - `m_toffset_bin`, `m_poffset_bin`: Offset bins for time and pitch dimensions

## Key Algorithms

### 1. Gaussian Sampling (GausDesc::sample)
```cpp
std::vector<double> sample(double start, double step, int nsamples) const {
    if (!sigma) {
        // Handle point source case
        return {1.0};
    }
    // Sample Gaussian at regular intervals
    for (int ind = 0; ind < nsamples; ++ind) {
        const double rel = (start + ind * step - center) / sigma;
        ret[ind] = exp(-0.5 * rel * rel);
    }
}
```

### 2. Bin Integration (GausDesc::binint)
- Uses error function (erf) to compute integrated charge in each bin
- More accurate than simple sampling for charge conservation
```cpp
std::vector<double> binint(double start, double step, int nbins) const {
    if (!sigma) {
        return {1.0};  // Point source case
    }
    // Compute erf differences for bin integration
    const double sqrt2 = sqrt(2.0);
    for (int ind = 0; ind <= nbins; ++ind) {
        double x = (start + step * ind - center) / (sqrt2 * sigma);
        erfs[ind] = 0.5 * std::erf(x);
    }
    // Calculate bin contents
    for (int ibin = 0; ibin < nbins; ++ibin) {
        bins[ibin] = erfs[ibin + 1] - erfs[ibin];
    }
}
```

### 3. Weight Calculation Algorithm (GausDesc::weight)

The weight calculation is a sophisticated algorithm designed to handle linear charge interpolation between impact positions. Here's the detailed breakdown:

#### Purpose
- Provides weights for linear interpolation of charge between adjacent wire positions
- Accounts for the continuous nature of the charge distribution

#### Algorithm Steps
1. For each bin:
   ```cpp
   double x2 = start;
   double x1 = 0;
   double gaus2 = exp(-0.5 * (start - center) / sigma * (start - center) / sigma);
   double gaus1 = 0;
   
   for (int ind = 0; ind < nbins; ind++) {
       x1 = x2;
       x2 = x1 + step;
       double rel = (x2 - center) / sigma;
       gaus1 = gaus2;
       gaus2 = exp(-0.5 * rel * rel);
   ```

2. Weight Calculation Formula:
   ```cpp
   wt[ind] = -1.0 * sigma / (x1 - x2) * (gaus2 - gaus1) / sqrt(2.0 * pi) / pvec[ind] 
           + (center - x2) / (x1 - x2);
   ```

#### Mathematical Explanation
The weight calculation combines two components:
1. Gaussian derivative term: `-1.0 * sigma / (x1 - x2) * (gaus2 - gaus1) / sqrt(2.0 * pi) / pvec[ind]`
   - Represents the rate of change of the Gaussian distribution
   - Normalized by the bin's total charge (pvec[ind])

2. Linear position term: `(center - x2) / (x1 - x2)`
   - Provides linear interpolation based on position relative to bin edges

### 4. Diffusion Sampling (set_sampling method)

The `set_sampling` method combines all these components to create the final diffusion model:

1. Calculate time and pitch ranges based on number of sigmas
2. Sample or integrate both dimensions
3. Create 2D charge distribution patch
4. Apply optional charge fluctuations
5. Normalize to preserve total charge

## Implementation Details

### Charge Conservation
- The implementation carefully preserves total charge through normalization
- Both bin integration and sampling methods are normalized
- Fluctuations (if applied) maintain the total charge through renormalization

### Coordinate Systems
- Uses two coordinate systems:
  1. Absolute coordinates (center, time)
  2. Bin-relative coordinates (offsets)
- Transforms between these systems using offset bins (m_toffset_bin, m_poffset_bin)

### Performance Considerations
- Pre-calculates error functions for efficiency
- Uses vectorized operations where possible
- Caches results in m_patch to avoid recalculation