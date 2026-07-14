# Algorithm Explanation - pyutil Code Examination

## General Overview

The pyutil module performs **3D vertex detection** on point cloud data using a
sparse convolutional neural network. Given a set of 3D points (x, y, z) with
associated charge values (q), it predicts the location of an interaction vertex
- the point in space where a particle interaction occurred.

This is used in the context of **liquid argon time projection chamber (LArTPC)**
detectors, where ionization charge deposits along particle trajectories form 3D
point clouds. Identifying the interaction vertex is a critical reconstruction
step for neutrino physics analysis.

The pipeline has four stages:

```
3D Point Cloud (x, y, z, q)
        |
        v
  [1] Voxelization
        |
        v
  [2] Sparse 3D UNet Inference
        |
        v
  [3] Score Differencing + Argmax
        |
        v
  [4] Coordinate Reconstruction
        |
        v
  Predicted Vertex (x, y, z)
```

---

## Detailed Algorithm

### Stage 0: Data Marshaling (C++ to Python)

**File:** `src/SCN_Vertex.cxx`

The C++ wrapper receives four vectors of floats: x-coordinates, y-coordinates,
z-coordinates, and charge values. It serializes each into raw byte buffers and
passes them to the Python function via the CPython API.

```
C++ vectors -> PyBytes objects -> Python function call -> PyBytes result -> C++ vector
```

The byte representation avoids Python object overhead (no per-element Python float
objects). Each float array is passed as a contiguous block of `N * sizeof(float)`
bytes.

### Stage 1: Voxelization

**File:** `python/SCN_Vertex.py`, function `voxelize()` (lines 29-52)

Voxelization converts a continuous 3D point cloud into a discrete 3D grid,
which is required input for the sparse convolutional network.

#### Step 1a: Coordinate normalization and discretization

```python
x = x - x.min(axis=0)          # Shift to origin
x = (x / resolution).astype(np.int64)   # Discretize
```

Each coordinate is shifted so the minimum is at the origin, then divided by the
resolution (default 0.5 cm) and truncated to integer voxel indices. This maps
continuous space into a regular 3D grid where each cell is `resolution` cm on
a side.

**Example:** A point at (105.3, 200.7, 50.1) with min=(100.0, 200.0, 50.0)
becomes: ((105.3-100.0)/0.5, (200.7-200.0)/0.5, (50.1-50.0)/0.5) = (10, 1, 0)
after integer truncation.

#### Step 1b: Voxel aggregation

Multiple points may map to the same voxel. These are merged:
- **Charge (feature dim 0):** averaged across all points in the voxel
- **Other features (dims 1+):** max-pooled (not used in current configuration)

```python
unique_coords, inverse = np.unique(x, axis=0, return_inverse=True)
np.add.at(ft_out[:, 0], inverse, y[:, 0])
ft_out[:, 0] /= counts
```

`np.unique` with `return_inverse` identifies unique voxels and provides a mapping
from each input point to its voxel index. The inverse mapping is used to scatter
features into voxel buckets.

The output is a sparse representation: only occupied voxels are stored, not the
full 3D grid (which would be enormous at fine resolution).

### Stage 2: Sparse 3D UNet Inference

**File:** `python/SCN/DeepVtx.py`

The neural network is a **sparse 3D UNet** built with the `sparseconvnet` library.
It operates directly on sparse data (only occupied voxels), avoiding the
computational cost of processing empty space.

#### Architecture

```
Input: (coords [N x 3], features [N x 1])
  |
  v
InputLayer (mode=3: sparse tensor, dimension=3)
  |
  v
SubmanifoldConvolution (1 -> 16 features, 3x3x3 kernel)
  |
  v
UNet (feature planes: [16, 32, 64, 128, 256])
  |  - Encoder: progressively downsample by factor 2, expand features
  |  - Bottleneck: 256 features at coarsest resolution
  |  - Decoder: progressively upsample, concatenate skip connections
  |  - Each level has 2 convolution blocks (reps=2)
  |
  v
BatchNormReLU (16 features)
  |
  v
OutputLayer (sparse -> dense per-voxel features)
  |
  v
Linear (16 -> 2 classes)
  |
  v
Sigmoid (output in [0, 1])
  |
  v
Output: [N x 2] per-voxel class probabilities
```

**Key architectural choices:**

- **Sparse convolutions:** Only compute at occupied voxels and their neighbors.
  Submanifold convolutions maintain the sparsity pattern (output is non-zero only
  where input is non-zero). Regular sparse convolutions (in the UNet
  downsampling) can create new active sites.

- **Mode 3 InputLayer:** Accumulates features from duplicate coordinates. Since
  voxelization already handles deduplication, this mode handles edge cases where
  the model receives pre-aggregated data.

- **UNet skip connections:** The encoder-decoder structure with skip connections
  allows the network to combine local fine-grained features with global context.
  This is important because vertex detection requires both local charge patterns
  and global topology.

- **2-class output:** The network outputs two scores per voxel: background (class 0)
  and vertex (class 1). After sigmoid activation, these are pseudo-probabilities.

- **Feature pyramid [16, 32, 64, 128, 256]:** Feature count doubles at each
  downsampling level, standard UNet practice. The spatial resolution halves while
  the feature richness doubles, maintaining roughly constant computational cost
  per level.

### Stage 3: Score Differencing and Vertex Selection

**File:** `python/SCN_Vertex.py`, lines 81-84

```python
pred_np = prediction.cpu().numpy()
pred_np = pred_np[:, 1] - pred_np[:, 0]     # vertex_score - background_score
pred_coord = coords_np[np.argmax(pred_np)]   # voxel with highest score
```

The network outputs two scores per voxel: P(background) and P(vertex). The
difference `P(vertex) - P(background)` serves as a confidence metric. The voxel
with the maximum confidence is selected as the predicted vertex location.

**Why difference instead of just class 1?** The difference penalizes voxels where
both scores are high (ambiguous) and rewards voxels where the network is
confident it's a vertex. Since sigmoid outputs are independent (not softmax),
a voxel could have high scores for both classes. The difference resolves this.

**Range:** After sigmoid, each score is in [0, 1], so the difference is in
[-1, 1]. A score of +1 means maximum vertex confidence; -1 means maximum
background confidence.

### Stage 4: Coordinate Reconstruction

**File:** `python/SCN_Vertex.py`, lines 88-90

```python
pred_coord = pred_coord.astype(dtype)
pred_coord *= resolution
pred_coord += coords_offset + 0.5 * resolution
```

The predicted voxel index is converted back to physical coordinates:

1. **Multiply by resolution:** Converts from voxel indices back to physical units
2. **Add offset + half resolution:** Reconstructs the original coordinate system.
   The `coords_offset` restores the original origin. The `0.5 * resolution` term
   places the prediction at the **center** of the voxel rather than at its corner.

**Example:** If voxel index is (10, 1, 0), resolution=0.5, offset=(100.0, 200.0, 50.0):
- 10 * 0.5 + 100.0 + 0.25 = 105.25
- 1 * 0.5 + 200.0 + 0.25 = 200.75
- 0 * 0.5 + 50.0 + 0.25 = 50.25

The predicted vertex is (105.25, 200.75, 50.25).

---

## Key Design Decisions

### Why sparse convolutions?

LArTPC point clouds are extremely sparse - ionization occurs only along particle
tracks, which occupy a tiny fraction of the detector volume. A dense 3D CNN
would waste >99% of computation on empty voxels. Sparse convolutions process
only occupied voxels, making the computation proportional to the signal, not
the detector volume.

### Why voxelization before inference?

Sparse convolutional networks require integer grid coordinates. The raw point
cloud has continuous coordinates with non-uniform density. Voxelization:
1. Maps to a regular grid (required by convolution)
2. Reduces redundancy (multiple hits in one voxel are aggregated)
3. Controls the spatial resolution (resolution parameter)

### Why single-point prediction?

The algorithm predicts exactly one vertex location. This is appropriate for
neutrino interaction reconstruction where there is typically one primary vertex
per event. The argmax selection is simple and robust when the network has high
confidence. For events with multiple vertices or ambiguous topologies, more
sophisticated selection (e.g., local maxima finding, clustering) would be needed.

### Why CPU inference with single thread?

The typical usage context is within a larger reconstruction framework that may
process many events in parallel (one per thread/process). Single-threaded CPU
inference avoids contention between the ML inference and the rest of the
reconstruction pipeline. GPU inference would provide a speedup but adds
deployment complexity and potential resource contention.

---

## Input/Output Specification

**Input:**
- `x, y, z`: 3D coordinates of charge deposits (float32 arrays, length N)
- `q`: Charge values at each point (float32 array, length N)
- `weights`: Path to trained model weights file (.pth)
- `resolution`: Voxel size in the same units as coordinates (default 0.5)

**Output:**
- 3 float32 values: predicted vertex (x, y, z) in the original coordinate system

**Constraints:**
- After voxelization, coordinates must fit within [0, 4096) in each dimension
  (the `spatialSize` of the InputLayer)
- At resolution=0.5, this means input range up to ~2048 units per dimension
- The model expects exactly 1 input feature (charge)
