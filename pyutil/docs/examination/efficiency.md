# Efficiency Examination - pyutil Code Examination

## Summary

This document examines runtime efficiency of the pyutil vertex detection pipeline.
The pipeline has three phases: C++ data marshaling, Python voxelization, and
neural network inference. The inference phase dominates wall time, but the other
phases have avoidable overhead.

---

## 1. C++ Wrapper (`src/SCN_Vertex.cxx`)

### 1.1 Unnecessary copies of input vectors

**Lines 28-31:**
```cpp
auto x = input[0];
auto y = input[1];
auto z = input[2];
auto q = input[3];
```
Each line copies the entire vector. For a point cloud with N points, this
allocates and copies `4 * N * sizeof(float)` bytes unnecessarily. Using
`const auto& x = input[0];` would eliminate this copy entirely.

**Impact:** For a typical point cloud (e.g., 100K points), this wastes ~1.6 MB
of allocation and copy time per call. Small compared to inference but trivially
avoidable.

### 1.2 Repeated `Py_Initialize()` calls

**Line 48:** `Py_Initialize()` is called on every invocation. While the Python
documentation states this is safe to call multiple times (it's a no-op after the
first), there is still a small overhead for the internal check each call.
A static flag or `Py_IsInitialized()` guard would be cleaner and marginally
faster.

### 1.3 Module import on every call

**Line 55:** `PyImport_Import` is called every time `SCN_Vertex` is invoked.
Python caches modules in `sys.modules`, so repeated imports are fast (dict
lookup), but there is still overhead: reference counting, dict lookup, creating
a new reference. Caching the module `PyObject*` across calls (with appropriate
lifetime management) would avoid this.

### 1.4 Data marshaling overhead

**Lines 90-94:** Input data is marshaled from C++ vectors to Python byte objects:
```cpp
size_t input_size = npts * sizeof(FLOAT);
auto pX = PyBytes_FromStringAndSize((const char *) x.data(), input_size);
```
This copies the data. Combined with the vector copy from 1.1, each input array
is copied twice before reaching Python. In Python, `np.frombuffer` (SCN_Vertex.py
line 56) creates a view without copying, but it's viewing the Python bytes object
(which is itself a copy).

**Total copies per input array: 2 (C++ vector copy + bytes object creation).**
Could be reduced to 1 by using `const auto&` references.

---

## 2. Python Voxelization (`python/SCN_Vertex.py`)

### 2.1 `np.unique` with `axis=0` is O(N log N)

**Line 38:**
```python
unique_coords, inverse = np.unique(x, axis=0, return_inverse=True)
```
`np.unique` with `axis=0` sorts the entire 2D array lexicographically. For N
input points in 3D, this is O(N log N) with significant constant factors due to
the multi-column sort. For large point clouds, this is the bottleneck of the
voxelization step.

**Alternative:** A hash-based approach (e.g., using `pandas.DataFrame.drop_duplicates`
or a custom hash on integer coordinates) would be O(N) expected time. However,
`np.unique` also provides sorted output which may be beneficial for memory
locality in downstream operations.

### 2.2 `np.add.at` and `np.maximum.at` are unbuffered and slow

**Lines 44-50:**
```python
np.add.at(ft_out[:, 0], inverse, y[:, 0])
...
np.maximum.at(ft_out[:, i], inverse, y[:, i])
```
The `np.*.at` ufunc methods are explicitly unbuffered, meaning they don't use
SIMD or other optimizations. For large arrays, they are significantly slower
than buffered alternatives.

**Alternative for `add.at`:** `np.bincount(inverse, weights=y[:, 0], minlength=n_voxels)`
is much faster for the summation (it uses optimized C loops internally).

**Alternative for `maximum.at`:** This is harder to replace efficiently in pure
numpy. Options include: `scipy.ndimage.maximum` by label, or a pandas groupby
approach: `pd.Series(y[:, i]).groupby(inverse).max()`.

In practice, for this codebase `n_feat=1` (only charge `q` is the feature), so
the `maximum.at` loop on lines 49-50 never executes (range(1, 1) is empty). Only
the `add.at` on line 44 runs, computing average charge per voxel.

### 2.3 Redundant min computation

**Lines 33, 66:**
```python
# In SCN_Vertex():
coords_offset = coords_np.min(axis=0)       # line 66
# Then inside voxelize():
x = x - x.min(axis=0)                       # line 33
```
The min is computed twice over the same data. The first computes the offset for
later reconstruction; the second (inside `voxelize`) computes and subtracts it.
Passing the pre-computed offset into `voxelize` would save one O(N) pass.

### 2.4 `torch.set_num_threads(1)` limits parallelism

**Line 10:** `torch.set_num_threads(1)` forces single-threaded inference. This
is likely intentional (the caller manages parallelism at a higher level, and
intra-op parallelism for small inference batches can hurt due to thread overhead).
However, for larger point clouds that produce many voxels, multi-threaded
inference could be faster. This is a tuning knob worth benchmarking.

---

## 3. Neural Network Inference

### 3.1 Model loading is cached (good)

**Lines 15-26:** The `_load_model` function caches models by weights path. This
avoids repeated disk I/O and `torch.load` deserialization on subsequent calls.
This is the most impactful optimization already present in the code.

### 3.2 CPU-only inference

**Line 72:** `device = 'cpu'` is hardcoded. If a GPU is available, inference would
be significantly faster (sparse convolutions are well-suited to GPU acceleration).
The data transfer overhead (CPU->GPU->CPU) would likely be amortized for
non-trivial point clouds.

However, GPU usage introduces deployment complexity (CUDA dependencies, memory
management), so CPU-only may be a deliberate choice for portability.

### 3.3 No batching

Each call processes a single point cloud. If multiple point clouds need vertex
detection, they are processed sequentially. Batching multiple inputs into a single
forward pass would amortize fixed overhead (Python call, model setup). However,
sparse convolution batching has its own complexity, and the C++ interface is
designed for single-input calls, so this would require API changes.

### 3.4 Sigmoid activation applied universally

**DeepVtx.py line 27:** `torch.sigmoid(x)` is applied to all output features.
The downstream code (SCN_Vertex.py line 82) computes `pred_np[:, 1] - pred_np[:, 0]`
and takes `argmax`. Since sigmoid is monotonic, the argmax result would be
identical without the sigmoid (i.e., on raw logits). Removing sigmoid would save
one elementwise pass over the output tensor. However, this is negligible compared
to convolution cost.

---

## 4. Overall Pipeline Efficiency Profile

For a typical point cloud of N points:

| Phase | Time Complexity | Dominant Cost |
|-------|----------------|---------------|
| C++ marshaling | O(N) | Memory copies (2x per array) |
| Voxelization | O(N log N) | `np.unique` sort |
| Model inference | O(V * model) | Sparse convolutions (V = voxel count) |
| Post-processing | O(V) | argmax |

**Inference dominates** for any non-trivial input. The voxelization and marshaling
phases are secondary but have room for constant-factor improvements.

### Recommended Priority for Optimization

1. **GPU inference** (if hardware available) - largest potential speedup
2. **Replace `np.add.at` with `np.bincount`** - easy win for voxelization
3. **Remove C++ vector copies** (`const auto&`) - trivial fix
4. **Hash-based voxelization** - O(N) vs O(N log N), matters for large clouds
5. **Batch inference** - amortizes fixed overhead, requires API redesign
