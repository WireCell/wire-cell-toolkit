# CUDA Package — Algorithm and Design Explanation

## 1. General Purpose

The `cuda/` package serves as a **build-system validation layer** for
optional CUDA support in the Wire-Cell Toolkit. It is not (yet) a signal
processing component. Its purpose is to:

1. Verify that the Waf build system can compile `.cu` files via `nvcc`
2. Verify that CUDA runtime API calls work on the target system
3. Provide a minimal integration point for future CUDA-accelerated components
4. Define the configuration interface for CUDA-based DFT services

---

## 2. Core Test Algorithm (`simplecudatest.cu`)

### 2.1 Overview

The test performs a GPU round-trip: allocate host data, copy to device, run a
trivial kernel, copy back, and verify. This validates the entire CUDA
toolchain from compilation through execution.

### 2.2 Step-by-step Data Flow

```
Host (CPU)                          Device (GPU)
──────────                          ─────────────
1. malloc foo[128] = {1,1,...,1}
                                    2. cudaMalloc recf[128]
3. cudaMemcpy foo → recf            ──────────────────────►
                                    4. truc kernel: recf[i] <<= 5
                                       (multiply each element by 32)
5. cudaMemcpy recf → foo            ◄──────────────────────
6. print foo[5]  (expect: 32)
7. cudaFree recf
```

### 2.3 The Kernel

```cuda
#define SIZ 128

__global__ void truc(unsigned int* buf)
{
    if (threadIdx.x < SIZ) {
        buf[threadIdx.x] = buf[threadIdx.x] << 5;
    }
    __syncthreads();
}
```

- **`__global__`**: CUDA keyword marking this as a kernel callable from host
- **`threadIdx.x`**: Built-in variable giving each thread its index within
  the block (0-127 here)
- **`<< 5`**: Left bit-shift by 5, equivalent to multiplying by 2^5 = 32
- **Guard `if (threadIdx.x < SIZ)`**: Prevents out-of-bounds access if the
  block is launched with more threads than SIZ (redundant here since the
  launch uses exactly SIZ threads)
- **`__syncthreads()`**: Block-level barrier — all threads in the block wait
  here before proceeding (unnecessary in this case since threads are
  independent)

The kernel is launched as:
```cuda
truc<<<1, SIZ>>>(recf);   // 1 block, 128 threads
```

### 2.4 Error Handling Macros

Two macros provide CUDA error checking:

- **`CUDA_SAFE_CALL(call)`**: Wraps any CUDA API call, checks return code,
  prints file/line on error, exits
- **`CHECKLASTERROR`**: Checks for asynchronous errors from the most recent
  kernel launch via `cudaGetLastError()`

These follow a common CUDA coding pattern. The separation exists because
kernel launches (`<<<>>>` syntax) do not return error codes — errors must be
polled separately.

---

## 3. Build System Integration

### 3.1 `wscript_build`

```python
bld.smplpkg('WireCellCuda', use='CUDA')
```

This single line registers the package with the Waf build system:
- Package name: `WireCellCuda`
- Dependency: `CUDA` (must be detected during `./wcb configure --with-cuda=...`)
- `smplpkg` is a Wire-Cell Waf helper that automatically discovers sources in
  `src/`, headers in `inc/`, and tests in `test/`

### 3.2 Build Command

From `README.org`, the build/test sequence is:

```bash
NVCCFLAGS="-O3" ./wcb configure --with-cuda=/path/to/cuda [...]
./wcb --target=test_simple_cuda -j1 -vv
```

---

## 4. Test Structure

### 4.1 `test_simple_cuda.cxx`

Minimal C++ test that calls `testcuda()` and returns 0. The Waf test
framework treats a zero exit code as pass. Since `testcuda()` calls `exit(1)`
on any CUDA error, a successful run means the full CUDA pipeline works.

### 4.2 `test_simple_nocuda.cxx`

```c
int main() { return 0; }
```

A no-op test that always passes. This exists so that the test suite has a
passing test even when CUDA is not available, preventing build-system errors
from missing test targets.

### 4.3 `test_idft_cufftdft.bats`

```bash
@test "test_idft with torchdft" {
    usepkg aux
    check test_idft cuFftDFT WireCellCuda
}
```

This BATS (Bash Automated Testing System) test runs the `test_idft` program
from the `aux` package with the `cuFftDFT` DFT implementation from
`WireCellCuda`. This tests that the CUDA-based cuFFT DFT service can be used
as a drop-in replacement for the default FFTW-based DFT in inverse-DFT
operations. The test name says "torchdft" but actually tests cuFFT — this
appears to be a copy-paste artifact from a torch-based test.

---

## 5. Configuration Layer (`cfg/layers/high/svcs.jsonnet`)

The jsonnet configuration defines a **service pack** abstraction that lets
users select DFT implementations by platform:

```jsonnet
local dfts = {
    default: self.cpu,
    cpu: { type: "FftwDFT" },
    cuda: { type: "CudaDFT" },
    torchcpu: { type: "TorchDFT", data: { device: "cpu" } },
    torchgpu: { type: "TorchDFT", data: { device: "gpu" } },
};
```

### 5.1 Platform Selection

When a user configures Wire-Cell with `platform="cuda"`, the service pack
returns `{ type: "CudaDFT" }` as the DFT service. This means:

- All signal processing components that need DFT (noise filtering, deconvolution,
  etc.) will use the CUDA-accelerated cuFFT backend
- The switch is transparent — components request a "dft" service and get
  whichever implementation the platform specifies
- No component code changes are needed to run on GPU vs CPU

### 5.2 Service Pack Structure

The function returns an object with:
- `platform`: The selected platform string
- `dft`: The DFT service configuration for that platform
- `random`: A random number generator service (currently CPU-only regardless
  of platform, with a comment noting future GPU implementations)

### 5.3 Relationship to Signal Processing

In the Wire-Cell signal processing pipeline, DFT is a core operation used in:
- **Noise filtering**: Transform waveforms to frequency domain, apply filters,
  transform back
- **Deconvolution**: Apply detector response deconvolution in frequency domain
- **Field response**: Convert field responses between time and frequency domains

The `CudaDFT` service would accelerate all of these operations by running
FFTs on the GPU via NVIDIA's cuFFT library, which is one of the most
performance-critical operations in the signal processing chain.

---

## 6. Architecture Summary

```
cfg/layers/high/svcs.jsonnet
  │
  ├── platform="cpu"  → FftwDFT (CPU, FFTW library)
  ├── platform="cuda" → CudaDFT (GPU, cuFFT library)  ← this package
  ├── platform="torchcpu" → TorchDFT on CPU
  └── platform="torchgpu" → TorchDFT on GPU

cuda/ package (current state):
  └── Build system test only
      ├── simplecudatest.cu  — validates CUDA compilation + runtime
      ├── test_idft_cufftdft.bats — validates cuFFT DFT integration
      └── wscript_build — registers WireCellCuda with build system
```

The actual `CudaDFT` component implementation (referenced in the jsonnet) is
not in this directory — it would be provided by a separate component or a
future expansion of this package. The current package validates that the
infrastructure exists to support such components.
