# ChargeSolving

Solves for blob charge values by decomposing the cluster into connected blob-measure subgraphs and applying one or more rounds of LASSO regression with configurable blob-weighting strategies and optional Cholesky whitening.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `ChargeSolving` |
| Concrete class | `WireCell::Img::ChargeSolving` |
| Node category | function |
| Primary interface | `IClusterFilter` |
| Input type(s) | `ICluster` |
| Output type(s) | `ICluster` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `meas_value_threshold` | lower bound on measure signal value for a measure to be included in solving (default: 10) |
| `meas_error_threshold` | upper bound on measure signal uncertainty for a measure to be included (default: 1e9) |
| `blob_value_threshold` | lower bound on solved blob charge for a blob to appear in the output (default: -1) |
| `blob_error_threshold` | upper bound on blob uncertainty for output inclusion (default: 0) |
| `lasso_tolerance` | convergence tolerance for the LASSO solver (default: 1e-3) |
| `lasso_minnorm` | minimum norm parameter for the LASSO solver (default: 1e-6) |
| `weighting_strategies` | ordered list of blob-weighting strategy names; choices are "uniform", "simple", "uboone" (default: ["uboone"]) |
| `solve_config` | high-level solver parameter preset; one of "uboone" or "simple" (default: "uboone") |
| `whiten` | if true, apply Cholesky whitening to the measurement covariance matrix before solving (default: true) |
