# NaiveStriper

Groups active channels in a slice into contiguous stripes per wire plane by connecting adjacent wires (within a configurable gap tolerance) and finding connected subgraph components.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `NaiveStriper` |
| Concrete class | `WireCell::Img::NaiveStriper` |
| Node category | function |
| Primary interface | `ISliceStriper` |
| Input type(s) | `ISlice` |
| Output type(s) | `IStripeSet` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `gap` | maximum number of intervening inactive wires still considered adjacent when forming stripes (default: 1) |
