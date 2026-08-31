# PointGraphProcessor

Loads a point-cloud graph from a tensor set, applies a sequence of IPointGraphVisitor plug-ins to it, and writes the processed graph back as a new tensor set.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `PointGraphProcessor` |
| Concrete class | `WireCell::Img::PointGraphProcessor` |
| Node category | function |
| Primary interface | `ITensorSetFilter` |
| Input type(s) | `ITensorSet` |
| Output type(s) | `ITensorSet` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `inpath` | tensor datapath from which the input point-cloud graph is read |
| `outpath` | tensor datapath to which the output point-cloud graph is written; supports "%" substitution with the tensor set ident |
| `visitors` | ordered array of IPointGraphVisitor type/name strings applied in sequence to the point graph |
