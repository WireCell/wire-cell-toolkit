# DepoMerger

Merges two time-ordered IDepo streams into one time-ordered output stream by consuming whichever head deposition has the smaller time; EOS is propagated only after both input streams signal end-of-stream.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `DepoMerger` |
| Concrete class | `WireCell::Gen::DepoMerger` |
| Node category | hydra |
| Primary interface | `IDepoMerger` |
| Input type(s) | `(IDepo, IDepo)` |
| Output type(s) | `IDepo` |
| Configurable | no |
