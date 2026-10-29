"""dv2max_vertexer: inference-only DV2-MAX (r7_cloud) for the in-chain DL-vertex step of Wire-Cell (SBND).

Vendored from sbnd-vertexer dv2max/dv2max and common/sbndvtx (see README.md).  Runtime dependencies: numpy, scipy,
torch (CPU), lightgbm.  Entry point for WCT: pyutil/python/DV2_Vertex.py; the work happens in vertexer.Vertexer.
"""
__version__ = "1.0"
