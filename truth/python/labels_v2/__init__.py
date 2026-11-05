"""labels_v2 -- the truth labeler v2 of the PR-chain ML training (pr_ml docs 50a / 50b secs 2.1-2.8; built in doc 54).

A separate, versioned python step beside the v1 labeler (PFMerger, TruthPointLabeler, train_targets.py), which it does not
touch.  It reads one campaign event (drifted deposits, truth table, PR point tree, wires, frames) and writes one
<event>.labels_v2.npz with per-point classes, instances, flags, typed key points and a per-instance table.

Modules: io (loading, reco-frame placement), pixels (the 2-D projection, presence, ring, clean share), objects (the truth
tree -> instances), rules (2.2, 2.3, 2.5, 2.6, Q12 A, Q15), points (deposit -> PR point), keypoints (G5), writer, and the
CLI labels_v2.py.  numpy + scipy only.
"""
VERSION = 1
