'''
Vertex finding using list of points
'''
import torch
import sparseconvnet as scn
import numpy as np
from SCN.DeepVtx import DeepVtx

# Set thread count once at module load, not per call
torch.set_num_threads(1)

# Cache loaded models keyed by weights path so we only load from disk once
_model_cache = {}

def _load_model(weights, device='cpu'):
    if weights in _model_cache:
        return _model_cache[weights]
    model = DeepVtx(dimension=3, nIn=1, device=device)
    model.eval()
    trained_dict = torch.load(weights, weights_only=True)
    for param_tensor in trained_dict:
        if trained_dict[param_tensor].shape != model.state_dict()[param_tensor].shape:
            trained_dict[param_tensor] = torch.squeeze(trained_dict[param_tensor], dim=1)
    model.load_state_dict(trained_dict)
    _model_cache[weights] = model
    return model


def voxelize(x, y, resolution=0.5):
    if len(x.shape) != 2:
        raise Exception('x should have 2 dims')

    x = x - x.min(axis=0)
    x = (x / resolution).astype(np.int64)

    # Vectorized: find unique voxels and accumulate features
    # For feature dim 0: average; for dims 1+: max (matches original logic)
    unique_coords, inverse = np.unique(x, axis=0, return_inverse=True)
    n_voxels = len(unique_coords)
    n_feat = y.shape[1]
    ft_out = np.zeros((n_voxels, n_feat), dtype=y.dtype)

    # Average first feature
    np.add.at(ft_out[:, 0], inverse, y[:, 0])
    counts = np.bincount(inverse, minlength=n_voxels).astype(y.dtype)
    ft_out[:, 0] /= counts

    # Max remaining features
    for i in range(1, n_feat):
        np.maximum.at(ft_out[:, i], inverse, y[:, i])

    return unique_coords, ft_out


def SCN_Vertex(weights, x, y, z, q, dtype='float32', resolution=0.5, verbose=False):
    x = np.frombuffer(x, dtype=dtype)
    y = np.frombuffer(y, dtype=dtype)
    z = np.frombuffer(z, dtype=dtype)
    q = np.frombuffer(q, dtype=dtype)
    coords_np = np.stack((x, y, z), axis=1)
    ft_np = np.expand_dims(q, axis=1)
    if verbose:
        print("in: coords: ", coords_np.shape, coords_np.dtype)
        print("in: ft: ", ft_np.shape, ft_np.dtype)

    coords_offset = coords_np.min(axis=0)
    coords_np, ft_np = voxelize(coords_np, ft_np, resolution=resolution)
    if verbose:
        print("vox: coords: ", coords_np.shape, coords_np.dtype)
        print("vox: ft: ", ft_np.shape, ft_np.dtype)

    device = 'cpu'
    model = _load_model(weights, device=device)

    coords = torch.LongTensor(coords_np)
    ft = torch.FloatTensor(ft_np).to(device)

    with torch.no_grad():
        prediction = model([coords, ft])

    pred_np = prediction.cpu().numpy()
    pred_np = pred_np[:, 1] - pred_np[:, 0]

    pred_coord = coords_np[np.argmax(pred_np)]
    if verbose:
        print('raw: pred_coord: ', pred_coord)

    pred_coord = pred_coord.astype(dtype)
    pred_coord *= resolution
    pred_coord += coords_offset + 0.5 * resolution
    if verbose:
        print('final: pred_coord', pred_coord)

    return pred_coord.tobytes()
