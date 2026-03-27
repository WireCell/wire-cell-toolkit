'''
Vertex finding using list of points
'''
import torch
import sparseconvnet as scn
import numpy as np
from SCN.DeepVtx import DeepVtx


def Test(weights, x, y, z, q, dtype):
    '''
    IO test
    '''
    print("python: Test")
    print("weights: ", weights)
    x = np.frombuffer(x, dtype=dtype)
    y = np.frombuffer(y, dtype=dtype)
    z = np.frombuffer(z, dtype=dtype)
    q = np.frombuffer(q, dtype=dtype)
    coords = np.stack((x, y, z), axis=1)
    ft = np.expand_dims(q, axis=1)
    print("coords: ", coords)
    print(" ft: ", ft)

    vx = np.mean(x)
    vy = np.mean(y)
    vz = np.mean(z)

    vtx = np.array([vx, vy, vz])
    return vtx.tobytes()

def voxelize(x, y, resolution=0.5) :
    if len(x.shape) != 2:
        raise Exception('x should have 2 dims')

    # all voxel indices needs to be non-negative
    x = x - x.min(axis=0)

    # in unit of resolution
    x = x/resolution

    # digitize
    x = x.astype(int)

    # filling histogram
    d = dict()
    w = dict()
    for idx in range(x.shape[0]) :
        key = tuple(x[idx,])
        if key in d :
            d[key][0] = d[key][0] + y[idx,][0]
            w[key][0] = w[key][0] + 1
            d[key][1:] = np.maximum(d[key][1:],y[idx,][1:])
            w[key][1:] = np.ones_like(y[idx,][1:])
        else :
            d[key] = np.copy(y[idx,])
            w[key] = np.ones_like(y[idx,])

    keys = []
    vals = []
    for key in d :
        keys.append(list(key))
        vals.append(list(d[key]/w[key]))

    coords = np.array(keys)
    ft = np.array(vals)

    return coords, ft


def SCN_Vertex(weights, x, y, z, q, dtype='float32', resolution=0.5, verbose=False):
    x = np.frombuffer(x, dtype=dtype)
    y = np.frombuffer(y, dtype=dtype)
    z = np.frombuffer(z, dtype=dtype)
    q = np.frombuffer(q, dtype=dtype)
    coords_np = np.stack((x, y, z), axis=1)
    ft_np = np.expand_dims(q, axis=1)
    if verbose :
        print("in: coords: ", coords_np.shape, coords_np.dtype)
        print("in: ft: ", ft_np.shape, ft_np.dtype)

    coords_offset = coords_np.min(axis=0)
    coords_np, ft_np = voxelize(coords_np, ft_np, resolution=resolution)
    if verbose :
        print("vox: coords: ", coords_np.shape, coords_np.dtype)
        print("vox: ft: ", ft_np.shape, ft_np.dtype)

    torch.set_num_threads(1)
    device = 'cpu'

    coords = torch.LongTensor(coords_np)
    ft = torch.FloatTensor(ft_np).to(device)

    nIn = 1
    model = DeepVtx(dimension=3, nIn=nIn, device=device)
    model.train()

    # torch 1.0.0 seems to have 3 dims for some tensors while 1.3.1 have 4 for them
    # in that case dim=1 is an unsqueezed dim with size only 1
    trained_dict = torch.load(weights, weights_only=True)
    for param_tensor in trained_dict:
        if trained_dict[param_tensor].shape !=  model.state_dict()[param_tensor].shape:
            trained_dict[param_tensor] = torch.squeeze(trained_dict[param_tensor], dim=1)

    model.load_state_dict(trained_dict)

    prediction = model([coords,ft])
    pred_np = prediction.cpu().detach().numpy()
    pred_np = pred_np[:,1] - pred_np[:,0]

    pred_coord = coords_np[np.argmax(pred_np)]
    if verbose :
        print('raw: pred_coord: ', pred_coord)

    pred_coord = pred_coord.astype(dtype)
    pred_coord *= resolution
    pred_coord += coords_offset + 0.5*resolution
    if verbose :
        print('final: pred_coord', pred_coord)

    return pred_coord.tobytes()
