import torch
import numpy as numpy
import matplotlib.pyplot as plt
import os,sys
# Usage : alias tpython=/nfs/data/1/abashyal/spng/spng_dev_08252025/.direnv/python-3.11.9/bin/python
# Usage: tpython plot_torch_tensors.py file1.pt file2.pt ...
# create a list of pt files with the path from a given path
files = os.listdir(sys.argv[1] if len(sys.argv)>1 else '.')
#files with path
pt_files = [os.path.join(sys.argv[1], f) for f in files if f.endswith('.pt')]
print(f"Found {len(pt_files)} .pt files in the current directory.")
for tensor_file in pt_files:
    tensor = torch.jit.load(tensor_file)
    print(f"Loaded tensor from {tensor_file}")
    ## get the param
    param = next(tensor.parameters())
    print(f"Parameter from {tensor_file} has shape {param.shape} and dtype {param.dtype}")
    # if arr in gpu, move to cpu
    if param.is_cuda:
        param = param.cpu()
    arr = param.detach().numpy()
    print(f"Numpy array shape: {arr.shape}, dtype: {arr.dtype}")
    # if the shape is 4D, squeeze to 2D for plotting
    if arr.ndim == 4:
        #if the first two dims are 1 squeeze them
        if arr.shape[0] == 1 and arr.shape[1] == 1:
            arr = arr.squeeze()
        else:
            continue
    if arr.ndim==3:
        if arr.shape[0]==1:
            arr = arr.squeeze()
        else:
            continue
        print(f"Squeezed array shape to {arr.shape} for plotting")
    # Make sure that the dim 0 is always channel dim (i.e., number of wires)
    if arr.shape[0] > arr.shape[1]:
        arr = arr.T
        print(f"Transposed array shape to {arr.shape} for plotting")
    # Plot the array as an image
    plt.imshow(arr, aspect='auto', cmap='viridis')
    plt.colorbar()
    plt.title(f'Tensor from {tensor_file}')
    plt.xlabel('Ticks')
    plt.ylabel('Channels')
    # save as pnt file
    plt.savefig(f'{tensor_file}_plot.png')
    plt.clf()
    print(f"Plot saved to {tensor_file}_plot.png")
