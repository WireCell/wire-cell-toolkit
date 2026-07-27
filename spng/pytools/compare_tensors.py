import torch
import numpy as np
import matplotlib.pyplot as plt
import os,sys

# Compare the tensors in tow .pt files
# Usage : alias tpython=/nfs/data/1/abashyal/spng/spng_dev_08252025/.direnv/python-3.11.9/bin/python
# Usage: tpython compare_tensors.py file1.pt file2.pt
if len(sys.argv)<3:
    print("Usage: tpython compare_tensors.py file1.pt file2.pt")
    sys.exit(1)
file1 = sys.argv[1]
file2 = sys.argv[2]
outdir = sys.argv[3] if len(sys.argv)>3 else '.'
#if outdir does not exist, create it
if not os.path.exists(outdir):
    os.makedirs(outdir)
tensor1 = torch.jit.load(file1)
tensor2 = torch.jit.load(file2)
param1 = next(tensor1.parameters())
param2 = next(tensor2.parameters())
if param1.is_cuda:
    param1 = param1.cpu()
if param2.is_cuda:
    param2 = param2.cpu()
#print shape and dtypes of each param
print(f"Parameter from {file1} has shape {param1.shape} and dtype {param1.dtype}")
print(f"Parameter from {file2} has shape {param2.shape} and dtype {param2.dtype}")
arr1 = param1.detach().numpy()
arr2 = param2.detach().numpy()
# if the dimensions are more than 2D, squeeze to 2D for comparison
def squeeze_to_2d(arr):
    if arr.ndim == 4:
        if arr.shape[0] == 1 and arr.shape[1] == 1:
            arr = arr.squeeze()
    if arr.ndim == 3:
        if arr.shape[0] == 1:
            arr = arr.squeeze()
    return arr
arr1 = squeeze_to_2d(arr1)
arr2 = squeeze_to_2d(arr2)
#Make sure that the dimensions of each array aligns
if arr1.shape[0] > arr1.shape[1]:
    arr1 = arr1.T
if arr2.shape[0] > arr2.shape[1]:
    arr2 = arr2.T
# Now compare the two arrays
if arr1.shape != arr2.shape:
    print(f"Arrays have different shapes: {arr1.shape} vs {arr2.shape}")
    sys.exit(1)
diff = arr1 - arr2
mse = np.mean(diff**2)
print(f"Mean Squared Error between the two tensors: {mse}")
# Plot the difference
plt.imshow(diff, aspect='auto', cmap='bwr')
plt.colorbar()
plt.title(f'Difference between {file1} and {file2}')
plt.xlabel('Ticks')
plt.ylabel('Channels')
# save as png file
outfile = os.path.join(outdir, f'comparison_{os.path.basename(file1)})_{os.path.basename(file2)}.png')
plt.savefig(outfile)
plt.clf()
print(f"Difference plot saved to {outfile}")

#Compare the first 20 values of each array
print("First 20 values of the first tensor:")
print(arr1.flatten()[:20])
print("First 20 values of the second tensor:")
print(arr2.flatten()[:20])
print("First 20 values of the difference:")
print(diff.flatten()[:20])