
import torch
import torch.nn as nn
import sparseconvnet as scn

reps = 2 # Conv block repetition factor
m = 16 # Unet number of features
nPlanes = [m, 2*m, 4*m, 8*m, 16*m] # UNet number of features per level

class DeepVtx(nn.Module):
    '''
    spatialSize seems need to be 2^n
    '''
    def __init__(self,  dimension = 3, device = 'cuda', spatialSize = 4096, nIn = 3, nClasses = 2):
        nn.Module.__init__(self)
        self.sparseModel = scn.Sequential().add(
            scn.InputLayer(dimension, torch.LongTensor([spatialSize]*3), mode=3)).add(
            scn.SubmanifoldConvolution(dimension, nIn, m, 3, False)).add(
            scn.UNet(dimension, reps, nPlanes, residual_blocks=False, downsample=[2,2])).add(
            scn.BatchNormReLU(m)).add(
            scn.OutputLayer(dimension)).to(device)
        self.inputLayer = scn.InputLayer(dimension, torch.LongTensor([spatialSize]*3), mode=3)
        self.linear = nn.Linear(m, nClasses).to(device)
    def forward(self,x):
        x=self.sparseModel(x)
        x=self.linear(x)
        x=torch.sigmoid(x)
        return x
