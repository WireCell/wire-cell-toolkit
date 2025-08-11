#!/usr/bin/env python
'''
A module to plot data in "tiling file format".

This file format is defined in and produced by check_generate_tiling.cxx.

It is a pickle file with a dictionary of tensors.

Note, after load_file() everything is a Numpy array, not a torch Tensor.

Example usage:

  $ ./build/spng/check_generate_tiling check_generate_tiling.pt 5

  $ python spng/test/tiling.py check_generate_tiling.{pt,pdf}

'''

# limit what we use from torch
from torch import Tensor as torch_Tensor
from torch import load as torch_load

# the rest is numpy
import numpy
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Polygon
def numpyify(toa):
    if isinstance(toa, dict):
        return {k:numpyify(v) for k,v in toa.items()}
    if isinstance(toa, tuple):
        return tuple(numpyify(list(toa)))
    if isinstance(toa,tuple):
        return [numpyify(one) for one in toa]
    if isinstance(toa, torch_Tensor):
        return toa.numpy()
    return toa


def load_file(tiling_file):
    '''
    Load a pickled .pt file and return a dict of Numpy arrays.

    Some keys:

    The arrays from the Coordinates object:

    - coords_pitch_mag: (Nview,) magnitude of pitch of each view
    - coords_pitch_dir: (Nview, 2) unit vector along the pitch direction of each view
    - coords_center: (Nview, 2) origin vector of each view
    - coords_zero_crossings: (Nview, Nview, 2) crossing point of a "ray zero" from a pair of views.
    - coords_ray_jump: (Nview, Nview, 2) displacement vector along ray direction between crossings
    - coords_ray_dir: (Nview, 2) direction of ray (wire) for view.
    - coords_a: The ray grid "A" tensor.  See raygrid.pdf doc.
    - coords_b: The ray grid "B" tensor.  See raygrid.pdf doc.

    The other arrays:

    - views : (5,2,2) the pitch rays for the 5 views.
    - points: (Ngroups, NptsPerGroup, 2) groups of 2D points.
    - activity_V: One in a list of size Nview of 1D activity tensors for each view V in [0,1,2,3,4].
    - blobs : (Nblobs, Nlayers=5, 2) blobs as 5 strips, each with 2 wire indices
    - crossings : (Nblobs, Npairs, 4, 2) the 4 pairs of strip edges crossings from Nlayer-choose-2 strips
    - insides : (Nblobs, Npairs, 4) boolean for which stip edges are "inside" all strips.

    '''
    ret = numpyify(torch_load(tiling_file, map_location='cpu', weights_only=False))
    for k,v in ret.items():
        print(f'{k}: {v.shape}')
    return ret


def get_activities(data):
    '''
    Return ordered tuple of activity tensors in the data.
    '''
    return tuple([data[f'activity_{layer}'] for layer in range(5)])


def threshold_1d(activity, minimum=0.0):
    '''
    Return pairs of indices of a 1D array of activity measured across one
    view.  Each pair of indices gives the begin/end bounds of a half-open range
    that has contiguous activity values above the minimum.

    Returns an array of shape (2,N) for N ranges.  The [0] element is the start
    indices for the half-open bounds and the [1] element is the end indices.
    '''
    mask = activity > minimum
    mask = numpy.hstack((numpy.array([False]), mask, numpy.array([False])))
    begs = numpy.array(((mask[1:]  == True) & (mask[:-1] == False)).nonzero()).squeeze()
    ends = numpy.array(((mask[:-1] == True) & (mask[1:]  == False)).nonzero()).squeeze()
    return numpy.vstack((begs, ends))


class Drawing:
    def __init__(self, file_or_data):
        if not isinstance(file_or_data, dict):
            file_or_data = load_file(file_or_data)
        for name, arr in file_or_data.items():
            setattr(self, name, arr)

        # normalize some arrays

        # fold activities into tuple
        self.activities = get_activities(file_or_data)
        # find lo/hi bounds for activities
        self.activity_bounds = tuple([threshold_1d(act) for act in self.activities])

        if len(self.points.shape) == 3:
            self.grouped_points = self.points
            self.points = self.points.reshape(-1, 2)
        else:
            self.grouped_points = self.points

        self.pitch_vec = self.coords_pitch_mag[:,None]*self.coords_pitch_dir

        # The canonical ordering for edge pairs as one runs over the dimension
        # of size 4 of crossings and insides arrays.
        self.strip_pair_edge_indices = numpy.array([ [0,0], [0,1], [1,0], [1,1] ])

        self.view_colors = ['grey', 'grey', 'r', 'g', 'b']

    # fixme: find corners, draw outline.



    def draw_points(self, color='k'):
        '''
        Draw the points.
        '''
        plt.scatter(self.points[:,0], self.points[:,1], s=10, c=color)
    

    def draw_strip_edges(self, lo, hi, view, scale=1000.0):
        '''
        Draw a strip edges with given wire-index bounds in the given view
        '''
        origin = self.coords_center[view]
        pvec = self.pitch_vec[view]

        # points on strip edges
        plo = origin + lo*pvec
        phi = origin + hi*pvec
        
        # A direction along the wire
        rdir = self.coords_ray_dir[view];

        # Extend to be "long enough"
        rvec = scale * rdir

        color = self.view_colors[view]
        plt.axline(plo-rvec, plo+rvec, c=color, linestyle='solid')
        plt.axline(phi-rvec, phi+rvec, c=color, linestyle='dashed')


    def draw_strip_body(self, lo, hi, view, scale=4000.0):
        '''
        Draw a strip by filling with given wire-index bounds in the given view
        '''
        origin = self.coords_center[view]
        pvec = self.pitch_vec[view]

        # points on strip edges
        plo = origin + lo*pvec
        phi = origin + hi*pvec
        
        # A direction along the wire
        rdir = self.coords_ray_dir[view];

        # Extend to be "long enough"
        rvec = scale * rdir

        ul = plo+rvec
        ll = plo-rvec
        ur = phi+rvec
        lr = phi-rvec
        verts = [ul, ur, lr, ll]

        color = self.view_colors[view]
        strip = Polygon(verts, closed=True, facecolor=color, alpha=0.3, edgecolor='none')

        plt.gca().add_patch(strip)


    def draw_blob_strips(self):
        '''
        Draw the blob strips as view rays
        '''

        for blob in self.blobs:
            for iview in range(2,5):
                self.draw_strip_edges(*blob[iview], iview)

        return

    def draw_activity_strips(self):
        '''
        Draw the activity strips.
        '''
        colors = ['k', 'grey', 'r', 'c', 'm']
        for iview in range(2, 5):
            bounds = self.activity_bounds[iview]
            for lo,hi in bounds.T:
                self.draw_strip_body(lo, hi, iview)

        return

    def draw_many(self, xlim=None, ylim=None):
        self.draw_points()
        self.draw_blob_strips()
        self.draw_activity_strips()
        if xlim is not None:
            plt.xlim(*xlim)
        if ylim is not None:
            plt.ylim(*ylim)


if "__main__" == __name__:
    import sys
    d = Drawing(sys.argv[1])
    #d.draw_many((1000,1200), (1000,1200))
    #d.draw_many((200,600), (0, 600))
    d.draw_many()

    #d.draw_many((3500, 4000), (1500, 1800)) # 4 40
    # d.draw_many((3300, 3400), (2800, 3000)) # 4 37
    print(f'writing {sys.argv[2]}')
    plt.savefig(sys.argv[2])
