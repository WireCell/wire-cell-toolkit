import torch
import matplotlib.pyplot as plt
import numpy as np


def draw_table(table_dict, shape=(800,1500), matched=None, verbose=False):
    frame = torch.zeros(shape, dtype=int)
    nrois = len(table_dict['row'])

    for roi in range(nrois):
        start = table_dict['start'][roi]
        roi_len = table_dict['length'][roi]
        row = table_dict['row'][roi]
        m = table_dict['matched'][roi]
        if matched is None: pass
        elif matched and not m: continue
        elif not matched and m: continue
        
        frame[row, start:start+roi_len] = 1
        mstr = 'matched' if m else 'unmatched'
        if verbose:
            print(
                f'ROI {roi} -- {mstr} -- row: {row}, start: {start}, len: {roi_len}'
            )


    return frame 

if __name__ == '__main__':
    import sys
    # plt.ion()
    table = torch.load(sys.argv[1])
    plt.subplot(3,1,1)
    plt.imshow(draw_table(table), cmap='bwr', vmin=-1, vmax=1, aspect='auto', interpolation='none')

    plt.subplot(3,1,2)
    plt.imshow(draw_table(table, matched=True, verbose=True), cmap='bwr', vmin=-1, vmax=1, aspect='auto', interpolation='none')

    print()

    plt.subplot(3,1,3)
    plt.imshow(draw_table(table, matched=False, verbose=True), cmap='bwr', vmin=-1, vmax=1, aspect='auto', interpolation='none')

    plt.show()