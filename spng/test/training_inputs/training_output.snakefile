import torch
import scripts.roi_metrics as roi_metrics
import numpy as np

def get_with_chan_range(y, labels, chan_range):
    res  = roi_metrics.roi_metrics(
        y[chan_range[0]:chan_range[1]],
        labels[chan_range[0]:chan_range[1]]
    )
    return res

def load_y_labels(filename):
    t = torch.load(filename)
    labels = t['labels'][0].detach().cpu()
    y = t['y'][0].detach().cpu()

    return y, labels

chan_ranges = {
    'u':(0,800),
    'v':(800,1600),
    'w':(1600, 2560),
    'w0':(1600, 2080),
    'w1':(2080, 2560),
}

rule roi_eff_pur:
    run:
        y, labels = load_y_labels(input[0])
        res = get_with_chan_range(y, labels, chan_ranges[wildcards.plane])
        eff, pur = res['efficiency'].item(), res['purity'].item()
        print(eff,pur)
        import numpy as np
        np.savez(output[0], effs=np.array([eff]), purs=np.array([pur]))

use rule roi_eff_pur as angled_roi_eff_pur with:
    input:
        "<results>/xvu-test_line_{plane}plane_{angles}.pt"
    output:
        temp("<results>/xvu-roi_effs_purs_{plane}plane_{angles}.npz")

def load(filename, k='y'):
    t = torch.load(filename)
    return t[k][0].detach().cpu()

rule roi_table:
    params:
    run:
        chan_range = chan_ranges[wildcards.plane]
        t = load(input[0], k=('y' if wildcards.type=='reco' else 'labels'))[chan_range[0]:chan_range[1]]
        res = roi_metrics.roi_table(
            (t == 1 if wildcards.type=='true' else t > float(wildcards.threshold))
        )
        
        torch.save(res, output[0])

use rule roi_table as xvu_angled_roi_table with:
    input:
        "<results>/xvu-test_line_{plane}plane_{angles}.pt"
    output:
        "<results>/xvu-{type}_roi_table_{plane}plane_{angles}-thresh-{threshold}.pt",


rule merge_roi_eff_pur:
    run:
        import numpy as np
        effs, purs = [], []
        for f in input:
            t = np.load(f)
            effs += [t['effs']]
            purs += [t['purs']]
        effs = np.array(effs)
        purs = np.array(purs)

        np.savez(output[0], effs=effs, purs=purs)
vangles = [ 
    "t1-75-t2-75", "t1-80-t2-80", "t1-82-t2-82", "t1-85-t2-85", 
    "t1-75-t2-87", "t1-85-t2-87", "t1-87-t2-87",
]
uangles = [ 
    "t1-75-t2-75", "t1-80-t2-80", "t1-82-t2-82", "t1-85-t2-85", 
    "t1-87-t2-75", "t1-87-t2-85", "t1-87-t2-87",
]
wangles=uangles
use rule merge_roi_eff_pur as merge_vplane_roi_eff_purs with:
    input:
        expand(
            '<results>/xvu-roi_effs_purs_vplane_{angles}.npz',
            angles=vangles
        )
    output:
        "<results>/xvu-merged_roi_effs_purs_vplane_high_end.npz"

use rule merge_roi_eff_pur as merge_uplane_roi_eff_purs with:
    input:
        expand(
            '<results>/xvu-roi_effs_purs_uplane_{angles}.npz',
            angles=uangles
        )
    output:
        "<results>/xvu-merged_roi_effs_purs_uplane_high_end.npz"

use rule merge_roi_eff_pur as merge_wplane_roi_eff_purs with:
    input:
        expand(
            '<results>/xvu-roi_effs_purs_wplane_{angles}.npz',
            angles=wangles
        )
    output:
        "<results>/xvu-merged_roi_effs_purs_wplane_high_end.npz"


rule plot_merged_roi_eff_pur:
    params:
        ylim=(-.4, .4)
    run:
        import matplotlib.pyplot as plt
        import numpy as np
        xlabels = params.xlabels
        xpos = [i for i in range(len(xlabels))]

        t = np.load(input[0])
        all_effs = t['effs']
        all_purs = t['purs']
    

        plt.subplot(2, 1, 1)
        plt.plot(1.-all_effs, label='1-Efficiency')
        plt.xticks(xpos, xlabels)
        plt.grid(True, axis='both')
        plt.ylim(*params.ylim)
        plt.subplot(2, 1, 2)
        plt.plot(1.-all_purs)
        plt.grid(True, axis='both')
        plt.yticks(np.arange(*params.ylim, .1))
        plt.xticks(xpos, xlabels, label='1-Purity')
        plt.ylim(*params.ylim)
        plt.yticks(np.arange(*params.ylim, .1))
        f = output[0]
        print(type(f))
        print(f)
        plt.savefig(f)

high_end_angles_uplane=([75, 75], [80, 80], [82,82], [85,85], [87, 75], [87, 85], [87, 87])
high_end_angles_wplane=([75, 75], [80, 80], [82,82], [85,85], [87, 75], [87, 85], [87, 87])
high_end_angles_vplane=([75, 75], [80, 80], [82,82], [85,85], [75, 87], [85, 87], [87, 87])

def get_angles(wildcards):
    if wildcards.plane == 'u': return high_end_angles_uplane
    elif wildcards.plane == 'v': return high_end_angles_vplane
    elif wildcards.plane == 'w': return high_end_angles_wplane


use rule plot_merged_roi_eff_pur as plot_merged_roi_eff_pur_uplane with:
    params:
        xlabels=lambda w : [f'{u},{v}' for u,v in get_angles(w)]
    input:
        "<results>/xvu-merged_roi_effs_purs_{plane}plane_high_end.npz"
    output:
        '<results>/xvu-all_roi_eff_pur_{plane}plane.png'

rule all_roi_tables:
    input:
        (
            expand("<results>/xvu-{type}_roi_table_uplane_{angles}-thresh-0.5.pt",
                type=['true','reco'], angles=uangles) +
            expand("<results>/xvu-{type}_roi_table_vplane_{angles}-thresh-0.5.pt",
                type=['true', 'reco'], angles=vangles) +
            expand("<results>/xvu-{type}_roi_table_wplane_{angles}-thresh-0.5.pt",
                type=['true', 'reco'], angles=wangles)
        )

rule roi_lengths:
    run:
        nrois = []
        lengths = []
        names = []
        for f in input:
            t = torch.load(f)
            lengths += t['length']
            nrois.append(len(t['length']))
            names.append(f)
        
        np.savez(output[0], lengths=lengths, nrois=nrois, names=names)

use rule roi_lengths as true_u_roi_lengths with:
    input:
        expand("<results>/xvu-{{type}}_roi_table_uplane_{angles}-thresh-0.5.pt", angles=uangles)
    output: "<results>/xvu-{type}_roi_table_info_uplane-thresh-0.5.npz"

use rule roi_lengths as true_v_roi_lengths with:
    input:
        expand("<results>/xvu-{{type}}_roi_table_vplane_{angles}-thresh-0.5.pt", angles=vangles)
    output: "<results>/xvu-{type}_roi_table_info_vplane-thresh-0.5.npz"

use rule roi_lengths as true_w_roi_lengths with:
    input:
        expand("<results>/xvu-{{type}}_roi_table_wplane_{angles}-thresh-0.5.pt", angles=wangles)
    output: "<results>/xvu-{type}_roi_table_info_wplane-thresh-0.5.npz"