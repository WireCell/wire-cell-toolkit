def get_with_chan_range(y, labels, chan_range):
    import scripts.roi_metrics as roi_metrics
    res  = roi_metrics.roi_metrics(
        y[chan_range[0]:chan_range[1]],
        labels[chan_range[0]:chan_range[1]]
    )
    return res

def load_y_labels(filename):
    import torch
    t = torch.load(filename)
    labels = t['labels'][0].detach().cpu()
    y = t['y'][0].detach().cpu()

    return y, labels

rule roi_eff_pur:
    params:
        chan_range=(0,800)
    run:

        y, labels = load_y_labels(input[0])
        res = get_with_chan_range(y, labels, params.chan_range)
        eff, pur = res['efficiency'].item(), res['purity'].item()
        print(eff,pur)
        import numpy as np
        np.savez(output[0], effs=np.array([eff]), purs=np.array([pur]))

use rule roi_eff_pur as angled_roi_eff_pur with:
    input:
        "test_line_{plane}plane_{angles}.pt"
    output:
        temp("roi_effs_purs_{plane}plane_{angles}.npz")


rule roi_tables:
    params:
        chan_range=(0,800)
    run:
        y, labels = load_y_labels(input[0])
        res = get_with_chan_range(y, labels, params.chan_range)
        true_t, reco_t = res['true_t'], res['reco_t']
        
        import torch
        torch.save(true_t, output.true)
        torch.save(reco_t, output.reco)

use rule roi_tables as angled_roi_tables with:
    input:
        "test_line_vplane_{angles}.pt"
    output:
        true="true_roi_table_vplane_{angles}.pt",
        reco="reco_roi_table_vplane_{angles}.pt"

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

use rule merge_roi_eff_pur as merge_vplane_roi_eff_purs with:
    input:
        "roi_effs_purs_vplane_t1-75-t2-75.npz",
        "roi_effs_purs_vplane_t1-80-t2-80.npz",
        "roi_effs_purs_vplane_t1-82-t2-82.npz",
        "roi_effs_purs_vplane_t1-85-t2-85.npz",
        "roi_effs_purs_vplane_t1-75-t2-87.npz",
        "roi_effs_purs_vplane_t1-85-t2-87.npz",
        "roi_effs_purs_vplane_t1-87-t2-87.npz",
    output:
        "merged_roi_effs_purs_vplane_high_end.npz"

use rule merge_roi_eff_pur as merge_uplane_roi_eff_purs with:
    input:
        "roi_effs_purs_uplane_t1-75-t2-75.npz",
        "roi_effs_purs_uplane_t1-80-t2-80.npz",
        "roi_effs_purs_uplane_t1-82-t2-82.npz",
        "roi_effs_purs_uplane_t1-85-t2-85.npz",
        "roi_effs_purs_uplane_t1-87-t2-75.npz",
        "roi_effs_purs_uplane_t1-87-t2-85.npz",
        "roi_effs_purs_uplane_t1-87-t2-87.npz",
    output:
        "merged_roi_effs_purs_uplane_high_end.npz"

use rule merge_roi_eff_pur as merge_wplane_roi_eff_purs with:
    input:
        "roi_effs_purs_wplane_t1-75-t2-75.npz",
        "roi_effs_purs_wplane_t1-80-t2-80.npz",
        "roi_effs_purs_wplane_t1-82-t2-82.npz",
        "roi_effs_purs_wplane_t1-85-t2-85.npz",
        "roi_effs_purs_wplane_t1-87-t2-75.npz",
        "roi_effs_purs_wplane_t1-87-t2-85.npz",
        "roi_effs_purs_wplane_t1-87-t2-87.npz",
    output:
        "merged_roi_effs_purs_wplane_high_end.npz"


rule plot_merged_roi_eff_pur:
    params:
        ylim=(0, 1.5)
    run:
        import matplotlib.pyplot as plt
        import numpy as np
        xlabels = params.xlabels
        xpos = [i for i in range(len(xlabels))]

        t = np.load(input[0])
        all_effs = t['effs']
        all_purs = t['purs']
    

        plt.subplot(2, 1, 1)
        plt.plot(all_effs, label='Efficiency')
        plt.xticks(xpos, xlabels)
        plt.grid(True, axis='both')
        plt.ylim(*params.ylim)
        plt.subplot(2, 1, 2)
        plt.plot(all_purs)
        plt.grid(True, axis='both')
        plt.yticks(np.arange(0.,1.5, .2))
        plt.xticks(xpos, xlabels, label='Purity')
        plt.ylim(*params.ylim)
        plt.yticks(np.arange(0.,1.5, .2))
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
        "merged_roi_effs_purs_{plane}plane_high_end.npz"
    output:
        'all_roi_eff_pur_{plane}plane.png'