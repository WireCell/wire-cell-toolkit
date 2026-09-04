import torch
import scripts.roi_metrics as roi_metrics
import numpy as np

def get_with_chan_range(y, labels, chan_range):
    res  = roi_metrics.roi_metrics(
        y[chan_range[0]:chan_range[1]],
        labels[chan_range[0]:chan_range[1]]
    )
    return res

def load_y_labels(filename, device='cpu'):
    t = torch.load(filename)
    labels = t['labels'][0].detach().to(device)
    y = t['y'][0].detach().to(device)

    return y, labels

chan_ranges = {
    'u':(0,800),
    'v':(800,1600),
    'w':(1600, 2560),
    'w0':(1600, 2080),
    'w1':(2080, 2560),
    'all_planes':(0,2560),
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

config.setdefault('threshold_step', 0.1)
rule threshold_scan:
    input:
        "run_one_{entry}.pt"
    output:
        "threshold_scan_roi_plane_{plane}_run_one_{entry}.pt"
    params:
        step=config['threshold_step']
    run:
        thresholds=np.arange(0., 1., params.step)
        y, labels = load_y_labels(input[0], config['device'])
        chan_range = chan_ranges[wildcards.plane]
        res = roi_metrics.threshold_scan(y[chan_range[0]:chan_range[1]], labels[chan_range[0]:chan_range[1]], thresholds, as_eff_pur=False)
        torch.save(res, output[0])

rule threshold_pixel_scan:
    input:
        "run_one_{entry}.pt"
    output:
        "threshold_scan_pixel_plane_{plane}_run_one_{entry}.pt"
    params:
        step=config['threshold_step']
    run:
        thresholds=torch.arange(0., 1., params.step).to(config['device'])
        chan_range = chan_ranges[wildcards.plane]
        y, labels = load_y_labels(input[0], config['device'])
        y, labels = y[chan_range[0]:chan_range[1]], labels[chan_range[0]:chan_range[1]]

        n_true = labels.sum()
        reco = (y.unsqueeze(-1).expand(-1,-1,len(thresholds)) > thresholds)
        n_reco = reco.sum((0,1))
        matched = ((labels>0).unsqueeze(-1).expand(-1,-1,reco.shape[-1]) == reco) * (reco > 0)

        n_reco_matched = matched.sum((0,1))
        n_true_matched = n_reco_matched #for compatibility with the other scan+aggregation
        torch.save({
            "threshold": thresholds,
            "n_reco": n_reco,
            "n_reco_matched": n_reco_matched,
            "n_true_matched": n_true_matched,
            "n_true": n_true,
            "efficiency": (n_true_matched/n_true),
            "purity": (n_reco_matched/n_reco),
        }, output[0])

config.setdefault('nentries', 1)
rule aggregate_scan:
    input:
        collect('threshold_scan_{{type}}_plane_{{plane}}_run_one_{entry}.pt', entry=[i for i in range(config['nentries'])])
    output:
        'aggregated_scan_{type}_plane_{plane}.npz'
    run:

        effs = []
        purs = []
        n_trues = []
        n_recos = []
        n_trues_matched = []
        n_recos_matched = []
        thresholds = [] #Should all be the same

        for input_file in input:
            t = torch.load(input_file)
            effs.append(t['efficiency'].cpu().numpy())
            purs.append(t['purity'].cpu().numpy())
            n_trues.append(t['n_true'].cpu().numpy())
            n_recos.append(t['n_reco'].cpu().numpy())
            n_recos_matched.append(t['n_reco_matched'].cpu().numpy())
            n_trues_matched.append(t['n_true_matched'].cpu().numpy())
            thresholds.append(t['threshold'].cpu().numpy())
        
        np.savez(
            output[0],
            effs=effs,
            purs=purs,
            n_trues=n_trues,
            n_recos=n_recos,
            n_recos_matched=n_recos_matched,
            n_trues_matched=n_trues_matched,
            n_recos_summed=np.array(n_recos).sum(axis=0),
            n_trues_summed=sum(n_trues),
            n_recos_matched_summed=np.array(n_recos_matched).sum(axis=0),
            n_trues_matched_summed=np.array(n_trues_matched).sum(axis=0),
            thresholds=thresholds[0]
        )
config.setdefault('device', 'cpu')
config.setdefault('app', 'xvunet')
config.setdefault('nevents', 10)
rule run_n:
    resources:
        gpu = 1 if ('gpu' in config['device'] or 'cuda' in config['device']) else 0
    output:
        "run_one_{entry}.pt"
    params:
        paths=config['paths'],
        model_file=config['model_file'],
        device=config['device'],
        config=config['cfg'],
        app=config['app'],
        nevents=config['nevents'],
    shell:
        """
        wcpy dnn run_one -l {params.model_file} -c {params.config} --manual-sigmoid -d {params.device} \
            -a {params.app} -o {output} -n {wildcards.entry} {params.paths}
        """




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