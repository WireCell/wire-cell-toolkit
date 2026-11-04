# -*- snakemake -*-
#
# Evaluate trained DNN-ROI models.  Unifies the evaluation half of ./Snakefile
# with ./training_output.snakefile.  Everything upstream of the h5 files
# (LArSoft, depo generation, wire-cell sigproc, tar -> h5) is out of scope:
# the line-sample h5 files must already exist (see config input_dir).

USAGE = '''
Usage:
  snakemake -s evaluate.snakefile --cores N --directory WORKDIR \\
      --config model=CHECKPOINT cfg=TRAIN_CFG [app=xvunet] [device=cpu|gpu] [KEY=VAL ...] \\
      TARGET ...

Config (put per-site paths in a --configfile):
  device          cpu (default), gpu or cuda
  app             wcpy dnn app (default xvunet; legacy alias: dnnapp)
  model           checkpoint for run_one / run_n (legacy alias: model_file);
                  also the xvunet model tagged 'xvu'
  cfg             training config passed to wcpy dnn run_one
  xvu_models      more xvunet models to evaluate side by side, as
                  {TAG: {model: ..., cfg: ..., app: ...}}; each TAG is a style
                  (letters, digits, _), with outputs results/TAG-test_{sample}.pt
  input_dir       directory holding the linedepos-*.h5 files (default: WORKDIR)
  lineseeds       seeds of the line samples, i.e. files linedepos-...-seedS-g4-*.h5
                  (Snakefile lineseeds); unset for unseeded files
  cosmics_dir     directory of N_M/ subdirectories holding cosmics_*-g4-*.h5
  cosmics_dirs    which N_M subdirectories to use (default: all of them)
  cosmics_events  events per cosmics file (default 10)
  view_app        per-view app run by run_view (default dnnroi_custom)
  view_models     per-view [{model: ..., cfg: ...}] for views 0, 1, 2
  view_widths     channels per view in the merged output (default [800, 800, 960])
  dnnroi_dir      where per-view and merged outputs go (default dnnroi)
  compare_styles  styles compared by compare_trio_rates (default [xvu, merged])
  cosmics_tune_dirs, cosmics_eval_dirs
                  cosmics directories to tune operating points on / evaluate
                  them on (default: first half / the rest)
  report_sets     sets op_report collates (default [cosmics_eval, uplane_high_end,
                  vplane_high_end, wplane_high_end])
  op_compare      [style, mode] pairs compared by compare_op_trio_rates, mode
                  perplane (best threshold per u, v, w0) or single (one over
                  u+v+w0) (default [[xvu, perplane], [xvu, single], [merged, perplane]])
  pathvars        e.g. {results: results_0_634} to keep each model's results apart
  test_tpcs       views scored by plot_eff_pur (default [0])
  cpt_dir, proc, cluster
                  locate checkpoint_{proc}_{cluster}_{epoch}.pt for the epoch_* rules
  paths, nentries input files and number of entries for run_n and the threshold scans
  threshold_start, threshold_step
                  first threshold and step of the threshold scans (default 0.1, 0.1)

Samples and styles:
  A {sample} is line_Pplane_t1-A-t2-B (one line-depo event) or cosmics_N_M_evE
  (event E of cosmics directory N_M).  The model output for a sample is
    xvu     results/xvu-test_{sample}.pt      model/cfg/app run on all three views
    merged  dnnroi/merged-test_{sample}.pt    view_models run per view, then
                                              concatenated (and zero-padded) to
                                              the same channel layout as xvu
    noxvu   results/test_{sample}.pt          (line samples only) model run on
                                              the sample's own view

Targets (rule names):
  all_plots_{u,v,w}[_xvu]   event displays + pixel eff/pur vs angle
  all_eff_pur_xvu           pixel eff/pur vs angle, all planes
  all_epochs_xvu            pixel eff/pur vs angle for every checkpoint in cpt_dir
  all_roi_tables            true/reco ROI tables at threshold 0.5
  all_trio_compare          compare_styles trio rates, cosmics and line samples
  all_trio_compare_op       op_compare trio rates at best ROI / pixel F1 thresholds
Targets (files):
  results/threshold-T-xvu-all_{roi,pixel}_eff_pur_Pplane.png
  results/xvu-{true,reco}_roi_table_info_Pplane-thresh-0.5.npz
  aggregated_scan_{roi,pixel}_plane_P_fbetas.png
  results/threshold-T-{xvu,noxvu,merged}-trio-merged-crossplane_Pplane_high_end.png
  results/threshold-T-{xvu,merged}-trio-merged-crossplane_cosmics.npz
  results/threshold-T-trio-compare_{cosmics,Pplane_high_end}.png
  results/op-{roi,pixel}F1-trio-compare_{cosmics_eval,Pplane_high_end}.png
  results/op-{roi,pixel}F1_{perplane,single}-{xvu,merged}.json   (thresholds)
  results/eff_pur_epochs_Pplane[-xvu].pdf
'''

import json
import os
import re
import sys
import warnings
from glob import glob

import numpy as np
import torch
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, workflow.basedir)
import scripts.roi_metrics as roi_metrics

# Force TrueType fonts and silence the missing-glyph warning that crashes
# some matplotlib builds.
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['font.family'] = 'DejaVu Sans'
warnings.filterwarnings("ignore", category=UserWarning, module="matplotlib")


#
# Configuration
#
if 'dnnapp' in config:
    config.setdefault('app', config['dnnapp'])
if 'model_file' in config:
    config.setdefault('model', config['model_file'])
config.setdefault('device', 'cpu')
config.setdefault('app', 'xvunet')
config.setdefault('model', '')
config.setdefault('cfg', '')
config.setdefault('input_dir', '')
config.setdefault('lineseeds', None)
config.setdefault('test_tpcs', [0])
config.setdefault('cpt_dir', '.')
config.setdefault('proc', '0')
config.setdefault('cluster', '0')
config.setdefault('paths', '')
config.setdefault('nentries', 1)
config.setdefault('threshold_step', 0.1)
# At threshold 0 every pixel passes, so each channel is one ROI holding some
# truth and the ROI scan reports R = P = 1.
config.setdefault('threshold_start', 0.1)
config.setdefault('cosmics_dir', '')
config.setdefault('cosmics_dirs', None)
config.setdefault('cosmics_events', 10)
config.setdefault('view_app', 'dnnroi_custom')
config.setdefault('view_models', [])
config.setdefault('view_widths', [800, 800, 960])
config.setdefault('dnnroi_dir', 'dnnroi')
config.setdefault('compare_styles', ['xvu', 'merged'])
config.setdefault('cosmics_tune_dirs', None)
config.setdefault('cosmics_eval_dirs', None)
config.setdefault('xvu_models', {})
config.setdefault('report_sets', ['cosmics_eval'] + [f'{p}plane_high_end' for p in 'uvw'])
config.setdefault('op_compare',[['xvu', 'perplane'], ['xvu', 'single'], ['merged', 'perplane']])

TORCH_DEVICE = 'cuda' if 'gpu' in config['device'] else config['device']
USE_GPU = 1 if ('gpu' in config['device'] or 'cuda' in config['device']) else 0


#
# Detector and sample definitions
#
DET = 'pdhd'
PLANES = ['u', 'v', 'w']
PLANE_TO_VIEW = {'u': 0, 'v': 1, 'w': 2}

# Channel ranges of one PDHD APA in the (channel, tick) images.
CHAN_RANGES = {
    'u': (0, 800),
    'v': (800, 1600),
    'w': (1600, 2560),
    'w0': (1600, 2080),
    'w1': (2080, 2560),
    'all_planes': (0, 2560),
    'uvw0': (0, 2080),          # u, v and the w face with simulated signal
}

# Ranges each scan_sample covers.
SCAN_PLANES = ['u', 'v', 'w0', 'uvw0']

# Operating points: per view (u, v, w), the SCAN_PLANES range whose best-F1
# threshold that view uses.
OP_PLANES = {
    'perplane': ['u', 'v', 'w0'],
    'single': ['uvw0', 'uvw0', 'uvw0'],
}
OP_TUNE_SET = 'cosmics_tune'

# High-end (theta_1, theta_2) line samples per plane, in degrees.
ANGLES = {
    'u': [(75, 75), (80, 80), (82, 82), (85, 85), (87, 75), (87, 85), (87, 87)],
    'v': [(75, 75), (80, 80), (82, 82), (85, 85), (75, 87), (85, 87), (87, 87)],
    'w': [(75, 75), (80, 80), (82, 82), (85, 85), (87, 75), (87, 85), (87, 87)],
}

FBETA_BETAS = [0.50, 0.75, 1.00, 1.25, 1.50, 2.00]

# Trio prediction patterns (U, V, W) counted by trios_cross_plane.
TRIO_COMBOS = {
    'uvw': (True, True, True),
    'uv': (True, True, False),
    'vw': (False, True, True),
    'wu': (True, False, True),
    'just_u': (True, False, False),
    'just_v': (False, True, False),
    'just_w': (False, False, True),
}

LINE_SAMPLE_RE = r'line_(?P<plane>[uvw])plane_t1-(?P<t1>\d+)-t2-(?P<t2>\d+)(?:_seed(?P<seed>\d+))?'
COSMICS_SAMPLE_RE = r'cosmics_(?P<dir>\d+_\d+)_ev(?P<event>\d+)'
# The same without named groups, for wildcard constraints.
LINE_SAMPLE = r'line_[uvw]plane_t1-\d+-t2-\d+(?:_seed\d+)?'
COSMICS_SAMPLE = r'cosmics_\d+_\d+_ev\d+'

# Seeds of the line samples (Snakefile lineseeds); None for unseeded files.
LINE_SEEDS = [str(s) for s in config['lineseeds']] if config['lineseeds'] else [None]

# Where each style's model output for a sample lives: STYLE_PREFIX + 'test_{sample}.pt'.
# xvunet-style models (all three views in, the whole APA out), by tag.  The
# top-level model/cfg/app keys give the model tagged 'xvu'.
XVU_MODELS = {}
if config['model']:
    XVU_MODELS['xvu'] = {'model': config['model'], 'cfg': config['cfg']}
XVU_MODELS.update(config['xvu_models'] or {})
for tag, spec in XVU_MODELS.items():
    if not re.fullmatch(r'[A-Za-z0-9_]+', tag) or tag in ('noxvu', 'merged'):
        raise ValueError(f'bad xvu_models tag {tag!r}: use letters, digits and _, not noxvu/merged')
    spec.setdefault('app', config['app'])
XVU_STYLES = '|'.join(XVU_MODELS) or 'xvu'

STYLE_PREFIX = {
    **{tag: f'<results>/{tag}-' for tag in XVU_MODELS},
    'noxvu': '<results>/',
    'merged': os.path.join(config['dnnroi_dir'], 'merged-'),
}

wildcard_constraints:
    plane='|'.join(CHAN_RANGES),
    t1=r'\d+',
    t2=r'\d+',
    angles=r't1-\d+-t2-\d+',
    threshold=r'[0-9.]+',
    beta=r'[0-9.]+',
    style='|'.join(STYLE_PREFIX),
    sample=f'{LINE_SAMPLE}|{COSMICS_SAMPLE}',
    view='[012]',
    dir=r'\d+_\d+',
    set=r'cosmics|cosmics_tune|cosmics_eval|[uvw]plane_high_end',
    mode='|'.join(OP_PLANES),
    type='roi|pixel|true|reco',
    epoch=r'\d+',
    entry=r'\d+',


#
# Helpers
#
def angle_tag(t1, t2):
    'Angle naming used by the run_one results, e.g. t1-75-t2-87.'
    return f't1-{t1}-t2-{t2}'

def depo_tag(t1, t2):
    'Angle naming used by the line-sample h5 files, e.g. t1_75-t2_87.'
    return f't1_{t1}-t2_{t2}'

def angle_tags(plane):
    return [angle_tag(t1, t2) for t1, t2 in ANGLES[plane]]

def angle_labels(plane):
    return [f'{t1},{t2}' for t1, t2 in ANGLES[plane]]

def h5_path(name):
    return os.path.join(config['input_dir'], name) if config['input_dir'] else name

def line_sample(plane, t1, t2, seed=None):
    return f'line_{plane}plane_{angle_tag(t1, t2)}' + (f'_seed{seed}' if seed is not None else '')

def parse_sample(sample):
    'Return (kind, fields) for a {sample} wildcard value.'
    for kind, pat in (('line', LINE_SAMPLE_RE), ('cosmics', COSMICS_SAMPLE_RE)):
        m = re.fullmatch(pat, sample)
        if m:
            return kind, m.groupdict()
    raise ValueError(f'unknown sample {sample!r}')

def sample_h5_base(sample):
    'Path of the sample\'s h5 files, up to and including "-g4-".'
    kind, s = parse_sample(sample)
    if kind == 'line':
        seed = f'-seed{s["seed"]}' if s['seed'] is not None else ''
        return h5_path(f'linedepos-{DET}-{s["plane"]}plane-{depo_tag(s["t1"], s["t2"])}{seed}-g4-')
    found = glob(os.path.join(config['cosmics_dir'], s['dir'], 'cosmics_*-g4-rec-0.h5'))
    if len(found) != 1:
        raise ValueError(f'expected one cosmics_*-g4-rec-0.h5 for {sample}, found {found}')
    return found[0][:-len('rec-0.h5')]

def sample_h5(sample, kind, views):
    'The sample\'s rec or tru h5 files for the given views.'
    base = sample_h5_base(sample)
    return [f'{base}{kind}-{v}.h5' for v in views]

def trio_group(sample):
    'Group of the sample\'s trio h5 file: the event for cosmics.'
    kind, s = parse_sample(sample)
    return s['event'] if kind == 'cosmics' else '0'

COSMICS_EVENTS = range(int(config['cosmics_events']))

def event_entry(event):
    '''
    Dataset entry of an event group.  wcpy sorts samples by the group name
    captured by path_res, as a string, so event 10 comes before event 2.
    '''
    return sorted(str(e) for e in COSMICS_EVENTS).index(str(event))

def cosmics_dir_names():
    dirs = config['cosmics_dirs']
    if dirs is None:
        dirs = sorted(d for d in os.listdir(config['cosmics_dir']) if re.fullmatch(r'\d+_\d+', d))
    return dirs

def cosmics_split():
    '''
    Directories of the cosmics_tune set (operating-point thresholds are tuned
    on it) and of the held-out cosmics_eval set.  Default: first half / rest.
    '''
    dirs = cosmics_dir_names()
    tune = config['cosmics_tune_dirs'] or dirs[:len(dirs) // 2]
    held_out = config['cosmics_eval_dirs'] or [d for d in dirs if d not in tune]
    overlap = set(tune) & set(held_out)
    if overlap:
        raise ValueError(f'cosmics directories in both tune and eval sets: {sorted(overlap)}')
    return {'cosmics_tune': tune, 'cosmics_eval': held_out}

def cosmics_samples(dirs=None):
    if dirs is None:
        dirs = cosmics_dir_names()
    return [f'cosmics_{d}_ev{e}' for d in dirs for e in COSMICS_EVENTS]

def run_n_events(rec, tru, app, model, cfg, outputs):
    '''
    Run every event of one cosmics directory through a model in a single wcpy
    process (startup dominates the run time) and name the outputs by event.
    outputs are in event order.
    '''
    by_entry = outputs[0].replace('_ev0.pt', '_entry{entry}.pt')
    device = TORCH_DEVICE
    n = len(outputs)
    shell("wcpy dnn run_n {rec} {tru} -a {app} -l {model} -c {cfg} --manual-sigmoid "
          "-d {device} -n {n} -o '{by_entry}'")
    for event, out in enumerate(outputs):
        os.replace(by_entry.replace('{entry}', str(event_entry(event))), out)

def line_file_re(kind, view):
    '''
    File regex for a run_one dataset over many line samples: the ID group is
    the whole sample part of the name (plane, angles and seed), so samples do
    not collide the way the default ID (digits before -g4) makes them.
    '''
    return rf'(?:.*/)?linedepos-{DET}-(.+)-g4-{kind}-{view}\.h5'

def line_id(sample):
    'The ID line_file_re captures for a line sample.'
    name = os.path.basename(sample_h5_base(sample))
    return name[len(f'linedepos-{DET}-'):-len('-g4-')]

def run_n_lines(rec, tru, app, model, cfg, dataset, samples, outputs):
    '''
    Run every line sample through a model in a single wcpy process and name
    the outputs by sample (outputs in samples order).  The model's cfg is
    copied with the [run_one_dataset] keys in dataset set, to select the
    samples by line_file_re.  wcpy orders samples by ID string.
    '''
    import configparser, shutil, tempfile
    tmp = tempfile.mkdtemp(prefix='.run_n-', dir=os.path.dirname(outputs[0]) or '.')
    parser = configparser.ConfigParser(interpolation=None, strict=False)
    parser.optionxform = str
    parser.read(cfg)
    if not parser.has_section('run_one_dataset'):
        parser.add_section('run_one_dataset')
    for key, value in dataset.items():
        parser.set('run_one_dataset', key, value)
    batch_cfg = os.path.join(tmp, 'run_one.cfg')
    with open(batch_cfg, 'w') as out:
        parser.write(out)

    by_entry = os.path.join(tmp, 'entry{entry}.pt')
    device = TORCH_DEVICE
    n = len(samples)
    shell("wcpy dnn run_n {rec} {tru} -a {app} -l {model} -c {batch_cfg} --manual-sigmoid "
          "-d {device} -n {n} -o '{by_entry}'")
    order = sorted(line_id(s) for s in samples)
    for sample, out in zip(samples, outputs):
        os.replace(by_entry.replace('{entry}', str(order.index(line_id(sample)))), out)
    shutil.rmtree(tmp)

def set_samples(name):
    'The samples of a set: cosmics, cosmics_tune, cosmics_eval or {plane}plane_high_end.'
    if name == 'cosmics':
        return cosmics_samples()
    if name.startswith('cosmics_'):
        return cosmics_samples(cosmics_split()[name])
    plane = name[0]
    return [line_sample(plane, t1, t2, seed) for t1, t2 in ANGLES[plane] for seed in LINE_SEEDS]

ALL_LINE_SAMPLES = [s for p in PLANES for s in set_samples(f'{p}plane_high_end')]

def sample_group(sample):
    '''
    (group, replica) of a sample: (angle "t1,t2", seed) for a line sample,
    (sample, sample) for anything else.
    '''
    kind, s = parse_sample(sample)
    if kind == 'line':
        return f'{s["t1"]},{s["t2"]}', s['seed'] or '0'
    return sample, sample

def summed_and_spread(k, n, groups, replicas):
    '''
    Rates of counts k and n, shaped (entries, samples), per group (in order of
    first appearance): the rate of the counts summed over the group, and the
    mean and standard deviation over the group's replicas of each replica's
    rate.  Returns (group names, summed, mean, std), arrays (entries, groups).
    '''
    k, n = np.asarray(k, dtype=float), np.asarray(n, dtype=float)
    names = list(dict.fromkeys(groups))
    summed, mean, std = [], [], []
    for g in names:
        sel = [i for i, gi in enumerate(groups) if gi == g]
        summed.append(ratio(k[:, sel].sum(-1), n[:, sel].sum(-1)))
        reps = list(dict.fromkeys(replicas[i] for i in sel))
        rates = np.stack([ratio(k[:, [i for i in sel if replicas[i] == r]].sum(-1),
                                n[:, [i for i in sel if replicas[i] == r]].sum(-1))
                          for r in reps], axis=-1)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', RuntimeWarning)   # all-nan groups
            mean.append(np.nanmean(rates, axis=-1))
            std.append(np.nanstd(rates, axis=-1))
    return names, np.stack(summed, -1), np.stack(mean, -1), np.stack(std, -1)

def ratio(num, den):
    num, den = np.asarray(num, dtype=float), np.asarray(den, dtype=float)
    return np.divide(num, den, out=np.full(np.broadcast(num, den).shape, np.nan), where=den > 0)

def pad_channels(x, width):
    'Zero-pad the channel (second to last) dimension of x up to width.'
    return torch.nn.functional.pad(x, (0, 0, 0, width - x.shape[-2]))

def results_by_plane(prefix):
    'Input function: the run_one results for every angle of the plane.'
    def inner(w):
        return [f'<results>/{prefix}test_line_{w.plane}plane_{a}.pt' for a in angle_tags(w.plane)]
    return inner

def load_y_labels(filename, device='cpu'):
    t = torch.load(filename, map_location=device)
    return t['y'][0].detach().to(device), t['labels'][0].detach().to(device)

def pixel_eff_pur(y, labels, chan_range, threshold=0.5):
    lo, hi = chan_range
    y, labels = y[lo:hi], labels[lo:hi]
    n_matched = (y[torch.where(labels > 0)] > threshold).sum()
    return n_matched / labels.sum(), n_matched / (y > threshold).sum()

def plot_two_panel(top, bottom, xlabels, ylabels, ylim, ystep, fname):
    xpos = range(len(xlabels))
    plt.figure()
    for i, (vals, ylabel) in enumerate(zip((top, bottom), ylabels)):
        plt.subplot(2, 1, i + 1)
        plt.plot(vals)
        plt.xticks(xpos, xlabels)
        plt.ylabel(ylabel)
        plt.grid(True, axis='both')
        plt.ylim(*ylim)
        plt.yticks(np.arange(*ylim, ystep))
    plt.savefig(fname)
    plt.close()


rule help:
    run:
        print(USAGE)


#
# Model inference
#
# Line sample through the model, using only the sample's own view.
rule run_one:
    input:
        rec=lambda w: sample_h5(line_sample(w.plane, w.t1, w.t2), 'rec', [PLANE_TO_VIEW[w.plane]]),
        tru=lambda w: sample_h5(line_sample(w.plane, w.t1, w.t2), 'tru', [PLANE_TO_VIEW[w.plane]]),
    output:
        '<results>/test_line_{plane}plane_t1-{t1}-t2-{t2}.pt'
    params:
        device=TORCH_DEVICE,
        model=config['model'],
        app=config['app'],
        cfg=config['cfg'],
        entry=0,
    resources:
        gpu=USE_GPU
    shell: """
    wcpy dnn run_one \
    {input.rec} {input.tru} \
    --app {params.app} \
    -l {params.model} \
    -c {params.cfg} \
    -n {params.entry} \
    --manual-sigmoid \
    -d {params.device} \
    -o {output}
    """

use rule run_one as run_one_epoch with:
    pathvars:
        results="epoch_{epoch}",
    params:
        model=f'{config["cpt_dir"]}/checkpoint_{config["proc"]}_{config["cluster"]}_{{epoch}}.pt'

# Line sample through the model, feeding all three views.  (Cosmics samples
# are run per directory by run_n_xvu_cosmics.)
use rule run_one as run_one_xvu with:
    input:
        rec=lambda w: sample_h5(w.sample, 'rec', range(3)),
        tru=lambda w: sample_h5(w.sample, 'tru', range(3)),
    output:
        '<results>/{style}-test_{sample}.pt'
    params:
        model=lambda w: XVU_MODELS[w.style]['model'],
        cfg=lambda w: XVU_MODELS[w.style]['cfg'],
        app=lambda w: XVU_MODELS[w.style]['app'],
    wildcard_constraints:
        sample=LINE_SAMPLE,
        style=XVU_STYLES,

# The 'xvu' model at each checkpoint in cpt_dir.
use rule run_one_xvu as run_one_xvu_epoch with:
    pathvars:
        results="epoch_{epoch}",
    params:
        model=f'{config["cpt_dir"]}/checkpoint_{config["proc"]}_{config["cluster"]}_{{epoch}}.pt'
    wildcard_constraints:
        sample=LINE_SAMPLE,
        style='xvu',

# One view of a line sample through that view's own model (view_models).
use rule run_one as run_view with:
    input:
        rec=lambda w: sample_h5(w.sample, 'rec', [int(w.view)]),
        tru=lambda w: sample_h5(w.sample, 'tru', [int(w.view)]),
    output:
        os.path.join(config['dnnroi_dir'], 'view{view}-test_{sample}.pt')
    params:
        model=lambda w: config['view_models'][int(w.view)]['model'],
        cfg=lambda w: config['view_models'][int(w.view)]['cfg'],
        app=config['view_app'],
    wildcard_constraints:
        sample=LINE_SAMPLE,

# All events of a cosmics directory through the model, feeding all three views.
rule run_n_xvu_cosmics:
    input:
        rec=lambda w: sample_h5(f'cosmics_{w.dir}_ev0', 'rec', range(3)),
        tru=lambda w: sample_h5(f'cosmics_{w.dir}_ev0', 'tru', range(3)),
    output:
        [f'<results>/{{style}}-test_cosmics_{{dir}}_ev{e}.pt' for e in COSMICS_EVENTS]
    params:
        model=lambda w: XVU_MODELS[w.style]['model'],
        cfg=lambda w: XVU_MODELS[w.style]['cfg'],
        app=lambda w: XVU_MODELS[w.style]['app'],
    wildcard_constraints:
        style=XVU_STYLES,
    resources:
        gpu=USE_GPU
    run:
        run_n_events(input.rec, input.tru, params.app, params.model, params.cfg, output)

# All events of a cosmics directory, one view, through that view's own model.
rule run_n_view_cosmics:
    input:
        rec=lambda w: sample_h5(f'cosmics_{w.dir}_ev0', 'rec', [int(w.view)]),
        tru=lambda w: sample_h5(f'cosmics_{w.dir}_ev0', 'tru', [int(w.view)]),
    output:
        [os.path.join(config['dnnroi_dir'], f'view{{view}}-test_cosmics_{{dir}}_ev{e}.pt')
         for e in COSMICS_EVENTS]
    params:
        model=lambda w: config['view_models'][int(w.view)]['model'],
        cfg=lambda w: config['view_models'][int(w.view)]['cfg'],
        app=config['view_app'],
    resources:
        gpu=USE_GPU
    run:
        run_n_events(input.rec, input.tru, params.app, params.model, params.cfg, output)

# Every line sample (ALL_LINE_SAMPLES, all planes, angles and seeds) through a
# model in one wcpy process, as for cosmics.  See run_n_lines.
rule run_n_xvu_lines:
    input:
        rec=[f for s in ALL_LINE_SAMPLES for f in sample_h5(s, 'rec', range(3))],
        tru=[f for s in ALL_LINE_SAMPLES for f in sample_h5(s, 'tru', range(3))],
    output:
        [f'<results>/{{style}}-test_{s}.pt' for s in ALL_LINE_SAMPLES]
    params:
        model=lambda w: XVU_MODELS[w.style]['model'],
        cfg=lambda w: XVU_MODELS[w.style]['cfg'],
        app=lambda w: XVU_MODELS[w.style]['app'],
    resources:
        gpu=USE_GPU
    wildcard_constraints:
        style=XVU_STYLES,
    run:
        run_n_lines(input.rec, input.tru, params.app, params.model, params.cfg, {
            'rec_file_res': str([line_file_re('rec', v) for v in range(3)]),
            'tru_file_res': str([line_file_re('tru', v) for v in range(3)]),
        }, ALL_LINE_SAMPLES, output)

rule run_n_view_lines:
    input:
        rec=lambda w: [f for s in ALL_LINE_SAMPLES for f in sample_h5(s, 'rec', [int(w.view)])],
        tru=lambda w: [f for s in ALL_LINE_SAMPLES for f in sample_h5(s, 'tru', [int(w.view)])],
    output:
        [os.path.join(config['dnnroi_dir'], f'view{{view}}-test_{s}.pt') for s in ALL_LINE_SAMPLES]
    params:
        model=lambda w: config['view_models'][int(w.view)]['model'],
        cfg=lambda w: config['view_models'][int(w.view)]['cfg'],
        app=config['view_app'],
    resources:
        gpu=USE_GPU
    run:
        view = int(wildcards.view)
        run_n_lines(input.rec, input.tru, params.app, params.model, params.cfg, {
            'rec_file_re': line_file_re('rec', view),
            'tru_file_re': line_file_re('tru', view),
        }, ALL_LINE_SAMPLES, output)

# The batched rules make the same files as the per-sample ones.
ruleorder: run_n_xvu_lines > run_one_xvu
ruleorder: run_n_view_lines > run_view

# Concatenate the per-view outputs into the xvu channel layout.  A view model
# that covers fewer channels than view_widths (e.g. one face of w) is
# zero-padded at the end.
rule merge_views:
    input:
        [os.path.join(config['dnnroi_dir'], f'view{v}-test_{{sample}}.pt') for v in range(3)]
    output:
        os.path.join(config['dnnroi_dir'], 'merged-test_{sample}.pt')
    run:
        parts = [torch.load(f, map_location='cpu') for f in input]
        merged = {
            key: torch.cat([pad_channels(p[key], width)
                            for p, width in zip(parts, config['view_widths'])], dim=-2)
            for key in parts[0]
        }
        torch.save(merged, output[0])

# Entry {entry} of the samples in config paths.
rule run_n:
    output:
        "run_one_{entry}.pt"
    params:
        paths=config['paths'],
        model=config['model'],
        device=TORCH_DEVICE,
        cfg=config['cfg'],
        app=config['app'],
    resources:
        gpu=USE_GPU
    shell: """
    wcpy dnn run_one -l {params.model} -c {params.cfg} --manual-sigmoid -d {params.device} \
        -a {params.app} -o {output} -n {wildcards.entry} {params.paths}
    """


#
# Pixel efficiency/purity vs angle
#
rule plot_eff_pur:
    input:
        results_by_plane('')
    output:
        png='<results>/all_eff_pur_{plane}plane.png',
        npz='<results>/all_eff_pur_{plane}plane.npz',
    params:
        test_tpcs=[int(i) for i in config['test_tpcs']], # views, not TPCs
        threshold=0.5,
    run:
        all_effs = []
        all_purs = []
        for f in input:
            y, labels = load_y_labels(f)
            for view in params.test_tpcs:
                eff, pur = pixel_eff_pur(y, labels, CHAN_RANGES[PLANES[view]], params.threshold)
                all_effs.append(eff)
                all_purs.append(pur)
        print(all_effs)
        print(all_purs)

        plot_two_panel(all_effs, all_purs, angle_labels(wildcards.plane),
                       ('Efficiency', 'Purity'), (0., 1.5), .2, output.png)
        np.savez(output.npz, effs=all_effs, purs=all_purs)

use rule plot_eff_pur as plot_eff_pur_epoch with:
    pathvars:
        results="epoch_{epoch}",

use rule plot_eff_pur as plot_eff_pur_xvu with:
    input:
        results_by_plane('xvu-')
    output:
        png='<results>/all_eff_pur_{plane}plane-xvu.png',
        npz='<results>/all_eff_pur_{plane}plane-xvu.npz',
    params:
        test_tpcs=lambda w: [PLANE_TO_VIEW[w.plane]]

use rule plot_eff_pur_xvu as plot_eff_pur_xvu_epoch with:
    pathvars:
        results="epoch_{epoch}",


#
# Event displays: feature, prediction and labels for one line sample
#
rule plot_feat_y_labels:
    input:
        '<results>/test_line_{plane}plane_t1-{t1}-t2-{t2}.pt'
    output:
        pickle=temp('line-display-{plane}plane-t1_{t1}-t2_{t2}.pickle'),
        png='line-display-{plane}plane-t1_{t1}-t2_{t2}.png'
    run:
        import pickle
        t = torch.load(input[0], map_location='cpu')
        fig, axs = plt.subplots(3, 1)
        for ax, f in zip(axs, [t['feat'], t['y'] > .5, t['labels']]):
            im = ax.imshow(f[0].detach().cpu(), aspect='auto')
            fig.colorbar(im, ax=ax)
        fig.savefig(output.png)
        with open(output.pickle, 'wb') as out:
            pickle.dump(fig, out)
        plt.close(fig)

use rule plot_feat_y_labels as plot_feat_y_labels_xvu with:
    input:
        '<results>/xvu-test_line_{plane}plane_t1-{t1}-t2-{t2}.pt'
    output:
        pickle=temp('xvu-line-display-{plane}plane-t1_{t1}-t2_{t2}.pickle'),
        png='xvu-line-display-{plane}plane-t1_{t1}-t2_{t2}.png'

def display_pickles(prefix):
    def inner(w):
        return [f'{prefix}line-display-{w.plane}plane-{depo_tag(t1, t2)}.pickle'
                for t1, t2 in ANGLES[w.plane]]
    return inner

rule collect_feat_y_labels_uvplane:
    input:
        display_pickles('')
    output:
        'line-displays_{plane}plane.tar'
    shell:
        'tar -cf {output} {input}'

use rule collect_feat_y_labels_uvplane as xvu_collect_feat_y_labels_uvplane with:
    input:
        display_pickles('xvu-')
    output:
        'xvu-line-displays_{plane}plane.tar'


#
# Threshold-dependent ROI and pixel metrics per line sample
#
rule angled_roi_eff_pur:
    input:
        "<results>/xvu-test_line_{plane}plane_{angles}.pt"
    output:
        temp("<results>/threshold-{threshold}-xvu-roi_effs_purs_{plane}plane_{angles}.npz")
    run:
        y, labels = load_y_labels(input[0])
        lo, hi = CHAN_RANGES[wildcards.plane]
        res = roi_metrics.roi_metrics(y[lo:hi], labels[lo:hi], threshold=float(wildcards.threshold))
        eff, pur = res['efficiency'].item(), res['purity'].item()
        print(eff, pur)
        np.savez(output[0], effs=np.array([eff]), purs=np.array([pur]))

rule pixel_eff_pur_xvu_line:
    input:
        '<results>/xvu-test_line_{plane}plane_{angles}.pt'
    output:
        '<results>/threshold-{threshold}-xvu-pixel_effs_purs_{plane}plane_{angles}.npz'
    run:
        y, labels = load_y_labels(input[0])
        eff, pur = pixel_eff_pur(y, labels, CHAN_RANGES[wildcards.plane], float(wildcards.threshold))
        np.savez(output[0], effs=np.array([eff]), purs=np.array([pur]))

rule merge_roi_eff_pur:
    input:
        lambda w: [f'<results>/threshold-{w.threshold}-xvu-{w.type}_effs_purs_{w.plane}plane_{a}.npz'
                   for a in angle_tags(w.plane)]
    output:
        "<results>/threshold-{threshold}-xvu-merged_{type}_effs_purs_{plane}plane_high_end.npz"
    run:
        effs, purs = [], []
        for f in input:
            t = np.load(f)
            effs.append(t['effs'])
            purs.append(t['purs'])
        np.savez(output[0], effs=np.array(effs), purs=np.array(purs))

rule plot_merged_roi_eff_pur:
    input:
        "<results>/threshold-{threshold}-xvu-merged_{type}_effs_purs_{plane}plane_high_end.npz"
    output:
        '<results>/threshold-{threshold}-xvu-all_{type}_eff_pur_{plane}plane.png'
    params:
        ylim=(-.4, .4)
    run:
        t = np.load(input[0])
        plot_two_panel(1. - t['effs'], 1. - t['purs'], angle_labels(wildcards.plane),
                       ('1-Efficiency', '1-Purity'), params.ylim, .1, output[0])


#
# ROI tables
#
rule xvu_angled_roi_table:
    input:
        "<results>/xvu-test_line_{plane}plane_{angles}.pt"
    output:
        "<results>/xvu-{type}_roi_table_{plane}plane_{angles}-thresh-{threshold}.pt",
    run:
        lo, hi = CHAN_RANGES[wildcards.plane]
        y, labels = load_y_labels(input[0])
        if wildcards.type == 'reco':
            mask = y[lo:hi] > float(wildcards.threshold)
        else:
            mask = labels[lo:hi] == 1
        torch.save(roi_metrics.roi_table(mask), output[0])

rule all_roi_tables:
    input:
        [f'<results>/xvu-{type}_roi_table_{p}plane_{a}-thresh-0.5.pt'
         for p in PLANES for type in ('true', 'reco') for a in angle_tags(p)]

rule roi_lengths:
    input:
        lambda w: [f'<results>/xvu-{w.type}_roi_table_{w.plane}plane_{a}-thresh-0.5.pt'
                   for a in angle_tags(w.plane)]
    output:
        "<results>/xvu-{type}_roi_table_info_{plane}plane-thresh-0.5.npz"
    run:
        nrois = []
        lengths = []
        for f in input:
            t = torch.load(f)
            lengths += t['length'].tolist()
            nrois.append(len(t['length']))
        np.savez(output[0], lengths=lengths, nrois=nrois, names=list(input))


#
# Threshold scans and F-beta
#
def roi_scan(y, labels, thresholds):
    return roi_metrics.threshold_scan(y, labels, thresholds, as_eff_pur=False)

def pixel_scan(y, labels, thresholds):
    thresholds = torch.as_tensor(thresholds, dtype=y.dtype, device=y.device)
    n_true = labels.sum()
    reco = (y.unsqueeze(-1).expand(-1, -1, len(thresholds)) > thresholds)
    n_reco = reco.sum((0, 1))
    matched = ((labels > 0).unsqueeze(-1).expand(-1, -1, reco.shape[-1]) == reco) * (reco > 0)

    n_reco_matched = matched.sum((0, 1))
    n_true_matched = n_reco_matched # pixels match one-to-one
    return {
        "threshold": thresholds,
        "n_reco": n_reco,
        "n_reco_matched": n_reco_matched,
        "n_true_matched": n_true_matched,
        "n_true": n_true,
        "efficiency": (n_true_matched / n_true),
        "purity": (n_reco_matched / n_reco),
    }

SCANS = {'roi': roi_scan, 'pixel': pixel_scan}

def aggregate_scans(results, output):
    'Save the per-sample scans in results and their sums over samples.'
    def as_np(x):
        return x.cpu().numpy() if isinstance(x, torch.Tensor) else x

    cols = {k: [] for k in ('efficiency', 'purity', 'n_true', 'n_reco',
                            'n_true_matched', 'n_reco_matched', 'threshold')}
    for t in results:
        for k in cols:
            cols[k].append(as_np(t[k]))

    np.savez(
        output,
        effs=cols['efficiency'],
        purs=cols['purity'],
        n_trues=cols['n_true'],
        n_recos=cols['n_reco'],
        n_recos_matched=cols['n_reco_matched'],
        n_trues_matched=cols['n_true_matched'],
        n_recos_summed=np.array(cols['n_reco']).sum(axis=0),
        n_trues_summed=sum(cols['n_true']),
        n_recos_matched_summed=np.array(cols['n_reco_matched']).sum(axis=0),
        n_trues_matched_summed=np.array(cols['n_true_matched']).sum(axis=0),
        thresholds=cols['threshold'][0], # same for every entry
    )

# One run_n entry (see run_n and config paths).
rule threshold_scan:
    input:
        "run_one_{entry}.pt"
    output:
        "threshold_scan_{type}_plane_{plane}_run_one_{entry}.pt"
    wildcard_constraints:
        type='roi|pixel',
    params:
        start=config['threshold_start'],
        step=config['threshold_step'],
    run:
        y, labels = load_y_labels(input[0], TORCH_DEVICE)
        lo, hi = CHAN_RANGES[wildcards.plane]
        thresholds = np.arange(params.start, 1., params.step)
        torch.save(SCANS[wildcards.type](y[lo:hi], labels[lo:hi], thresholds), output[0])

rule aggregate_scan:
    input:
        collect('threshold_scan_{{type}}_plane_{{plane}}_run_one_{entry}.pt',
                entry=range(config['nentries']))
    output:
        'aggregated_scan_{type}_plane_{plane}.npz'
    run:
        aggregate_scans([torch.load(f, map_location='cpu') for f in input], output[0])

rule compute_fbeta:
    input:
        'aggregated_scan_{type}_plane_{plane}.npz'
    output:
        'fbeta_aggregated_scan_{type}_plane_{plane}_{beta}.npz'
    run:
        t = np.load(input[0])
        R = t['n_trues_matched_summed'] / t['n_trues_summed']
        P = t['n_recos_matched_summed'] / t['n_recos_summed']
        beta = float(wildcards.beta)
        fbeta = (1 + beta * beta) * (R * P) / (beta * beta * P + R)
        np.savez(
            output[0],
            thresholds=t['thresholds'],
            beta=beta,
            fbeta=fbeta,
            maxloc=t['thresholds'][np.argmax(fbeta)],
            maxval=np.max(fbeta),
        )

rule plot_fbeta_scan:
    input:
        expand('fbeta_aggregated_scan_{{type}}_plane_{{plane}}_{beta}.npz', beta=FBETA_BETAS)
    output:
        png='aggregated_scan_{type}_plane_{plane}_fbetas.png',
    run:
        plt.figure()
        for f in input:
            t = np.load(f)
            maxloc = t['maxloc']
            plt.plot(t['thresholds'], t['fbeta'],
                     label=r"$\beta$" + f"={t['beta']:.2f} | max = {maxloc:.2f}")
            plt.scatter(maxloc, t['maxval'])
        plt.xlabel('Threshold')
        plt.ylabel(r'$F_\beta$')
        plt.legend()
        plt.savefig(output.png)
        plt.close()

# The same scans for any style's output of a {sample}, over every SCAN_PLANES
# range at once, aggregated over a {set} of samples.
rule scan_sample:
    input:
        lambda w: f'{STYLE_PREFIX[w.style]}test_{w.sample}.pt'
    output:
        '<results>/scan-{style}-test_{sample}.pt'
    params:
        start=config['threshold_start'],
        step=config['threshold_step'],
    run:
        y, labels = load_y_labels(input[0], TORCH_DEVICE)
        thresholds = np.arange(params.start, 1., params.step)
        res = {}
        for plane in SCAN_PLANES:
            lo, hi = CHAN_RANGES[plane]
            for kind, scan in SCANS.items():
                res[f'{kind}_{plane}'] = scan(y[lo:hi], labels[lo:hi], thresholds)
        torch.save(res, output[0])

rule aggregate_set_scan:
    input:
        lambda w: [f'<results>/scan-{w.style}-test_{s}.pt' for s in set_samples(w.set)]
    output:
        '<results>/aggregated_scan_{type}_plane_{plane}_{style}_{set}.npz'
    wildcard_constraints:
        type='roi|pixel',
    run:
        key = f'{wildcards.type}_{wildcards.plane}'
        aggregate_scans([torch.load(f, map_location='cpu')[key] for f in input], output[0])

use rule compute_fbeta as compute_fbeta_set with:
    input:
        '<results>/aggregated_scan_{type}_plane_{plane}_{style}_{set}.npz'
    output:
        '<results>/fbeta_aggregated_scan_{type}_plane_{plane}_{style}_{set}_{beta}.npz'

use rule plot_fbeta_scan as plot_fbeta_scan_set with:
    input:
        expand('<results>/fbeta_aggregated_scan_{{type}}_plane_{{plane}}_{{style}}_{{set}}_{beta}.npz',
               beta=FBETA_BETAS)
    output:
        png='<results>/aggregated_scan_{type}_plane_{plane}_{style}_{set}_fbetas.png',

# Operating point of a style: the thresholds maximising {type} F1 on the
# OP_TUNE_SET samples, one per view (see OP_PLANES).
rule op_thresholds:
    input:
        unpack(lambda w: {p: f'<results>/fbeta_aggregated_scan_{w.type}_plane_{p}_{w.style}_{OP_TUNE_SET}_1.0.npz'
                          for p in set(OP_PLANES[w.mode])})
    output:
        '<results>/op-{type}F1_{mode}-{style}.json'
    run:
        planes = OP_PLANES[wildcards.mode]
        best = {p: np.load(input[p]) for p in set(planes)}
        with open(output[0], 'w') as out:
            json.dump({
                'tune_set': OP_TUNE_SET,
                'planes': planes,
                'thresholds': [float(best[p]['maxloc']) for p in planes],
                'f1': [float(best[p]['maxval']) for p in planes],
            }, out, indent=1)


#
# Cross-plane agreement on true (u, v, w, tick) trios
#
def trio_rates(run_one, trio_file, group, threshold, output):
    '''
    Fraction of the true trios in group of trio_file that the model output
    run_one predicts in each TRIO_COMBOS pattern of views.  threshold is one
    value or one per view.
    '''
    import h5py
    with h5py.File(trio_file, 'r') as ftrio:
        trios = ftrio[group]['trio_uvwt'][...].astype(np.int64)
    y, labels = load_y_labels(run_one)
    y, labels = y.numpy(), labels.numpy()
    thresholds = np.broadcast_to(np.asarray(threshold, dtype=float), (3,))

    # Keep trios whose pixel is true in all three views.
    tick = trios[:, -1]
    true_uvw = np.stack([labels[trios[:, i], tick] for i in range(3)], axis=-1)
    valid_trios = trios[np.all(true_uvw > 0., axis=1)]

    tick = valid_trios[:, -1]
    ys = np.stack([y[valid_trios[:, i], tick] > thresholds[i] for i in range(3)], axis=-1)
    n_valid = len(ys)

    rates = {}
    for name, pattern in TRIO_COMBOS.items():
        rates[name] = np.sum(np.all(ys == np.array(pattern), axis=-1)) * 1. / n_valid
        print(f'{name} -- {rates[name]:.3f}')

    np.savez(output, ys=ys, valid_trios=valid_trios, n_valid=n_valid,
             thresholds=thresholds, **rates)

rule sample_trios_cross_plane:
    input:
        run_one=lambda w: f'{STYLE_PREFIX[w.style]}test_{w.sample}.pt',
        trios=lambda w: sample_h5_base(w.sample) + 'trio.h5',
    output:
        "<results>/threshold-{threshold}-{style}-trio-crossplane-{sample}.npz"
    run:
        trio_rates(input.run_one, input.trios, trio_group(wildcards.sample),
                   float(wildcards.threshold), output[0])

# The same at a style's operating point (op_thresholds).
rule op_trios_cross_plane:
    input:
        run_one=lambda w: f'{STYLE_PREFIX[w.style]}test_{w.sample}.pt',
        trios=lambda w: sample_h5_base(w.sample) + 'trio.h5',
        op='<results>/op-{type}F1_{mode}-{style}.json',
    output:
        "<results>/op-{type}F1_{mode}-{style}-trio-crossplane-{sample}.npz"
    run:
        with open(input.op) as f:
            thresholds = json.load(f)['thresholds']
        trio_rates(input.run_one, input.trios, trio_group(wildcards.sample),
                   thresholds, output[0])

# The same for a line sample, under its older name.
rule trios_cross_plane:
    input:
        run_one=lambda w: f'{STYLE_PREFIX[w.style]}test_line_{w.plane}plane_{w.angles}.pt',
        trios=lambda w: sample_h5_base(f'line_{w.plane}plane_{w.angles}') + 'trio.h5',
    output:
        "threshold-{threshold}-{style}-trio-crossplane-{plane}plane-{angles}-g4-trio.npz"
    run:
        trio_rates(input.run_one, input.trios, '0', float(wildcards.threshold), output[0])

rule merge_trio_cross_plane:
    input:
        lambda w: [f"<results>/threshold-{w.threshold}-{w.style}-trio-crossplane-{s}.npz"
                   for s in set_samples(w.set)]
    output:
        "<results>/threshold-{threshold}-{style}-trio-merged-crossplane_{set}.npz"
    run:
        merged = {name: [] for name in list(TRIO_COMBOS) + ['n_valid']}
        for f in input:
            t = np.load(f)
            for name in merged:
                merged[name].append(t[name])
        np.savez(output[0], samples=set_samples(wildcards.set),
                 **{name: np.array(v) for name, v in merged.items()})

use rule merge_trio_cross_plane as merge_op_trio_cross_plane with:
    input:
        lambda w: [f"<results>/op-{w.type}F1_{w.mode}-{w.style}-trio-crossplane-{s}.npz"
                   for s in set_samples(w.set)]
    output:
        "<results>/op-{type}F1_{mode}-{style}-trio-merged-crossplane_{set}.npz"

def set_labels(name):
    return None if name.startswith('cosmics') else angle_labels(name[0])

rule plot_trio_rates:
    input:
        "<results>/threshold-{threshold}-{style}-trio-merged-crossplane_{set}.npz"
    output:
        "<results>/threshold-{threshold}-{style}-trio-merged-crossplane_{set}.png"
    run:
        t = np.load(input[0])
        xpos = np.arange(len(t['uvw']))
        fig, ax = plt.subplots()
        bottom = np.zeros(len(xpos))
        for name in TRIO_COMBOS:
            ax.bar(xpos, t[name], label=name, bottom=bottom)
            bottom += t[name]
        labels = set_labels(wildcards.set)
        if labels:
            ax.set_xticks(xpos, labels)
        else:
            ax.set_xlabel('Sample')
        ax.legend()
        fig.savefig(output[0])
        plt.close(fig)

def replica_kind(set_name):
    return 'seeds' if re.fullmatch(LINE_SET, set_name) else 'events'

def grouped_bars(ax, groups, mean, std, labels, colors, zoom=True):
    '''
    Bars of mean, shaped (entries, groups), side by side per group, with
    error bars of std.  zoom: start the y axis just below the data.
    '''
    n = len(labels)
    width = 0.8 / n
    xpos = np.arange(len(groups))
    for i in range(n):
        x = xpos + (i - (n - 1) / 2) * width
        ax.bar(x, mean[i], width, yerr=std[i], color=colors[i], label=labels[i],
               error_kw=dict(ecolor='black', elinewidth=0.8, capsize=1.5))
    ax.set_xticks(xpos, groups)
    ax.grid(True, axis='y')
    if zoom:
        low = np.nanmin(mean - std)
        high = np.nanmax(mean + std)
        ax.set_ylim(max(0., low - 0.02), max(1.0, high) + 0.005)

def plot_trio_compare(files, labels, title, set_name, png, npz):
    '''
    Trio rates of several merged-crossplane files side by side.  Bars are the
    mean and error bars the standard deviation over replicas (seeds for line
    samples, events for cosmics, see sample_group).  Panels: U & V & W per entry, the other
    combinations, and U & V & W per angle (line samples) or per event.
    '''
    files = [np.load(f) for f in files]
    names = list(TRIO_COMBOS)
    samples = [str(s) for s in files[0]['samples']]
    groups, replicas = (list(x) for x in zip(*(sample_group(s) for s in samples)))
    n_valid = np.array([t['n_valid'] for t in files])        # (entries, samples)
    # Trios in each combination; a sample with no valid trios has nan rates.
    counts = {name: np.nan_to_num(np.array([t[name] for t in files]) * n_valid) for name in names}

    overall = {name: summed_and_spread(counts[name], n_valid, ['all'] * len(samples), replicas)
               for name in names}
    pooled, mean, std = (np.stack([overall[name][j][:, 0] for name in names], axis=-1)
                         for j in (1, 2, 3))                   # (entries, combos)
    group_names, uvw_sum, uvw_mean, uvw_std = summed_and_spread(counts['uvw'], n_valid, groups, replicas)
    for label, row, err in zip(labels, pooled, std):
        print(label, ' '.join(f'{n}={r:.3f}+-{e:.3f}' for n, r, e in zip(names, row, err)))

    n = len(labels)
    colors = plt.cm.tab20(np.arange(n) % 20)
    fig, (top, middle, bottom) = plt.subplots(
        3, 1, figsize=(13, 4 + 0.35 * n + 8), height_ratios=[0.35 * n + 1, 4, 4])
    spread = replica_kind(set_name)
    fig.suptitle(f'{title}\nmean $\\pm$ std over {spread}')

    # U & V & W, the headline number: one labelled bar per entry, zoomed in.
    iuvw = names.index('uvw')
    ypos = np.arange(n)[::-1]
    top.barh(ypos, mean[:, iuvw], xerr=std[:, iuvw], color=colors,
             error_kw=dict(ecolor='black', capsize=2))
    for y, value, err in zip(ypos, mean[:, iuvw], std[:, iuvw]):
        top.text(min(value + err, 1.0), y, f' {value:.3f} $\\pm$ {err:.3f}', va='center', fontsize=8)
    top.set_yticks(ypos, labels, fontsize=8)
    low = np.nanmin(mean[:, iuvw] - std[:, iuvw])
    top.set_xlim(max(0., low - 0.03), 1.02)
    top.set_xlabel('U & V & W fraction of true trios')
    top.grid(True, axis='x')

    # The rest: trios seen in two views or only one.
    others = [i for i, name in enumerate(names) if name != 'uvw']
    grouped_bars(middle, [names[i] for i in others], mean[:, others], std[:, others],
                 labels, colors, zoom=False)
    middle.set_ylabel('Fraction of true trios')
    middle.legend(fontsize=7, loc='upper left', bbox_to_anchor=(1.01, 1.0))

    if spread == 'seeds':
        grouped_bars(bottom, group_names, uvw_mean, uvw_std, labels, colors)
        bottom.set_xlabel(f'{set_name[0]}-plane line angles (t1, t2), '
                          f'{len(set(replicas))} seed(s) each')
    else:
        for i, t in enumerate(files):
            bottom.plot(t['uvw'], marker='o', color=colors[i])
        bottom.set_xlabel('Sample')
        bottom.grid(True, axis='both')
    bottom.set_ylabel('U & V & W fraction')

    fig.tight_layout()
    fig.savefig(png)
    plt.close(fig)

    np.savez(npz, labels=labels, combos=names, pooled=pooled, mean=mean, std=std,
             samples=samples, uvw=np.array([t['uvw'] for t in files]), n_valid=n_valid,
             groups=group_names, uvw_sum=uvw_sum, uvw_mean=uvw_mean, uvw_std=uvw_std)

# compare_styles side by side at one fixed threshold.
rule compare_trio_rates:
    input:
        lambda w: [f"<results>/threshold-{w.threshold}-{style}-trio-merged-crossplane_{w.set}.npz"
                   for style in config['compare_styles']]
    output:
        png="<results>/threshold-{threshold}-trio-compare_{set}.png",
        npz="<results>/threshold-{threshold}-trio-compare_{set}.npz",
    run:
        plot_trio_compare(input, config['compare_styles'],
                          f'{wildcards.set}, threshold {wildcards.threshold}',
                          wildcards.set, output.png, output.npz)

# op_compare (style, mode) pairs side by side, each at its own operating point.
rule compare_op_trio_rates:
    input:
        rates=lambda w: [f"<results>/op-{w.type}F1_{mode}-{style}-trio-merged-crossplane_{w.set}.npz"
                         for style, mode in config['op_compare']],
        ops=lambda w: [f'<results>/op-{w.type}F1_{mode}-{style}.json'
                       for style, mode in config['op_compare']],
    output:
        png="<results>/op-{type}F1-trio-compare_{set}.png",
        npz="<results>/op-{type}F1-trio-compare_{set}.npz",
    run:
        labels = []
        for (style, mode), f in zip(config['op_compare'], input.ops):
            with open(f) as fp:
                thresholds = json.load(fp)['thresholds']
            labels.append(f'{style} {mode} ({"/".join(f"{t:.2f}" for t in thresholds)})')
        plot_trio_compare(input.rates, labels,
                          f'{wildcards.set}, best {wildcards.type} F1 on {OP_TUNE_SET}',
                          wildcards.set, output.png, output.npz)

OP_KINDS = ('roi', 'pixel')

#
# Efficiency and purity per line sample at the operating points
#
# Channel range scored for a line sample of view 0, 1, 2 (its own plane; w0 is
# the w face with simulated signal).
LINE_METRIC_RANGES = ['u', 'v', 'w0']
LINE_SET = r'[uvw]plane_high_end'

# Counts behind efficiency (n_true_matched / n_true) and purity
# (n_reco_matched / n_reco).
COUNT_KEYS = ('n_true', 'n_true_matched', 'n_reco', 'n_reco_matched')

def line_counts(kind, y, labels, threshold):
    'COUNT_KEYS counts of ROIs or pixels in (channel, tick) images.'
    if kind == 'roi':
        res = roi_metrics.roi_metrics(y, labels, threshold=threshold)
        counts = [res[key] for key in COUNT_KEYS]
    else:
        truth = labels > 0
        reco = y > threshold
        n_matched = (truth & reco).sum()
        counts = [truth.sum(), n_matched, reco.sum(), n_matched]
    return {key: int(count) for key, count in zip(COUNT_KEYS, counts)}


# {type} efficiency and purity on the sample's own plane, at the style's best
# {type} F1 threshold for that view.
rule op_line_effpur:
    input:
        run_one=lambda w: f'{STYLE_PREFIX[w.style]}test_{w.sample}.pt',
        op='<results>/op-{type}F1_{mode}-{style}.json',
    output:
        '<results>/op-{type}F1_{mode}-{style}-effpur-{sample}.npz'
    wildcard_constraints:
        sample=LINE_SAMPLE,
        type='roi|pixel',
    run:
        view = PLANE_TO_VIEW[parse_sample(wildcards.sample)[1]['plane']]
        with open(input.op) as f:
            threshold = json.load(f)['thresholds'][view]
        lo, hi = CHAN_RANGES[LINE_METRIC_RANGES[view]]
        y, labels = load_y_labels(input.run_one)
        counts = line_counts(wildcards.type, y[lo:hi], labels[lo:hi], threshold)
        np.savez(output[0], threshold=threshold, chan_range=(lo, hi),
                 eff=ratio(counts['n_true_matched'], counts['n_true']),
                 pur=ratio(counts['n_reco_matched'], counts['n_reco']), **counts)

rule merge_op_line_effpur:
    input:
        lambda w: [f'<results>/op-{w.type}F1_{w.mode}-{w.style}-effpur-{s}.npz'
                   for s in set_samples(w.set)]
    output:
        '<results>/op-{type}F1_{mode}-{style}-effpur_{set}.npz'
    wildcard_constraints:
        set=LINE_SET,
        type='roi|pixel',
    run:
        parts = [np.load(f) for f in input]
        np.savez(output[0], samples=set_samples(wildcards.set),
                 effs=np.array([float(t['eff']) for t in parts]),
                 purs=np.array([float(t['pur']) for t in parts]),
                 thresholds=np.array([float(t['threshold']) for t in parts]),
                 **{key: np.array([int(t[key]) for t in parts]) for key in COUNT_KEYS})

# op_compare entries' efficiency and purity vs angle: mean and std over seeds.
rule compare_op_line_effpur:
    input:
        lambda w: [f'<results>/op-{w.type}F1_{mode}-{style}-effpur_{w.set}.npz'
                   for style, mode in config['op_compare']]
    output:
        png='<results>/op-{type}F1-effpur-compare_{set}.png',
        npz='<results>/op-{type}F1-effpur-compare_{set}.npz',
    wildcard_constraints:
        set=LINE_SET,
        type='roi|pixel',
    run:
        parts = [np.load(f) for f in input]
        labels = [f'{style} {mode} ({float(t["thresholds"][0]):.2f})'
                  for (style, mode), t in zip(config['op_compare'], parts)]
        counts = {key: np.array([t[key] for t in parts]) for key in COUNT_KEYS}   # (entries, samples)
        samples = [str(s) for s in parts[0]['samples']]
        groups, replicas = (list(x) for x in zip(*(sample_group(s) for s in samples)))
        nseeds = len(set(replicas))

        res = {}
        for name, num, den in (('eff', 'n_true_matched', 'n_true'), ('pur', 'n_reco_matched', 'n_reco')):
            angles, res[f'{name}_sum'], res[f'{name}_mean'], res[f'{name}_std'] = \
                summed_and_spread(counts[num], counts[den], groups, replicas)
            _, s, m, d = summed_and_spread(counts[num], counts[den], ['all'] * len(samples), replicas)
            res[f'{name}_all_sum'], res[f'{name}_all_mean'], res[f'{name}_all_std'] = s[:, 0], m[:, 0], d[:, 0]
        # Counts summed over seeds per angle, (entries, angles).
        sums = {key: np.stack([counts[key][:, [i for i, g in enumerate(groups) if g == a]].sum(-1)
                               for a in angles], axis=-1) for key in COUNT_KEYS}

        colors = plt.cm.tab20(np.arange(len(labels)) % 20)
        plane = wildcards.set[0]
        # The truth does not depend on the model: show its size under each angle.
        xlabels = [f'{a}\n{n} true' for a, n in zip(angles, sums['n_true'][0])]

        fig, axs = plt.subplots(2, 1, figsize=(15, 10), sharex=True)
        for ax, name in zip(axs, ('efficiency', 'purity')):
            key = name[:3]
            grouped_bars(ax, xlabels, res[f'{key}_mean'], res[f'{key}_std'], labels, colors)
            ax.set_ylabel(f'{wildcards.type} {name}')
        axs[0].legend(fontsize=7, loc='upper left', bbox_to_anchor=(1.01, 1.0))
        axs[1].set_xlabel(f'{plane}-plane line angles (t1, t2) and number of true {wildcards.type}s '
                          f'summed over {nseeds} seed(s)')
        fig.suptitle(f'{wildcards.set}: {wildcards.type} efficiency and purity on the '
                     f'{LINE_METRIC_RANGES[PLANE_TO_VIEW[plane]]} channels, '
                     f'at best {wildcards.type} F1 on {OP_TUNE_SET}\n'
                     f'mean $\\pm$ std over seeds')
        fig.tight_layout()
        fig.savefig(output.png)
        plt.close(fig)
        np.savez(output.npz, labels=labels, angles=angles, nseeds=nseeds, samples=samples,
                 **res, **{f'{key}_sum': value for key, value in sums.items()},
                 **{f'{key}_sample': value for key, value in counts.items()})

def op_styles():
    'The styles in op_compare, in order of first appearance.'
    return list(dict.fromkeys(style for style, mode in config['op_compare']))

def describe_set(name):
    if name == 'cosmics':
        return 'cosmics: ' + ', '.join(cosmics_dir_names())
    if name.startswith('cosmics_'):
        return f'{name}: ' + ', '.join(cosmics_split()[name])
    return name

def text_page(pdf, lines):
    fig = plt.figure(figsize=(11, 8.5))
    fig.text(0.04, 0.96, '\n'.join(lines), va='top', family='monospace', fontsize=8)
    pdf.savefig(fig, bbox_inches='tight')
    plt.close(fig)

def table_page(pdf, title, col_labels, row_labels, cells):
    fig, ax = plt.subplots(figsize=(11, 8.5))
    ax.axis('off')
    ax.set_title(title)
    table = ax.table(cellText=cells, colLabels=col_labels, rowLabels=row_labels, loc='center')
    table.auto_set_font_size(False)
    table.set_fontsize(7)
    table.auto_set_column_width(list(range(-1, len(col_labels))))
    table.scale(1, 1.3)
    # tight: wide tables extend past the axes; grow the page to fit them.
    pdf.savefig(fig, bbox_inches='tight')
    plt.close(fig)

def image_page(pdf, png):
    img = plt.imread(png)
    height, width = img.shape[:2]
    fig = plt.figure(figsize=(11, 11 * height / width))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.imshow(img)
    ax.axis('off')
    pdf.savefig(fig)
    plt.close(fig)

# Everything op_compare produces, in one PDF: the setup; a summary of the
# U & V & W rate per entry and report_sets set; per set, a table of operating
# points and pooled trio rates per F1 kind and the comparison plots; and the
# F-beta scans the operating points were picked from.
rule op_report:
    input:
        png=[f'<results>/op-{kind}F1-trio-compare_{s}.png'
             for s in config['report_sets'] for kind in OP_KINDS],
        npz=[f'<results>/op-{kind}F1-trio-compare_{s}.npz'
             for s in config['report_sets'] for kind in OP_KINDS],
        effpur=[f'<results>/op-{kind}F1-effpur-compare_{s}.{ext}'
                for s in config['report_sets'] if re.fullmatch(LINE_SET, s)
                for kind in OP_KINDS for ext in ('png', 'npz')],
        ops=[f'<results>/op-{kind}F1_{mode}-{style}.json'
             for kind in OP_KINDS for style, mode in config['op_compare']],
        fbeta=[f'<results>/fbeta_aggregated_scan_{kind}_plane_{p}_{style}_{OP_TUNE_SET}_{beta}.npz'
               for kind in OP_KINDS for style in op_styles() for p in SCAN_PLANES
               for beta in FBETA_BETAS],
    output:
        '<results>/op-trio-report.pdf'
    run:
        from matplotlib.backends.backend_pdf import PdfPages
        res = os.path.dirname(output[0])
        sets = config['report_sets']
        entries = [f'{style} {mode}' for style, mode in config['op_compare']]
        compare = {(s, kind): np.load(f'{res}/op-{kind}F1-trio-compare_{s}.npz')
                   for s in sets for kind in OP_KINDS}
        ops = {}
        for kind in OP_KINDS:
            for style, mode in config['op_compare']:
                with open(f'{res}/op-{kind}F1_{mode}-{style}.json') as f:
                    ops[kind, style, mode] = json.load(f)

        lines = ['Trio rates at best-F1 operating points', '',
                 f'Tuned on    {describe_set(OP_TUNE_SET)}']
        lines += [f'Evaluated   {describe_set(s)}' for s in sets]
        lines += [f'Thresholds  {config["threshold_start"]} to 1 in steps of {config["threshold_step"]}',
                  '', 'Models:']
        for style in op_styles():
            if style in XVU_MODELS:
                lines += [f'  {style}: {XVU_MODELS[style]["model"]}',
                          f'  {" " * len(style)}  {XVU_MODELS[style]["cfg"]}']
            elif style == 'merged':
                lines += [f'  merged view {v}: {m["model"]}' for v, m in enumerate(config['view_models'])]

        with PdfPages(output[0]) as pdf:
            text_page(pdf, lines)

            def uvw(key, i):
                t = compare[key]
                j = list(t['combos']).index('uvw')
                return f"{t['mean'][i, j]:.3f} $\\pm$ {t['std'][i, j]:.3f}"

            table_page(pdf, 'U & V & W fraction of true trios at each best-F1 operating point:\n'
                       'mean $\\pm$ std over seeds (lines) or events (cosmics)',
                       [f'{s}\n{kind} F1' for s in sets for kind in OP_KINDS], entries,
                       [[uvw((s, kind), i) for s in sets for kind in OP_KINDS]
                        for i in range(len(entries))])

            for s in sets:
                for kind in OP_KINDS:
                    t = compare[s, kind]
                    cells = []
                    juvw = list(t['combos']).index('uvw')
                    for (style, mode), mean, std, n_valid in zip(config['op_compare'], t['mean'],
                                                                 t['std'], t['n_valid']):
                        op = ops[kind, style, mode]
                        cells.append(['/'.join(f'{x:.2f}' for x in op['thresholds']),
                                      '/'.join(f'{x:.3f}' for x in op['f1'])]
                                     + [f'{x:.3f}' for x in mean] + [f'{std[juvw]:.3f}']
                                     + [str(int(np.sum(n_valid)))])
                    table_page(pdf, f'{s}: best {kind} F1 on {OP_TUNE_SET}, '
                               f'trio rates (mean over {replica_kind(s)})',
                               ['thresh u/v/w', 'tune F1 u/v/w'] + list(t['combos'])
                               + [f'uvw std\n({replica_kind(s)})', 'n_valid'],
                               entries, cells)
                for kind in OP_KINDS:
                    image_page(pdf, f'{res}/op-{kind}F1-trio-compare_{s}.png')

                if not re.fullmatch(LINE_SET, s):
                    continue
                effpur = {kind: np.load(f'{res}/op-{kind}F1-effpur-compare_{s}.npz') for kind in OP_KINDS}

                nseeds = int(effpur[OP_KINDS[0]]['nseeds'])
                cells = []
                for i in range(len(entries)):
                    row = []
                    for kind in OP_KINDS:
                        t = effpur[kind]
                        row += [f'{t["eff_all_mean"][i]:.3f} $\\pm$ {t["eff_all_std"][i]:.3f}',
                                f'{t["pur_all_mean"][i]:.3f} $\\pm$ {t["pur_all_std"][i]:.3f}',
                                f'{np.nanmin(t["eff_mean"][i]):.3f}', f'{np.nanmin(t["pur_mean"][i]):.3f}']
                    cells.append(row)
                table_page(pdf, f'{s}: efficiency and purity, each kind at its own best F1 on {OP_TUNE_SET}\n'
                           f'mean $\\pm$ std over {nseeds} seed(s) (each seed over all angles); '
                           'min over angles of the per-angle mean',
                           [f'{kind}\n{name}' for kind in OP_KINDS
                            for name in ('eff', 'pur', 'min eff', 'min pur')],
                           entries, cells)

                for kind in OP_KINDS:
                    t = effpur[kind]
                    table_page(pdf, f'{s}: {kind} counts per angle summed over {nseeds} seed(s), as '
                               f'efficiency found/true and purity real/reco, at best {kind} F1 on {OP_TUNE_SET}',
                               list(t['angles']), entries,
                               [[f'{t["n_true_matched_sum"][i][j]}/{t["n_true_sum"][i][j]}  '
                                 f'{t["n_reco_matched_sum"][i][j]}/{t["n_reco_sum"][i][j]}'
                                 for j in range(len(t['angles']))]
                                for i in range(len(entries))])
                for kind in OP_KINDS:
                    image_page(pdf, f'{res}/op-{kind}F1-effpur-compare_{s}.png')

            for kind in OP_KINDS:
                for style in op_styles():
                    fig, axs = plt.subplots(2, 2, figsize=(11, 8.5), sharex=True, sharey=True)
                    for ax, plane in zip(axs.flat, SCAN_PLANES):
                        for beta in FBETA_BETAS:
                            f = np.load(f'{res}/fbeta_aggregated_scan_{kind}_plane_{plane}_{style}_{OP_TUNE_SET}_{beta}.npz')
                            ax.plot(f['thresholds'], f['fbeta'],
                                    label=rf'$\beta$={beta:.2f} max at {float(f["maxloc"]):.2f}')
                        ax.set_title(plane)
                        ax.grid(True)
                        ax.legend(fontsize=6)
                    fig.suptitle(rf'{style}: {kind} $F_\beta$ vs threshold on {OP_TUNE_SET}')
                    fig.supxlabel('Threshold')
                    pdf.savefig(fig)
                    plt.close(fig)


#
# Comparison across training epochs
#
def glob_epochs(w):
    cpts = glob(f'{config["cpt_dir"]}/checkpoint_{config["proc"]}_{config["cluster"]}_*.pt')
    # Assumes the epoch number is the last _-separated field of the name.
    epochs = sorted({os.path.basename(c)[:-len('.pt')].split('_')[-1] for c in cpts})
    return [f'epoch_{e}/all_eff_pur_{p}plane-xvu.npz' for e in epochs for p in PLANES]

rule all_epochs_xvu:
    input: glob_epochs

def epoch_results(suffix):
    'Input function: existing per-epoch eff/pur npz files, in epoch order.'
    def inner(w):
        files = glob(f'epoch_*/all_eff_pur_{w.plane}plane{suffix}.npz')
        files = [f for f in files if re.fullmatch(r'epoch_\d+', f.split('/')[0])]
        return sorted(files, key=lambda s: int(s.split('/')[0].split('_')[-1]))
    return inner

rule compare_effpur_epochs:
    input: epoch_results('')
    output:
        pdf='<results>/eff_pur_epochs_{plane}plane.pdf',
    run:
        from matplotlib.backends.backend_pdf import PdfPages
        print(input)
        files = [np.load(f) for f in input]
        epoch_nums = [int(f.split('/')[0].split('_')[-1]) for f in input]
        # [angle, epoch]
        effs = np.array([f['effs'] for f in files]).T
        purs = np.array([f['purs'] for f in files]).T

        with PdfPages(output.pdf) as pdf:
            for i in range(len(effs)):
                plt.figure()
                plt.plot(epoch_nums, effs[i], label='Efficiency')
                plt.plot(epoch_nums, purs[i], label='Purity')
                plt.plot(epoch_nums, effs[i] * purs[i], label='eff x pur')
                plt.ylim(0, 1.2)
                plt.grid(True, axis='both')
                plt.title(f'Angle #{i}')
                plt.xlabel('Epoch')
                plt.legend()
                pdf.savefig()
                plt.close()

            for name, vals in (('Efficiency', effs), ('Purity', purs)):
                plt.figure()
                plt.title(f'{name} -- all epochs')
                plt.xlabel('Angle')
                plt.ylim(0, 1.2)
                plt.ylabel(name)
                plt.grid(True, axis='both')
                for i, v in enumerate(vals.T):
                    plt.plot(v, label=f'Epoch {epoch_nums[i]}')
                plt.legend()
                pdf.savefig()
                plt.close()

            for name, vals in (('Purity', purs), ('Efficiency', effs)):
                plt.figure()
                plt.title(f'{name} -- all epochs')
                img = plt.imshow(vals, aspect='auto', interpolation='none', vmin=0.5, vmax=1.0)
                plt.xlabel('Epoch')
                plt.ylabel('Angle')
                cbar = plt.colorbar(img)
                cbar.ax.set_title(name, pad=10)
                pdf.savefig()
                plt.close()

use rule compare_effpur_epochs as compare_effpur_epochs_xvu with:
    input: epoch_results('-xvu')
    output:
        pdf='<results>/eff_pur_epochs_{plane}plane-xvu.pdf',


#
# Aggregate targets
#
for p in PLANES:
    rule:
        name: f'all_plots_{p}'
        input:
            tar=f'line-displays_{p}plane.tar',
            effpur=f'<results>/all_eff_pur_{p}plane.png',
            npz=f'<results>/all_eff_pur_{p}plane.npz',
    rule:
        name: f'all_plots_{p}_xvu'
        input:
            tar=f'xvu-line-displays_{p}plane.tar',
            effpur=f'<results>/all_eff_pur_{p}plane-xvu.png',
            npz=f'<results>/all_eff_pur_{p}plane-xvu.npz',

rule all_eff_pur_xvu:
    input:
        [f'<results>/all_eff_pur_{p}plane-xvu.npz' for p in PLANES]

rule all_trio_compare:
    input:
        [f'<results>/threshold-0.5-trio-compare_{s}.png'
         for s in [f'{p}plane_high_end' for p in PLANES] + (['cosmics'] if config['cosmics_dir'] else [])]

# Trio comparisons at the best ROI and pixel F1 operating points, tuned on
# cosmics_tune and evaluated on the report_sets, collated by op_report.
rule all_trio_compare_op:
    input:
        '<results>/op-trio-report.pdf'
