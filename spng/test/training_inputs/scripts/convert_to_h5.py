import h5py as h5, torch, numpy as np

subsets = ['frame', 'tickinfo', 'channels'] #Need these for a check?

infile = snakemake.input[0]
if '.pt' in infile:
    t = torch.load(infile, weights_only=False)
    print(t)
elif '.npz' in infile:
    #Open the file
    t = np.load(infile)
    print(t)
    requested_tags = snakemake.params.tags
    all_events = dict()
    with h5.File(snakemake.output[0], 'w') as fout:
        for k in t.keys():
            print(k)
            name, tag, event_num = k.split('_')
            if tag not in requested_tags.keys(): continue
            # event_num = k.split('_')[-1]
            if event_num not in all_events:
                all_events[event_num] = fout.create_group(event_num)
            event_group = all_events[event_num]
            tag_out = requested_tags[tag]
            # name = k.split('_')[0]

            label = f'{name}_{tag_out}0'
            event_group.create_dataset(label, data=t[k])