#!/usr/bin/env python3
"""One-time extraction of the PDHD run28368 v1 per-channel SPE
templates, the channel->template map and the run27950 FFT noise
templates into the JSON config data consumed by Flash::OpDecon.

This is the template set actually used by the LArSoft production that
made the reference deconv in pdhd/example_light_data (dunesw
v10_20_09d00 protodunehd_deconvolution -> protodunehd_pds_channels_data_v1
-> protodunehd_template_list_v1).

Sources (a directory holding copies of):
  v1/run28368_*_Change_Scale_PDHD_Jun2025.txt     113 SPE templates from
      /cvmfs/dune.osgstorage.org/pnfs/fnal.gov/usr/dune/persistent/stash/ProtoDUNE/HD/opdetresponse/v1/
  FFT_Noise_Template_run27950_PDHD_VGain2318_sample513_{fbk,hpk}_Dec2024.txt
      from .../stash/ProtoDUNE/HD/opdetresponse/v0/
  protodunehd_template_list_v1.fcl                SPETemplateFiles order
  protodunehd_first_template_list.fcl             channel maps (SPE + noise)
      from duneopdet v10_20_09d00 fcl/

Usage:
  python3 extract_pdhd_spe_templates_v1.py <source_dir> <spe_out.json> <noise_out.json>
"""
import json
import re
import sys


def fcl_string_list(text, key):
    block = re.search(key + r"\s*:\s*\[(.*?)\]", text, re.S).group(1)
    return re.findall(r'"([^"]+)"', block)


def fcl_int_list(text, key):
    block = re.search(key + r"\s*:\s*\[(.*?)\]", text, re.S).group(1)
    block = re.sub(r"#[^\n]*", "", block)
    return [int(tok) for tok in block.replace(",", " ").split()]


def read_column(path):
    with open(path) as fp:
        return [float(line.split()[0]) for line in fp if line.strip()]


def main(srcdir, spe_out, noise_out):
    with open(f"{srcdir}/protodunehd_template_list_v1.fcl") as fp:
        v1 = fp.read()
    with open(f"{srcdir}/protodunehd_first_template_list.fcl") as fp:
        first = fp.read()

    # v1 overrides SPETemplateFiles only; the channel map (same file
    # order) and the noise tables come from the first_template_list.
    spe_files = fcl_string_list(v1, "SPETemplateFiles")
    spe_channels = fcl_int_list(first, "SPETemplateMapChannels")
    spe_index = fcl_int_list(first, "SPETemplateMapTemplates")
    assert len(spe_channels) == len(spe_index) == len(spe_files)

    templates = []
    for fname in spe_files:
        values = read_column(f"{srcdir}/v1/{fname}")
        templates.append(dict(name=fname.replace(".txt", ""), values=values))
    print(f"SPE: {len(templates)} per-channel templates, "
          f"max amplitudes {min(max(t['values']) for t in templates):.2f}"
          f"..{max(max(t['values']) for t in templates):.2f} ADC")

    out = {
        "comment": "PDHD run28368 v1 per-channel single-p.e. response templates (ADC, "
                   "16 ns ticks) and the offline-OpChannel -> template map, for "
                   "Flash::OpDecon. The production set of protodunehd_template_list_v1 "
                   "(duneopdet v10_20_09d00); files from the DUNE StashCache "
                   "ProtoDUNE/HD/opdetresponse/v1/. Unmapped channels (dead 86,87,97,"
                   "107,116,117, noisy 3, full-stream 120-159) are skipped, as in "
                   "LArSoft IgnoreChannels.",
        "source": "ProtoDUNE/HD/opdetresponse/v1 via protodunehd_template_list_v1.fcl",
        "templates": templates,
        "channels": spe_channels,
        "template_index": spe_index,
    }
    with open(spe_out, "w") as fp:
        json.dump(out, fp, indent=1)
    print(f"wrote {spe_out} ({len(spe_channels)} channels, {len(templates)} templates)")

    noise_files = fcl_string_list(first, "NoiseTemplateFiles")
    noise_channels = fcl_int_list(first, "NoiseTemplateMapChannels")
    noise_index = fcl_int_list(first, "NoiseTemplateMapTemplates")
    assert len(noise_channels) == len(noise_index)
    noise_templates = []
    for fname in noise_files:
        values = read_column(f"{srcdir}/{fname}")
        noise_templates.append(dict(name=fname.replace(".txt", ""), values=values))
        print(f"noise {fname}: {len(values)} bins")

    out = {
        "comment": "PDHD run27950 noise power spectra (|FFT|^2 per half-spectrum bin, "
                   "same normalization as LineNoiseRMS^2*Samples) and the "
                   "offline-OpChannel -> template map, for Flash::OpDecon N^2. From "
                   "protodunehd_first_template_list (duneopdet v10_20_09d00); files "
                   "from the DUNE StashCache ProtoDUNE/HD/opdetresponse/v0/. Read to "
                   "Samples/2+1 bins, missing tail bins zero-padded as in LArSoft.",
        "source": "ProtoDUNE/HD/opdetresponse/v0 via protodunehd_first_template_list.fcl",
        "templates": noise_templates,
        "channels": noise_channels,
        "template_index": noise_index,
    }
    with open(noise_out, "w") as fp:
        json.dump(out, fp, indent=1)
    print(f"wrote {noise_out} ({len(noise_channels)} channels, {len(noise_templates)} templates)")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3])
