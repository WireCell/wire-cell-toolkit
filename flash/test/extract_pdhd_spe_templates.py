#!/usr/bin/env python3
"""One-time extraction of the PDHD SPE templates and the
channel->template map from the duneopdet product on cvmfs into JSON
config data consumed by Flash::OpDecon.

Sources (duneopdet v10_20_09d00):
  config_data/SPE_NP04_FBK_2024_without_pretrigger.dat  (template 0)
  config_data/SPE_NP04_HPK_2024_without_pretrigger.dat  (template 1)
  fcl/dune_opdet_channels.fcl protodunehd_pds_channels_mc channel lists
  (the 2024-era FBK/HPK split of the 160 offline OpChannels).

Usage:
  python3 extract_pdhd_spe_templates.py <duneopdet_dir> <outfile>
"""
import json
import sys

# protodunehd_pds_channels_mc, dune_opdet_channels.fcl
FBK_CHANNELS = [
    4, 14, 24, 34, 40, 42, 45, 46, 47, 49, 50, 52, 55, 56, 57, 59, 60, 62, 65, 66,
    67, 69, 70, 72, 75, 76, 77, 79, 84, 85, 86, 87, 94, 95, 96, 97, 104, 105, 106,
    107, 114, 115, 116, 117, 120, 121, 124, 125, 127, 129, 130, 131, 134, 135, 137,
    139, 140, 141, 144, 145, 147, 149, 150, 151, 154, 155, 157, 159,
]
HPK_CHANNELS = [
    0, 1, 2, 3, 5, 6, 7, 8, 9, 10, 11, 12, 13, 15, 16, 17, 18, 19, 20, 21, 22, 23,
    25, 26, 27, 28, 29, 30, 31, 32, 33, 35, 36, 37, 38, 39, 41, 43, 44, 48, 51, 53,
    54, 58, 61, 63, 64, 68, 71, 73, 74, 78, 80, 81, 82, 83, 88, 89, 90, 91, 92, 93,
    98, 99, 100, 101, 102, 103, 108, 109, 110, 111, 112, 113, 118, 119, 122, 123,
    126, 128, 132, 133, 136, 138, 142, 143, 146, 148, 152, 153, 156, 158,
]

TEMPLATES = [
    "SPE_NP04_FBK_2024_without_pretrigger.dat",
    "SPE_NP04_HPK_2024_without_pretrigger.dat",
]


def main(duneopdet_dir, outfile):
    templates = []
    for fname in TEMPLATES:
        with open(f"{duneopdet_dir}/config_data/{fname}") as fp:
            values = [float(line.split()[0]) for line in fp if line.strip()]
        templates.append(dict(name=fname.replace(".dat", ""), values=values))
        print(f"{fname}: {len(values)} samples, max {max(values):.3f}")

    channels = FBK_CHANNELS + HPK_CHANNELS
    template_index = [0] * len(FBK_CHANNELS) + [1] * len(HPK_CHANNELS)
    assert sorted(channels) == list(range(160))

    out = {
        "comment": "PDHD single-p.e. response templates (ADC*us, 16 ns ticks) and the "
                   "offline-OpChannel -> template map, for Flash::OpDecon. Extracted from "
                   "duneopdet config_data (2024 NP04 FBK/HPK templates) with the channel "
                   "split of protodunehd_pds_channels_mc (dune_opdet_channels.fcl).",
        "source": duneopdet_dir,
        "templates": templates,
        "channels": channels,
        "template_index": template_index,
    }
    with open(outfile, "w") as fp:
        json.dump(out, fp, indent=1)
    print(f"wrote {outfile} ({len(channels)} channels, {len(templates)} templates)")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
