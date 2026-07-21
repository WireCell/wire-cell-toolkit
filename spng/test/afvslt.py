#!/usr/bin/env python

import sys
import json

from collections import defaultdict

infinite_defaultdict = lambda: defaultdict(infinite_defaultdict)

one = infinite_defaultdict()

try:
    fname = sys.argv[1]
except IndexError:
    fname = "/dev/stdin"

for tst in json.loads(open(fname).read()):
    inp = tst["input"]

    one[inp["kind"]][inp["device"]][inp["tech"]] = tst["output"]["dt"]

print('|test\t|dev\t|tech\t|time (ms)|')
for l1, l1_rest in one.items():
    for l2, l2_rest in l1_rest.items():
        for l3, l3_rest in l2_rest.items():
            print(f'|{l1}\t|{l2.upper()}\t|{l3.upper()}\t|{l3_rest:8.2f} |')



