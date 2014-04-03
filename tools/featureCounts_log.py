from __future__ import print_function
import re
import sys
from os import path


def parse_log(infn='-', dct=None):
    if infn in {"-", "stdin"}:
        infh = sys.stdin
    else:
        infh = open(infn)
    if dct is not None:
        results = dct
    else:
        results = {}
    for line in infh:
        if line.startswith("+ sample="):
            sample = line.split('=')[1].strip()
        if line.startswith("||    Number of features"):
            features = re.search(r'is (\S+)', line).groups()[0]
        if line.startswith("||    Total number of fragments"):
            frag = re.search(r' : (\S+)', line).groups()[0]
        if line.startswith("||    Number of successfully assigned fragments"):
            mapd, pct = re.search(r' : (\S+) \((\S+)%\)', line).groups()
            results[sample] = {'sample': sample, 'frag': frag, "mapped":
                    mapd, 'percent_mapped': pct}
    return results

def print_dict(dct, sep=','):
    header = ['sample', 'frag', 'mapped', 'percent_mapped']
    print(sep.join(header))
    for fq, stats in dct.items():
        line_cells = [stats[st] for st in header]
        print(sep.join(line_cells))


if __name__ == "__main__":
    res = {}
    for fle in sys.argv[1:]:
        res = parse_log(fle, dct=res)
    print_dict(res)
