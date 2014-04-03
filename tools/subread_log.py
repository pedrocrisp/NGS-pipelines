from __future__ import print_function
import re
import sys
from os import path


def parse_log(infn='-', dct=None):
    # input and output files, accept stdin for both.
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
        if re.search("==== Summary ====", line):
            z = next(infh) # gap
            pairs = re.search(r'Processed : (\S+)', next(infh)).groups()[0]
            mapd, pct = re.search(r'Mapped : (\S+) \S+ \((\S+)%\)', next(infh)).groups()
            results[sample] = {'sample': sample, 'pairs': pairs, "mapped":
                    mapd, 'percent_mapped': pct}
    return results

def print_dict(dct, sep=','):
    header = ['sample', 'pairs', 'mapped', 'percent_mapped']
    print(sep.join(header))
    for fq, stats in dct.items():
        line_cells = [stats[st] for st in header]
        print(sep.join(line_cells))


if __name__ == "__main__":
    res = {}
    for fle in sys.argv[1:]:
        res = parse_log(fle, dct=res)
    print_dict(res)
