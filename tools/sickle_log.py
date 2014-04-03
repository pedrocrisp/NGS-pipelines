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
        fq = None
        if line.startswith("+ sickle"):
            fq = re.search(r'(\S+\.f(q|astq)(\.gz)?)', line).groups()[0]
            sample = path.basename(path.dirname(fq))
            z = next(infh) # gap
            pe_keep = re.search(r'paired.+ \((\S+) pairs\)', next(infh)).groups()[0]
            se_keep = re.search(r'single.+: (\S+) \(', next(infh)).groups()[0]
            pe_shit = re.search(r'paired.+ \((\S+) pairs\)', next(infh)).groups()[0]
            se_shit = re.search(r'single.+: (\S+) \(', next(infh)).groups()[0]
            results[sample] = {'sample': sample, 'pe_keep': pe_keep,
                    'se_keep': se_keep,  'pe_shit': pe_shit,
                    'se_shit': se_shit}
    return results

def print_dict(dct, sep=','):
    header = ['sample', 'pe_keep', 'se_keep', 'pe_shit', 'se_shit', ]
    print(sep.join(header))
    for fq, stats in dct.items():
        line_cells = [stats[st] for st in header]
        print(sep.join(line_cells))


if __name__ == "__main__":
    res = {}
    for fle in sys.argv[1:]:
        res = parse_log(fle, dct=res)
    print_dict(res)
