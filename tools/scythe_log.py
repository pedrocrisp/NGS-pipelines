from __future__ import print_function
import re
import sys


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
        if line.startswith("+ scythe"):
            fq = re.search(r'(\S+\.f(q|astq)(\.gz)?)', line).groups()[0]
            prior = re.search(r'prior: (\S+)', next(infh)).groups()[0]
            _ = next(infh) # gap
            _ = next(infh) # adap
            line = next(infh)
            cont = re.search(r'contaminated: (\S+),', line).groups()[0]
            notcont = re.search(r'uncontaminated: (\S+),', line).groups()[0]
            total = re.search(r'total: (\S+)', line).groups()[0]
            rate = re.search(r'rate: (\S+)', next(infh)).groups()[0]
            results[fq] = {'fq': fq, 'cont': cont, 'uncont': notcont,
                    'total': total, 'rate':rate}
    return results


def print_dict(dct, sep=','):
    header = ['fq', 'cont', 'uncont', 'total', 'rate']
    print(sep.join(header))
    for fq, stats in dct.items():
        line_cells = [stats[st] for st in header]
        print(sep.join(line_cells))


if __name__ == "__main__":
    res = {}
    for fle in sys.argv[1:]:
        res = parse_log(fle, dct=res)
    print_dict(res)
