import re
import sys

def parse_log(infn='-'):
    # input and output files, accept stdin for both.
    if infn in {"-", "stdin"}:
        infh = sys.stdin
    else:
        infh = open(infn)
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

def print_dict(dct, outfn='-', sep=','):
    if outfn in {"-", "stdout"}:
        outfh = sys.stdout
    else:
        outfh = open(outfn)
    header = ['fq', 'cont', 'uncont', 'total', 'rate']
    outfh.write(sep.join(header) + '\n')
    for fq, stats in dct.items():
        line_cells = [stats[st] for st in header]
        outfh.write(sep.join(line_cells) + '\n')

if __name__ == "__main__":
    res = parse_log(sys.argv[1])
    try:
        print_dict(res, sys.argv[2], sys.argv[3])
    except IndexError:
        try:
            print_dict(res, sys.argv[2])
        except IndexError:
            try:
                print_dict(res)
            except:
                pass
