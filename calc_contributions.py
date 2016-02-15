#!/usr/bin/env python
# author          : Pavel
# date            : 15.02.16
# version         : 0.1
# python_version  : 3
# copyright       : Pavel 2016
# license         : GPL3
#==============================================================================

import argparse
from collections import OrderedDict


def main_params(in_fname, out_fname, sep):

    # read input data (2 columns, each line: mol_name###frag_name 1.25)
    d = OrderedDict()
    with open(in_fname) as f:
        for line in f:
            tmp = line.strip().split()
            names = tmp[0].split(sep)
            if names[0] not in d.keys():
                d[names[0]] = dict()
            if len(names) == 1:
                d[names[0]][''] = float(tmp[1])
            else:
                d[names[0]][names[1]] = float(tmp[1])

    with open(out_fname, 'wt') as f:
        for mol_name, v in d.items():
            for frag_name, pred in v.items():
                if frag_name == '':
                    continue
                else:
                    f.write(mol_name + sep + frag_name + '\t' + str(v[''] - pred) + '\n')


def main():

    parser = argparse.ArgumentParser(description='Calculate contributions of molecular fragments based predicted '
                                                 'values from QSAR models.')
    parser.add_argument('-i', '--input', metavar='predictions.txt', required=True,
                        help='input text file with predicted values for initial compounds and compounds with '
                             'removed fragments.')
    parser.add_argument('-o', '--output', metavar='contributions.txt', required=True,
                        help='output text file with fragments contributions.')
    parser.add_argument('-s', '--sep', metavar='', required=False, default='###',
                        help='separator used for separation of molecule name and fragment name. Default: ###.')

    args = vars(parser.parse_args())
    for o, v in args.items():
        if o == "input": in_fname = v
        if o == "output": out_fname = v
        if o == "sep": sep = v

    main_params(in_fname, out_fname, sep)


if __name__ == '__main__':
    main()
