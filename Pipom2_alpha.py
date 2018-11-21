# -*- coding: utf-8 -*-
from __future__ import print_function
import argparse
import regex


def read_fastaK(files):
    x = open(files)
    key = None
    seq = ""
    out = dict()
    for line in x:
        line = line.strip()
        if key is None and line.startswith(">"):
            key = line[1:]
        elif not(key is None) and line.startswith(">"):
            out[key] = seq
            seq = ""
            key = line[1:]
        else:
            seq += line.upper()
    out[key] = seq
    return(out)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--fasta_ref', dest='fasta_ref', required=True, help='subject')
    parser.add_argument('--fasta_target', dest='fasta_target', required=True, help='subject')
    parser.add_argument('--n_mismatch', dest='n_mismatch', default='1', help='output file')
    parser.add_argument('--out', dest='out_file', default='Pipom_result.txt', help='output file')

    args = parser.parse_args()    

    fasta_ref = read_fastaK(args.fasta_ref)
    fasta_target = read_fastaK(args.fasta_target)
    n_mismatch = str(args.n_mismatch)

    with open('out.txt', 'w') as out:
        for target_annot, target_seq in fasta_target.items():
            pattern = "(" + target_seq + "){e<=" + n_mismatch + "}"
            for ref_annot, ref_seq in fasta_ref.items():
                match = regex.findall(pattern, ref_seq)
                if match:
                    print(target_annot, target_seq, ref_annot.split(' ')[0], file=out)


if __name__ == '__main__':
    main()
