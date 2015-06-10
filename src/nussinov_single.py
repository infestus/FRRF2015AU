#! /usr/bin/python3

from Bio import SeqIO
# import argparse
import numpy
import sys

debug = 0
numpy.set_printoptions(suppress=False, threshold = 5000, linewidth = 2000)
def nussinov(seq, verbose=False):
    mat = numpy.full(shape=(len(seq), len(seq)), fill_value=-1)
    back = numpy.zeros(shape=(len(seq), len(seq)))
    def score(i, j):
        if (j-2 < i):
            return 0
        if (mat[i, j] == -1):
            best = 0
            source = 0
            for k in range(i+1, j-1):
                if score(i, k) + score(k+1, j) > best:
                    best = score(i, k) + score(k+1, j)
                    source = k
            if permitted_match(seq[i], seq[j]) and score(i+1, j-1)+1 > best:
                global debug
                debug = debug +1
                best = score(i+1, j-1)+1
                source = -1
            if score(i+1, j) > best:
                best = score(i+1, j)
                source = -2
            if score(i, j-1) > best:
                best = score(i, j-1)
                source = -3
            mat[i, j] = best
            back[i, j] = source
        return mat[i, j]
    def backtrack(i, j):
        if (j < i):
            return ''
        source = back[i, j]
        if source > 0:
            return backtrack(i, source) + backtrack(source+1, j)
        if source == -1:
            return '(%s)' % backtrack(i+1, j-1)
        if source == -2:
            return '.%s' % backtrack(i+1, j)
        if source == -3:
            return '%s.' % backtrack(i, j-1)
        return '.'*(j-i+1)

    s = score(0, len(seq)-1)
    if verbose:
        print(mat)
    return (s, backtrack(0, len(seq)-1))


def permitted_match(a, b):
    if b < a:
        b, a = a, b
    return (a, b) in [('A', 'U'), ('C', 'G'), ('G', 'U')]


if __name__ == '__main__':
    # parser = argparse.ArgumentParser(description='Nussinov implementations')
    # parser.add_argument('fasta', type=str, help='Fasta file location')
    # parser.add_argument("--verbose", action='store_true', help="Pretty print of matrix")
    # parser.add_argument("--vienna", "-v", action='store_true', help="Print the calculated 2nd structure")
    # args = parser.parse_args()
    # seqs = SeqIO.parse(open(args.fasta), "fasta")
    # for seq in seqs:
    if (len(sys.argv) <= 1):
        seq = "GAUAUGCACUUUUGUACGUCGCUUUGGAUAAAAGCGUCUGCGAAAUAAAUGUAA"
        score, alignment = nussinov(seq.upper(), True)
        print('>%s score %d\n%s\n%s' % ("Sequence", score, seq, alignment))
    else:
        seq = sys.argv[1] 
        seq = seq.upper()
        score = nussinov(seq)
        print(int(score[0]))
    # if args.vienna:
    # else:
    #     print('%-15s %d' % ("Sequence", score))
