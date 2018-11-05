#!/usr/bin/env python

"""
The purpose of this code is to walk out a random string of DNA
and then analyze the resulting k-mer frequency distribution

python Random_String_Analyzer.py <string length> <kmer length>
"""

import sys
import random

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


class Analyzer():
    def __init__(self):
        self.kmers = dict()

    def __enter__(self):
        return self

    def __exit__(self, exctype, value, traceback):
        if exctype is None:
            return

    def generate_dna(self, L):
        return ''.join(random.choice('ATGC') for _ in range(L))

    def count_kmers(self, S, k):
        for i in range(len(S) - k + 1):
            s = S[i:i + k]
            try:
                self.kmers[s] += 1
            except KeyError:
                self.kmers[s] = 1


def main():
    L = int(sys.argv[1])
    k = int(sys.argv[2])
    with Analyzer() as A:
        S = A.generate_dna(L)
        A.count_kmers(S, k)
        D = A.kmers.values()
        N = L - k + 1
        O = len(A.kmers)
        P = 4 ** k
        plt.hist(D, bins=max(max(D), 10) - min(D), range=(min(D), max(max(D), 10)))
        plt.figtext(0.60, 0.85, str(N)+' total '+str(k)+'-mers')
        plt.figtext(0.60, 0.80, str(O)+' distinct '+str(k)+'-mers')
        plt.figtext(0.60, 0.75, str(P)+' possible '+str(k)+'-mers')
        plt.xlabel(str(k)+'-mer Count')
        plt.ylabel('Distinct '+str(k)+'-mers')
        plt.savefig('Analyzer_Output_L'+str(L)+'_k'+str(k)+'.png')

if __name__ == '__main__':
    main()





