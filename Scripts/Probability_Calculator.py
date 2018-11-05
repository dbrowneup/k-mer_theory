#!/usr/bin/env python

"""
The purpose of this code is to walk out a random string of DNA
and then analyze the resulting k-mer frequency distribution

python Probability_Calculator.py <string length>
"""

import sys
import random
from math import log

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


class Analyzer():
    def __init__(self):
        self.kmers = dict()
        self.probability = list()

    def __enter__(self):
        return self

    def __exit__(self, exctype, value, traceback):
        if exctype is None:
            return

    def generate_dna(self, L):
        return ''.join(random.choice('ATGC') for _ in range(L))

    def clear_kmers(self):
        self.kmers = dict()
        return

    def count_kmers(self, S, k):
        for i in range(len(S) - k + 1):
            s = S[i:i + k]
            try:
                self.kmers[s] += 1
            except KeyError:
                self.kmers[s] = 1

    def kmer_scan(self, S):
        for i in range(len(S)):
            k = i + 1
            self.count_kmers(S, k)
            O = len(self.kmers)
            Pr  = log(float(1) / O, 10)
            self.probability.append((k, Pr))
            self.clear_kmers()
            N = len(S) - k + 1
            if N == O:
                return k

def main():
    L = int(sys.argv[1])
    with Analyzer() as A:
        S = A.generate_dna(L)
        R = A.kmer_scan(S)
        k, Pr = zip(*A.probability)
        plt.plot(k, Pr)
        plt.figtext(0.70, 0.75, 'R = '+str(R))
        plt.xlabel('k')
        plt.ylabel('Pr(S, k)')
        plt.savefig('Probability_Output_L'+str(L)+'.png')

if __name__ == '__main__':
    main()





