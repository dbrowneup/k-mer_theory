#!/usr/bin/env python

"""
The purpose of this code is to walk out a random string of DNA
and then analyze the resulting k-mer frequency distribution

python Bio_Probability_Scanner_Limited_Log.py <Fasta File> <Window Length> <Output Name>

This version is designed to scan all the windows of given length in a Fasta file
containing one or more sequences, using a single-threaded approach

Requires:
pyfaidx
NumPy
matplotlib
"""

import sys
import random
from math import log

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from pyfaidx import Fasta


class Analyzer():
    def __init__(self):
        self.kmers = dict()
        self.probability = list()

    def __enter__(self):
        return self

    def __exit__(self, exctype, value, traceback):
        if exctype is None:
            return

    def clear_kmers(self):
        self.kmers = dict()
        return

    def count_kmers(self, S, k):
        U = 0
        for i in range(len(S) - k + 1):
            s = S[i:i + k]
            if 'N' in s:
                U += 1
                continue
            try:
                self.kmers[s] += 1
            except KeyError:
                self.kmers[s] = 1
        return U

    def kmer_scan(self, S):
        for i in range(len(S)):
            k = i + 1
            U = self.count_kmers(S, k)
            O = len(self.kmers)
            if O == 0:
                self.probability.append((k, 0))
                return (k, 'gap_stop')
            Pr  = float(1) / O
            self.probability.append((k, log(Pr, 10)))
            self.clear_kmers()
            N = len(S) - k + 1
            if N == O + U:
                return (k - 1, 'complete')

def main():
    F = Fasta(sys.argv[1])
    W = int(sys.argv[2])
    OUT = str(sys.argv[3])
    R_list = list()
    R_dict = {'complete': list(), 'gap_stop': list()}
    P_dict = dict()
    with Analyzer() as A:
        for f in F:
            S = str(f)
            if len(S) < W:
                continue
            for i in range(len(S[::W])):
                x = i * W
                y = (i + 1) * W
                s = S[x:y]
                if len(s) < W:
                    continue
                R = A.kmer_scan(s)
                R_list.append(R)
                for k, Pr in A.probability:
                    try:
                        P_dict[k].append(Pr)
                    except KeyError:
                        P_dict[k] = [Pr]
                A.probability = list()
        
        x, y = zip(*P_dict.iteritems())
        y = np.array([[min(Pr), np.mean(Pr) - np.std(Pr), np.mean(Pr), np.mean(Pr) + np.std(Pr), max(Pr)] for Pr in y])
        plt.semilogx(x, y[:,0], 'k-')
        plt.semilogx(x, y[:,1], 'r-')
        plt.semilogx(x, y[:,2], 'b-')
        plt.semilogx(x, y[:,3], 'r-')
        plt.semilogx(x, y[:,4], 'k-')
        plt.xlabel('k')
        plt.ylabel('log10(Pr(S, k))')
        plt.savefig('Pr_Distribution_'+OUT+'_W'+str(W)+'.png')
        plt.clf()
        
        for score, status in R_list:
            R_dict[status].append(score)
        mn = min(R_dict['complete'])
        mx = max(R_dict['complete'])
        plt.hist(R_dict['complete'], bins=mx - mn, range=(mn, mx), alpha=0.5, label='complete')
        plt.legend(loc='upper right')
        plt.xlabel('R')
        plt.ylabel('Number of Windows')
        plt.savefig('R_Distribution_'+OUT+'_W'+str(W)+'.png')


if __name__ == '__main__':
    main()
