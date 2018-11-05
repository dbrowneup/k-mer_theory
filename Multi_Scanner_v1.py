#!/usr/bin/env python

"""
The purpose of this code is to:
1) walk out a random string of DNA of length L
2) determine k-mer frequencies in windows of length W
3) load k-mer frequencies into a dictionary
4) analyze k-mer distributions across windows

python Multi_Scanner_v1.py <k-mer size> <string length> <window length>
"""

import sys
import random

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import numpy as np
import pandas as pd

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

    def clear_kmers(self):
        self.kmers = dict()
        return

    def count_kmers(self, S, k):
        for i in range(len(S) - k + 1):
            kmer = S[i:i + k]
            try:
                self.kmers[kmer] += 1
            except KeyError:
                self.kmers[kmer] = 1

def main():
    k = int(sys.argv[1])
    L = int(sys.argv[2])
    W = int(sys.argv[3])
    D = dict()
    with Analyzer() as A:
        DNA = A.generate_dna(L)
        for i in range(len(DNA[::W])):
            x = i * W
            y = (i + 1) * W
            S = DNA[x:y]
            A.count_kmers(S, k)





if __name__ == '__main__':
    main()
