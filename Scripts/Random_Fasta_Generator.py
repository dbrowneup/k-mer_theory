#!/usr/bin/env python

"""
The purpose of this code is to walk out a random string of DNA
multiple times to construct a Fasta file

python Random_Fasta_Generator.py <Number of Entries> <Entry Length> <Output Name>
"""

import sys
import random

class Analyzer():
    def __init__(self):
        return None

    def __enter__(self):
        return self

    def __exit__(self, exctype, value, traceback):
        if exctype is None:
            return

    def generate_dna(self, L):
        return ''.join(random.choice('ATGC') for _ in range(L))

def main():
    N = int(sys.argv[1])
    L = int(sys.argv[2])
    OUT = open(str(sys.argv[3])+'.fa', 'w')
    with Analyzer() as A:
        for i in range(N):
            S = A.generate_dna(L)
            OUT.write('>'+str(i+1)+'\n'+S+'\n')

if __name__ == '__main__':
    main()





