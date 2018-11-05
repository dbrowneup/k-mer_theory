import decimal
import math
import sys

import matplotlib
import matplotlib.pyplot as plt

"""
What is the probability of a k-mer m occuring C times in a random sample of N k-mers?

lambda is the average number of events in an interval, otherwise known as the mean 
count of all k-mers, roughly approximated by 1 + N / P for a random sequence.

K is the number of times a k-mer m is observed

Thus for each value of k, the Poisson model should be able to approximate the 
empirical distribution of k-mer counts in random sequences. This can be calculated 
for each value of k and fitted to experimental data to determine the predictive 
power of the model.

usage: python Poisson_Calculator.py <L (int)> <k (int)>
"""

def Poisson(lam, K):
    return decimal.Decimal((lam ** K)) * decimal.Decimal((math.e ** -lam)) / decimal.Decimal(math.factorial(K))

def main():
    L = int(sys.argv[1])
    k = int(sys.argv[2])
    N = L - k + 1
    P = 4 ** k
    prob = list()
    lam = 1 + N / P
    if lam > 20:
        for K in range(int(math.ceil(2 * lam))):
            prob.append((K, Poisson(lam, K)))
    else:
        for K in range(20):
            prob.append((K, Poisson(lam, K)))
    prob = [(x, y) for x, y in prob if y > 0]
    x, y = zip(*prob)
    plt.plot(x, y, 'bo')
    plt.savefig('Poisson_Probability_L'+str(L)+'_k'+str(k)+'_v1.png')

if __name__ == '__main__':
	main()