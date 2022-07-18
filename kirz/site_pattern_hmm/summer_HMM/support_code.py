import sys
import numpy as np
try:
    logaddexp = np.logaddexp
except AttributeError:
    def logaddexp(logx, logy):
        if logy - logx > 100:
            return logy
        elif logx - logy > 100:
            return logx
        minxy = min(logx, logy)
        return minxy + np.log(np.exp(logx - minxy) + np.exp(logy - minxy))


# SUPPORT CODE
# Stochastic Parameters
# def initialize():
#     std = 0.05
#
#     A = np.array([[np.random.normal(0.7, std), np.random.normal(0.3, std)],
#                   [np.random.normal(0.4, std), np.random.normal(0.6, std)]])
#
#     B = np.array([[np.random.normal(0.3, std), np.random.normal(0.2, std),
#                    np.random.normal(0.2, std), np.random.normal(0.3, std)],
#                   [np.random.normal(0.2, std), np.random.normal(0.3, std),
#                    np.random.normal(0.3, std), np.random.normal(0.2, std)]])
#     pi = [np.random.normal(0.8, std), np.random.normal(0.2, std)]
#
#     A /= np.tile(np.sum(A, axis=1), (2, 1)).T
#     B /= np.tile(np.sum(B, axis=1), (4, 1)).T
#     pi /= np.sum(pi)
#     return A, B, pi

# Viterbi Algorithm (Most likely sequence)
def viterbi(A, B, pi, Ob, N, T):
    V = np.zeros((T, N))
    tb = np.zeros((T, N), dtype='int')

    for i in range(N):
        V[0, i] = np.log(pi[i]) + np.log(B[i, Ob[0]])
        tb[0, i] = -1

    for t in range(1, T):
        for j in range(N):
            V[t, j] = np.NINF
            for i in range(N):
                p = V[t - 1, i] + np.log(A[i, j]) + np.log(B[j, Ob[t]])
                if p > V[t, j]:
                    V[t, j] = p
                    tb[t, j] = i

    logP_max = np.NINF

    for j in range(N):
        if V[-1, j] > logP_max:
            logP_max = V[-1, j]
            Q = [j]
            i = j
    for t in range(T - 2, -1, -1):
        i = tb[t + 1, i]
        Q.insert(0, i)

    return Q, logP_max


def compute_logP(alpha):
    if alpha[-1][0] >= 0:
        return np.log(np.sum(alpha[-1]))
    else:
        return np.log(np.sum(np.exp(alpha[-1])))


def print_results(A, B, pi, Ob, N, T, logP_new):
    Q, logP_max = viterbi(A, B, pi, Ob, N, T)
    print('logP(O|lambda): {0:.2f}'.format(logP_new),
          ''.join([str(i) for i in Q]),
          'logP(O|Q*): {0:.2f}'.format(logP_max))

def logsum(matrix):
    #Implement logsum for a matrix object
    if len(matrix.shape) > 1:
        vec = np.reshape(matrix, (np.product(matrix.shape),))
    else:
        vec = matrix
    # TODO: SUM IS WHATEVER LOG 0 IS
    sum = np.NINF
    for num in vec:
        sum = logaddexp(sum, num)
    return sum

def exp_logsum(numbers):
    #Return the exponential of a logsum
    sum = logsum(numbers)
    return np.exp(sum)