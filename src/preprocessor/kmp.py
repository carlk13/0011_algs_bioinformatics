import numpy as np


def compute_kmp_pi(P: str) -> np.ndarray:
    """
    Computes the KMP failure table (pi table) of size n + 1
    """
    n = len(P)
    pi = np.zeros(n + 1, dtype=int)
    k = 0

    for i in range(1, n):
        while k > 0 and P[k] != P[i]:
            k = pi[k]
        if P[k] == P[i]:
            k = k + 1
        pi[i + 1] = k

    return pi


def kmp_suffix_prefix(read: str, P: str, pi: np.ndarray) -> int:
    """
    Suffix-prefix matcher using the curriculum's pi table of size n + 1.
    """
    n = len(P)
    q = 0  # state

    for c in read:
        while q > 0 and (q == n or P[q] != c):
            q = pi[q]
        if P[q] == c:
            q = q + 1

    return q
