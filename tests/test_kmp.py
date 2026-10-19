import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parent.parent))


import numpy as np

from src.preprocessor.kmp import compute_kmp_pi, kmp_suffix_prefix


def test_compute_kmp_pi_basic():
    """Test pi table construction for a standard repeating pattern."""
    pattern = "ABABA"
    pi = compute_kmp_pi(pattern)

    # Expected pi table of size n + 1 (length 6) for "ABABA"
    # pi[0] = 0, pi[1]=0, pi[2]=0, pi[3]=1, pi[4]=2, pi[5]=3 (depending on implementation indexing)
    assert isinstance(pi, np.ndarray)
    assert len(pi) == len(pattern) + 1
    # Check specific failure transitions
    assert pi[1] == 0  # 'A'
    assert pi[3] == 1  # prefix 'AB' in 'ABA' -> overlap 1 ('A')
    assert pi[5] == 3  # prefix 'ABA' in 'ABABA' -> overlap 3 ('ABA')


def test_kmp_suffix_prefix_matching():
    """Test direct suffix-prefix matching using KMP functions."""
    pattern = "TGGAATTCTC"
    pi = compute_kmp_pi(pattern)

    # Case A: Read ending in a 4-bp prefix of the pattern
    read_a = "ACGTACGT" + "TGGA"
    assert kmp_suffix_prefix(read_a, pattern, pi) == 4

    # Case B: Full match (read equals pattern)
    assert kmp_suffix_prefix(pattern, pattern, pi) == len(pattern)

    # Case C: No overlap at all
    read_c = "CCCCCCCC"
    assert kmp_suffix_prefix(read_c, pattern, pi) == 0

    # Case D: Empty read
    assert kmp_suffix_prefix("", pattern, pi) == 0
