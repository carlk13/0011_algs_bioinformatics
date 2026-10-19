import sys
from pathlib import Path

import pytest

sys.path.append(str(Path(__file__).resolve().parent.parent))

from src.preprocessor.exact_matcher import ExactAdapterTrimmer

# Standard test adapter
TEST_ADAPTER = "TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACACAGTGATCTCGTATGCCGTCTTCTGCTTG"


@pytest.fixture
def trimmer():
    return ExactAdapterTrimmer(adapter_sequence=TEST_ADAPTER, min_overlap=1)


def test_1_empty_read(trimmer):
    """Test 1: Empty read returns overlap 0 and is not trimmed."""
    trimmed_seq, L, is_trimmed = trimmer.trim_read("")
    assert L == 0
    assert is_trimmed is False
    assert trimmed_seq == ""


def test_2_no_overlap(trimmer):
    """Test 2: Read completely unrelated to the adapter returns overlap 0."""
    read = "AAAAAA"
    trimmed_seq, L, is_trimmed = trimmer.trim_read(read)
    assert L == 0
    assert is_trimmed is False
    assert trimmed_seq == read


def test_3_one_bp_overlap(trimmer):
    """Test 3: Read ending in 'T' (the first character of the adapter) yields 1-bp overlap."""
    # Since TEST_ADAPTER starts with 'T', a read ending in 'T' has a 1-bp suffix-prefix overlap
    read = "ACGTACGTACGT"  # Ends with 'T'
    trimmed_seq, L, is_trimmed = trimmer.trim_read(read)
    assert L == 1
    assert is_trimmed is True
    assert trimmed_seq == read[:-1]


def test_4_full_read_is_adapter_prefix(trimmer):
    """Test 4: Read that matches a prefix of the adapter perfectly."""
    # Take the first 20 characters of the adapter
    prefix_len = 20
    read = TEST_ADAPTER[:prefix_len]

    trimmed_seq, L, is_trimmed = trimmer.trim_read(read)
    assert L == prefix_len
    assert is_trimmed is True
    assert trimmed_seq == ""  # Entire read is trimmed away (remaining length 0)


def test_5_partial_adapter_suffix(trimmer):
    """Test 5: Biological sequence followed by an adapter fragment as a suffix."""
    adapter_part = TEST_ADAPPED_PART = TEST_ADAPTER[:10]  # First 10 chars: "TGGAATTCTC"
    read = "ACGT" + adapter_part

    trimmed_seq, L, is_trimmed = trimmer.trim_read(read)
    assert L == 10
    assert is_trimmed is True
    assert trimmed_seq == "ACGT"


def test_6_internal_adapter_no_suffix(trimmer):
    """Test 6: Adapter sequence present in the middle, but NOT as a suffix, yields 0 overlap."""
    # Contains the adapter in the middle, but ends with 'AAAA'
    read = "ACGT" + TEST_ADAPTER[:15] + "AAAA"

    trimmed_seq, L, is_trimmed = trimmer.trim_read(read)
    assert L == 0
    assert is_trimmed is False
    assert trimmed_seq == read


# ==========================================
# TODO: Implement Test Case 7 yourself!
# ==========================================
def test_7_equivalence_check_naive(trimmer):
    """
    Test 7: Compare KMP results against naive string slicing across
    5,000 synthetic random sequences to guarantee 100% equivalence.
    """
    # Hint: Generate random biological reads (e.g., using random.choices('ACGT', k=...)),
    # append various lengths of the adapter as suffixes, and assert that your
    # KMP-based L matches the expected suffix length found by naive iteration or slicing.
