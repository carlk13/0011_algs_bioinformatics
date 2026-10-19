import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parent.parent / "src"))

from preprocessor.approx_matcher import HammingAdapterTrimmer, SemiGlobalDPTrimmer
from preprocessor.error_analyzer import AlignmentRecord, ErrorProfiler

TEST_ADAPTER = "TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACACAGTGATCTCGTATGCCGTCTTCTGCTTG"


# ---------------------------------------------------------------------------
# 1. Exact match through approx matcher — returns 0 mismatches
# ---------------------------------------------------------------------------
def test_hamming_exact_match():
    """Exact adapter prefix should match with 0 mismatches."""
    trimmer = HammingAdapterTrimmer(TEST_ADAPTER, max_error_rate=0.10, min_overlap=1)
    prefix_len = 20
    read = "ACGTACGT" + TEST_ADAPTER[:prefix_len]
    trimmed, L, mismatches, is_trimmed = trimmer.trim_read(read)
    assert is_trimmed is True
    assert L == prefix_len
    assert mismatches == 0
    assert trimmed == "ACGTACGT"


# ---------------------------------------------------------------------------
# 2. 10% mismatch boundary — 1 mismatch in 10 bp accepted
# ---------------------------------------------------------------------------
def test_hamming_10pct_threshold_acceptance():
    """1 mismatch in 10 bp is 10% (allowed under 0.10 error rate)."""
    trimmer = HammingAdapterTrimmer(TEST_ADAPTER, max_error_rate=0.10, min_overlap=1)
    # Adapter prefix 10: "TGGAATTCTC". Mutate index 7: 'T' -> 'A' -> "TGGAATACTC"
    mutated_prefix = TEST_ADAPTER[:7] + "A" + TEST_ADAPTER[8:10]
    read = "NNNN" + mutated_prefix
    trimmed, L, mismatches, is_trimmed = trimmer.trim_read(read)
    assert is_trimmed is True
    assert L == 10
    assert mismatches == 1
    assert trimmed == "NNNN"


# ---------------------------------------------------------------------------
# 3. 10% rejection / 25% acceptance — 2 mismatches in 10 bp
# ---------------------------------------------------------------------------
def test_hamming_10pct_threshold_rejection():
    """2 mismatches in 10 bp is 20% (rejected under 0.10, accepted under 0.25)."""
    trimmer_10 = HammingAdapterTrimmer(TEST_ADAPTER, max_error_rate=0.10, min_overlap=1)
    trimmer_25 = HammingAdapterTrimmer(TEST_ADAPTER, max_error_rate=0.25, min_overlap=1)
    # Mutate 2 positions in first 10 bp:
    mutated_prefix = TEST_ADAPTER[:6] + "AA" + TEST_ADAPTER[8:10]
    read = "NNNN" + mutated_prefix

    # 10% trimmer should not accept L=10
    trimmed_10, L_10, mismatches_10, is_trimmed_10 = trimmer_10.trim_read(read)
    assert L_10 < 10

    # 25% trimmer should accept L=10 with 2 mismatches
    trimmed_25, L_25, mismatches_25, is_trimmed_25 = trimmer_25.trim_read(read)
    assert is_trimmed_25 is True
    assert L_25 == 10
    assert mismatches_25 == 2
    assert trimmed_25 == "NNNN"


# ---------------------------------------------------------------------------
# 4. Longest match priority — greedy longest-first strategy
# ---------------------------------------------------------------------------
def test_hamming_longest_match_priority():
    """If a read matches at L=30 with 2 mismatches (6.7%) the trimmer must
    select L=30 rather than falling through to a shorter perfect match."""
    trimmer = HammingAdapterTrimmer(TEST_ADAPTER, max_error_rate=0.10, min_overlap=1)

    # Build a read whose last 30 characters align to adapter[:30] with 2 mismatches.
    # Adapter[:30] = "TGGAATTCTCGGGTGCCAAGGAACTCCAGT"
    adapter_prefix_30 = list(TEST_ADAPTER[:30])
    # Mutate positions 5 and 15 (well within 10% of 30 = 3 allowed)
    adapter_prefix_30[5] = "C" if adapter_prefix_30[5] != "C" else "A"
    adapter_prefix_30[15] = "T" if adapter_prefix_30[15] != "T" else "A"
    mutated_suffix = "".join(adapter_prefix_30)

    biological_insert = "GGCCAATTGG"  # 10 bp of insert
    read = biological_insert + mutated_suffix  # total 40 bp

    trimmed, L, mismatches, is_trimmed = trimmer.trim_read(read)
    assert is_trimmed is True
    assert L == 30, f"Expected longest match L=30, got L={L}"
    assert mismatches == 2
    assert trimmed == biological_insert


# ---------------------------------------------------------------------------
# 5a. Indel test — single-base DELETION in read suffix
# ---------------------------------------------------------------------------
def test_semiglobal_dp_single_deletion():
    """A read suffix that is the adapter prefix with one base deleted in the
    middle should be detected by the Semi-Global DP aligner with 1 deletion.

    We use a longer prefix (20 bp) and delete a base from position 10, creating
    a 19-character suffix that aligns to the 20-char adapter prefix with one gap.
    The long flanking context on both sides of the gap forces the DP to prefer
    the gapped alignment over a shorter clean match.
    """
    biological_insert = "ACGTACGTACGT"
    adapter_prefix_20 = TEST_ADAPTER[:20]  # "TGGAATTCTCGGGTGCCAAG"
    # Delete position 10 (the 'G') → "TGGAATTCTCGGTGCCAAG" (19 chars)
    read_suffix = adapter_prefix_20[:10] + adapter_prefix_20[11:]
    read = biological_insert + read_suffix

    trimmer = SemiGlobalDPTrimmer(TEST_ADAPTER, min_score_threshold=5)
    trimmed, L, score, subs, ins, dels, err_pos, is_trimmed = trimmer.trim_read(
        read, detailed=True
    )

    assert is_trimmed is True, "DP should detect the deletion-containing overlap"
    assert L > 0
    # The aligner should report at least 1 deletion (gap in read relative to adapter)
    assert dels >= 1, f"Expected >=1 deletion, got {dels}"


# ---------------------------------------------------------------------------
# 5b. Indel test — single-base INSERTION in read suffix
# ---------------------------------------------------------------------------
def test_semiglobal_dp_single_insertion():
    """A read suffix that is the adapter prefix with one extra base inserted should
    be detected by the Semi-Global DP aligner with 1 insertion."""
    # Adapter[:10] = "TGGAATTCTC"
    # Insert an extra 'A' after position 4 -> "TGGAAATTCTC" (11 bp -> 10 bp adapter prefix)
    biological_insert = "ACGTACGTACGT"
    adapter_prefix = TEST_ADAPTER[:10]
    read_suffix = adapter_prefix[:4] + "A" + adapter_prefix[4:]  # "TGGAAATTCTC"
    read = biological_insert + read_suffix

    trimmer = SemiGlobalDPTrimmer(TEST_ADAPTER, min_score_threshold=5)
    trimmed, L, score, subs, ins, dels, err_pos, is_trimmed = trimmer.trim_read(
        read, detailed=True
    )

    assert is_trimmed is True, "DP should detect the insertion-containing overlap"
    assert L > 0
    # The aligner should report at least 1 insertion (extra base in read)
    assert ins >= 1, f"Expected >=1 insertion, got {ins}"


# ---------------------------------------------------------------------------
# Existing: Semi-Global DP exact match + traceback
# ---------------------------------------------------------------------------
def test_semi_global_dp_match_and_traceback():
    """Semi-Global DP aligns adapter with substitutions and indels."""
    trimmer = SemiGlobalDPTrimmer(TEST_ADAPTER, min_score_threshold=8)
    # Exact 15 bp prefix
    read = "ACGT" + TEST_ADAPTER[:15]
    trimmed, L, score, subs, ins, dels, err_pos, is_trimmed = trimmer.trim_read(
        read, detailed=True
    )
    assert is_trimmed is True
    assert L == 15
    assert subs == 0
    assert ins == 0
    assert dels == 0
    assert err_pos == []
    assert trimmed == "ACGT"


# ---------------------------------------------------------------------------
# Existing: Error profiler positional tallying (single-pass)
# ---------------------------------------------------------------------------
def test_error_profiler_single_pass():
    """ErrorProfiler accurately aggregates per-nt, per-seq, and cycle metrics."""
    profiler = ErrorProfiler(max_cycles=50, min_overlap_cutoff=10)
    # Record 1: 50 bp read, overlap 20, 1 substitution at cycle 40
    r1 = AlignmentRecord(
        read_length=50,
        overlap_length=20,
        num_substitutions=1,
        error_positions=[40],
        is_trimmed=True,
    )
    # Record 2: 50 bp read, overlap 15, 0 errors
    r2 = AlignmentRecord(
        read_length=50,
        overlap_length=15,
        num_substitutions=0,
        error_positions=[],
        is_trimmed=True,
    )
    profiler.update(r1)
    profiler.update(r2)

    summary = profiler.get_summary()
    assert summary["total_adapter_reads"] == 2
    assert summary["reads_with_errors"] == 1
    assert summary["total_aligned_bases"] == 35  # 20 + 15
    assert summary["per_sequence_error_rate"] == 0.5
    assert summary["per_nucleotide_error_rate"] == round(1 / 35, 6)

    # Check cycle 40 recorded 1 error and 2 observations
    cycle_40 = next(c for c in summary["cycle_distribution"] if c["cycle"] == 40)
    assert cycle_40["errors"] == 1
    assert cycle_40["observations"] == 2
