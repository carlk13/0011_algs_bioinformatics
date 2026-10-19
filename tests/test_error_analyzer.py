import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parent.parent / "src"))

import pytest

from preprocessor.error_analyzer import AlignmentRecord, ErrorProfiler

# ════════════════════════════════════════════════════════════════════════════
# AlignmentRecord unit tests
# ════════════════════════════════════════════════════════════════════════════


class TestAlignmentRecord:
    """Tests for the AlignmentRecord dataclass."""

    def test_total_errors_all_types(self):
        """total_errors should sum substitutions + insertions + deletions."""
        rec = AlignmentRecord(
            read_length=50,
            overlap_length=20,
            num_substitutions=3,
            num_insertions=1,
            num_deletions=2,
            error_positions=[30, 35, 40, 42, 45, 46],
            is_trimmed=True,
        )
        assert rec.total_errors == 6

    def test_total_errors_zero(self):
        """A perfect alignment should report 0 total errors."""
        rec = AlignmentRecord(
            read_length=50,
            overlap_length=20,
            num_substitutions=0,
            num_insertions=0,
            num_deletions=0,
            error_positions=[],
            is_trimmed=True,
        )
        assert rec.total_errors == 0

    def test_defaults(self):
        """Optional fields should default correctly."""
        rec = AlignmentRecord(
            read_length=50,
            overlap_length=10,
            num_substitutions=1,
        )
        assert rec.num_insertions == 0
        assert rec.num_deletions == 0
        assert rec.error_positions == []
        assert rec.is_trimmed is True


# ════════════════════════════════════════════════════════════════════════════
# ErrorProfiler — per-nucleotide error rate
# ════════════════════════════════════════════════════════════════════════════


class TestErrorProfilerNucleotideRate:
    """Tests for per-nucleotide error rate calculation."""

    def test_per_nucleotide_single_record(self):
        """Single record: error_rate_nt = total_errors / aligned_bases."""
        profiler = ErrorProfiler(max_cycles=50, min_overlap_cutoff=1)
        rec = AlignmentRecord(
            read_length=50,
            overlap_length=20,
            num_substitutions=2,
            num_insertions=0,
            num_deletions=0,
            error_positions=[40, 45],
            is_trimmed=True,
        )
        profiler.update(rec)
        summary = profiler.get_summary()
        assert summary["total_aligned_bases"] == 20
        assert summary["per_nucleotide_error_rate"] == round(2 / 20, 6)

    def test_per_nucleotide_mixed_error_types(self):
        """Substitutions + insertions + deletions all contribute to per-nt rate."""
        profiler = ErrorProfiler(max_cycles=50, min_overlap_cutoff=1)
        rec = AlignmentRecord(
            read_length=50,
            overlap_length=25,
            num_substitutions=1,
            num_insertions=1,
            num_deletions=1,
            error_positions=[30, 35, 40],
            is_trimmed=True,
        )
        profiler.update(rec)
        summary = profiler.get_summary()
        expected = round(3 / 25, 6)
        assert summary["per_nucleotide_error_rate"] == expected

    def test_per_nucleotide_no_errors(self):
        """Zero errors → 0.0 per-nt error rate."""
        profiler = ErrorProfiler(max_cycles=50, min_overlap_cutoff=1)
        rec = AlignmentRecord(
            read_length=50,
            overlap_length=30,
            num_substitutions=0,
            error_positions=[],
            is_trimmed=True,
        )
        profiler.update(rec)
        assert profiler.get_summary()["per_nucleotide_error_rate"] == 0.0


# ════════════════════════════════════════════════════════════════════════════
# ErrorProfiler — per-sequence error rate
# ════════════════════════════════════════════════════════════════════════════


class TestErrorProfilerSequenceRate:
    """Tests for per-sequence error rate calculation."""

    def test_per_sequence_all_clean(self):
        """All reads with 0 errors → per-sequence error rate = 0."""
        profiler = ErrorProfiler(max_cycles=50, min_overlap_cutoff=1)
        for _ in range(10):
            profiler.update(
                AlignmentRecord(
                    read_length=50,
                    overlap_length=20,
                    num_substitutions=0,
                    error_positions=[],
                    is_trimmed=True,
                )
            )
        summary = profiler.get_summary()
        assert summary["reads_with_errors"] == 0
        assert summary["per_sequence_error_rate"] == 0.0

    def test_per_sequence_all_erroneous(self):
        """All reads have ≥1 error → per-sequence error rate = 1.0."""
        profiler = ErrorProfiler(max_cycles=50, min_overlap_cutoff=1)
        for i in range(5):
            profiler.update(
                AlignmentRecord(
                    read_length=50,
                    overlap_length=15,
                    num_substitutions=1,
                    error_positions=[40 + i],
                    is_trimmed=True,
                )
            )
        summary = profiler.get_summary()
        assert summary["reads_with_errors"] == 5
        assert summary["per_sequence_error_rate"] == 1.0

    def test_per_sequence_mixed(self):
        """3 reads: 1 with errors, 2 without → per-seq rate = 1/3."""
        profiler = ErrorProfiler(max_cycles=50, min_overlap_cutoff=1)
        profiler.update(
            AlignmentRecord(
                read_length=50,
                overlap_length=20,
                num_substitutions=2,
                error_positions=[40, 45],
                is_trimmed=True,
            )
        )
        for _ in range(2):
            profiler.update(
                AlignmentRecord(
                    read_length=50,
                    overlap_length=20,
                    num_substitutions=0,
                    error_positions=[],
                    is_trimmed=True,
                )
            )
        summary = profiler.get_summary()
        assert summary["total_adapter_reads"] == 3
        assert summary["reads_with_errors"] == 1
        assert summary["per_sequence_error_rate"] == round(1 / 3, 6)


# ════════════════════════════════════════════════════════════════════════════
# ErrorProfiler — cycle / positional error distribution
# ════════════════════════════════════════════════════════════════════════════


class TestErrorProfilerCycleDistribution:
    """Tests for positional (cycle) error tracking."""

    def test_error_at_specific_cycle(self):
        """An error at cycle 40 should be recorded in errors_at_cycle[40]."""
        profiler = ErrorProfiler(max_cycles=50, min_overlap_cutoff=1)
        profiler.update(
            AlignmentRecord(
                read_length=50,
                overlap_length=20,
                num_substitutions=1,
                error_positions=[40],
                is_trimmed=True,
            )
        )
        summary = profiler.get_summary()
        cycle_40 = summary["cycle_distribution"][40]
        assert cycle_40["errors"] == 1
        assert cycle_40["observations"] >= 1

    def test_observations_span_adapter_region(self):
        """Observations should be incremented for every cycle in [m-L, m-1]."""
        profiler = ErrorProfiler(max_cycles=50, min_overlap_cutoff=1)
        # read_length=50, overlap=10 → adapter covers cycles 40..49
        profiler.update(
            AlignmentRecord(
                read_length=50,
                overlap_length=10,
                num_substitutions=0,
                error_positions=[],
                is_trimmed=True,
            )
        )
        summary = profiler.get_summary()
        for cycle in range(40, 50):
            assert summary["cycle_distribution"][cycle]["observations"] == 1
        # Cycles before the adapter region should have 0 observations
        for cycle in range(40):
            assert summary["cycle_distribution"][cycle]["observations"] == 0

    def test_multiple_reads_cycle_accumulation(self):
        """Two reads overlapping different cycle ranges accumulate correctly."""
        profiler = ErrorProfiler(max_cycles=50, min_overlap_cutoff=1)
        # Read 1: overlap=20 → cycles 30..49
        profiler.update(
            AlignmentRecord(
                read_length=50,
                overlap_length=20,
                num_substitutions=1,
                error_positions=[35],
                is_trimmed=True,
            )
        )
        # Read 2: overlap=10 → cycles 40..49
        profiler.update(
            AlignmentRecord(
                read_length=50,
                overlap_length=10,
                num_substitutions=1,
                error_positions=[45],
                is_trimmed=True,
            )
        )
        summary = profiler.get_summary()
        # Cycle 35: 1 observation (only read 1), 1 error
        c35 = summary["cycle_distribution"][35]
        assert c35["observations"] == 1
        assert c35["errors"] == 1

        # Cycle 45: 2 observations (both reads), 1 error
        c45 = summary["cycle_distribution"][45]
        assert c45["observations"] == 2
        assert c45["errors"] == 1

    def test_cycle_error_rate_computed_correctly(self):
        """error_rate at a cycle = errors / observations."""
        profiler = ErrorProfiler(max_cycles=50, min_overlap_cutoff=1)
        # 4 reads each covering cycle 45, only 1 has an error there
        for i in range(4):
            err_pos = [45] if i == 0 else []
            subs = 1 if i == 0 else 0
            profiler.update(
                AlignmentRecord(
                    read_length=50,
                    overlap_length=10,
                    num_substitutions=subs,
                    error_positions=err_pos,
                    is_trimmed=True,
                )
            )
        c45 = profiler.get_summary()["cycle_distribution"][45]
        assert c45["observations"] == 4
        assert c45["errors"] == 1
        assert c45["error_rate"] == round(1 / 4, 6)


# ════════════════════════════════════════════════════════════════════════════
# ErrorProfiler — error type breakdown
# ════════════════════════════════════════════════════════════════════════════


class TestErrorProfilerTypeBreakdown:
    """Tests for substitution / insertion / deletion percentage breakdown."""

    def test_substitution_only(self):
        """100% substitutions when no indels are present."""
        profiler = ErrorProfiler(max_cycles=50, min_overlap_cutoff=1)
        profiler.update(
            AlignmentRecord(
                read_length=50,
                overlap_length=20,
                num_substitutions=5,
                num_insertions=0,
                num_deletions=0,
                error_positions=[30, 35, 40, 42, 45],
                is_trimmed=True,
            )
        )
        breakdown = profiler.get_summary()["error_type_breakdown"]
        assert breakdown["substitutions"]["count"] == 5
        assert breakdown["substitutions"]["percentage"] == 100.0
        assert breakdown["insertions"]["count"] == 0
        assert breakdown["deletions"]["count"] == 0

    def test_mixed_error_types_percentages(self):
        """Percentages should sum to 100% for mixed error types."""
        profiler = ErrorProfiler(max_cycles=50, min_overlap_cutoff=1)
        profiler.update(
            AlignmentRecord(
                read_length=50,
                overlap_length=25,
                num_substitutions=6,
                num_insertions=3,
                num_deletions=1,
                error_positions=list(range(30, 40)),
                is_trimmed=True,
            )
        )
        breakdown = profiler.get_summary()["error_type_breakdown"]
        total = 6 + 3 + 1
        assert breakdown["substitutions"]["percentage"] == pytest.approx(
            6 / total * 100
        )
        assert breakdown["insertions"]["percentage"] == pytest.approx(3 / total * 100)
        assert breakdown["deletions"]["percentage"] == pytest.approx(1 / total * 100)
        # Sum to 100%
        pct_sum = sum(v["percentage"] for v in breakdown.values())
        assert pct_sum == pytest.approx(100.0)

    def test_no_errors_breakdown(self):
        """All percentages should be 0 when there are no errors."""
        profiler = ErrorProfiler(max_cycles=50, min_overlap_cutoff=1)
        profiler.update(
            AlignmentRecord(
                read_length=50,
                overlap_length=20,
                num_substitutions=0,
                error_positions=[],
                is_trimmed=True,
            )
        )
        breakdown = profiler.get_summary()["error_type_breakdown"]
        for error_type in ("substitutions", "insertions", "deletions"):
            assert breakdown[error_type]["count"] == 0
            assert breakdown[error_type]["percentage"] == 0.0


# ════════════════════════════════════════════════════════════════════════════
# ErrorProfiler — min_overlap_cutoff filtering
# ════════════════════════════════════════════════════════════════════════════


class TestErrorProfilerOverlapCutoff:
    """Records with overlap < min_overlap_cutoff should be ignored."""

    def test_below_cutoff_ignored(self):
        """A record with overlap_length < min_overlap_cutoff is skipped entirely."""
        profiler = ErrorProfiler(max_cycles=50, min_overlap_cutoff=10)
        profiler.update(
            AlignmentRecord(
                read_length=50,
                overlap_length=5,  # below cutoff of 10
                num_substitutions=3,
                error_positions=[45, 46, 47],
                is_trimmed=True,
            )
        )
        summary = profiler.get_summary()
        assert summary["total_adapter_reads"] == 0
        assert summary["total_aligned_bases"] == 0

    def test_at_cutoff_accepted(self):
        """A record with overlap_length == min_overlap_cutoff is accepted."""
        profiler = ErrorProfiler(max_cycles=50, min_overlap_cutoff=10)
        profiler.update(
            AlignmentRecord(
                read_length=50,
                overlap_length=10,  # exactly at cutoff
                num_substitutions=1,
                error_positions=[42],
                is_trimmed=True,
            )
        )
        summary = profiler.get_summary()
        assert summary["total_adapter_reads"] == 1
        assert summary["total_aligned_bases"] == 10

    def test_untrimmed_record_ignored(self):
        """A record with is_trimmed=False should be skipped."""
        profiler = ErrorProfiler(max_cycles=50, min_overlap_cutoff=1)
        profiler.update(
            AlignmentRecord(
                read_length=50,
                overlap_length=20,
                num_substitutions=2,
                error_positions=[40, 45],
                is_trimmed=False,
            )
        )
        summary = profiler.get_summary()
        assert summary["total_adapter_reads"] == 0


# ════════════════════════════════════════════════════════════════════════════
# ErrorProfiler — process_stream (batch ingestion)
# ════════════════════════════════════════════════════════════════════════════


class TestErrorProfilerProcessStream:
    """Tests for the process_stream convenience method."""

    def test_process_stream_returns_summary(self):
        """process_stream should ingest all records and return the summary dict."""
        records = [
            AlignmentRecord(
                read_length=50,
                overlap_length=20,
                num_substitutions=1,
                error_positions=[40],
                is_trimmed=True,
            ),
            AlignmentRecord(
                read_length=50,
                overlap_length=15,
                num_substitutions=0,
                error_positions=[],
                is_trimmed=True,
            ),
            AlignmentRecord(
                read_length=50,
                overlap_length=25,
                num_substitutions=2,
                num_insertions=1,
                error_positions=[30, 35, 40],
                is_trimmed=True,
            ),
        ]
        profiler = ErrorProfiler(max_cycles=50, min_overlap_cutoff=1)
        summary = profiler.process_stream(iter(records))

        assert summary["total_adapter_reads"] == 3
        assert summary["reads_with_errors"] == 2
        assert summary["total_aligned_bases"] == 20 + 15 + 25
        total_errs = 1 + 0 + 3
        expected_nt_rate = round(total_errs / 60, 6)
        assert summary["per_nucleotide_error_rate"] == expected_nt_rate

    def test_process_stream_empty(self):
        """An empty stream should produce zeroed-out summary."""
        profiler = ErrorProfiler(max_cycles=50, min_overlap_cutoff=1)
        summary = profiler.process_stream(iter([]))
        assert summary["total_adapter_reads"] == 0
        assert summary["per_nucleotide_error_rate"] == 0.0
        assert summary["per_sequence_error_rate"] == 0.0


# ════════════════════════════════════════════════════════════════════════════
# ErrorProfiler — edge-case: error position at boundary
# ════════════════════════════════════════════════════════════════════════════


class TestErrorProfilerEdgeCases:
    """Edge cases for the cycle tracking arrays."""

    def test_error_at_cycle_49(self):
        """An error at the very last cycle (49) should be recorded."""
        profiler = ErrorProfiler(max_cycles=50, min_overlap_cutoff=1)
        profiler.update(
            AlignmentRecord(
                read_length=50,
                overlap_length=5,
                num_substitutions=1,
                error_positions=[49],
                is_trimmed=True,
            )
        )
        c49 = profiler.get_summary()["cycle_distribution"][49]
        assert c49["errors"] == 1
        assert c49["observations"] == 1

    def test_error_at_cycle_0(self):
        """An error at cycle 0 (the very first cycle) should be recorded."""
        profiler = ErrorProfiler(max_cycles=50, min_overlap_cutoff=1)
        profiler.update(
            AlignmentRecord(
                read_length=50,
                overlap_length=50,
                num_substitutions=1,
                error_positions=[0],
                is_trimmed=True,
            )
        )
        c0 = profiler.get_summary()["cycle_distribution"][0]
        assert c0["errors"] == 1
        assert c0["observations"] == 1

    def test_error_position_beyond_max_cycles_ignored(self):
        """Error positions >= max_cycles should be silently ignored."""
        profiler = ErrorProfiler(max_cycles=50, min_overlap_cutoff=1)
        profiler.update(
            AlignmentRecord(
                read_length=100,
                overlap_length=20,
                num_substitutions=1,
                error_positions=[55],  # beyond max_cycles=50
                is_trimmed=True,
            )
        )
        # No cycle should record this error
        summary = profiler.get_summary()
        total_cycle_errors = sum(c["errors"] for c in summary["cycle_distribution"])
        assert total_cycle_errors == 0
        # But the per-nt rate should still reflect the error
        assert summary["per_nucleotide_error_rate"] == round(1 / 20, 6)
