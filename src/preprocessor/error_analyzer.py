from collections.abc import Iterable
from dataclasses import dataclass, field


@dataclass
class AlignmentRecord:
    """Standardized record returned by approximate / semi-global matchers."""

    read_length: int
    overlap_length: int
    num_substitutions: int
    num_insertions: int = 0
    num_deletions: int = 0
    # Error positions mapped to read cycle coordinates (0 to read_length - 1)
    error_positions: list[int] = field(default_factory=list)
    is_trimmed: bool = True

    @property
    def total_errors(self) -> int:
        return self.num_substitutions + self.num_insertions + self.num_deletions


class ErrorProfiler:
    """
    Profiles sequencing error rates per nucleotide, per sequence,
    per cycle position, and error type breakdowns for Task 3.
    """

    def __init__(self, max_cycles: int = 50, min_overlap_cutoff: int = 10):
        self.max_cycles = max_cycles
        self.min_overlap_cutoff = min_overlap_cutoff

        # Single-pass accumulators
        self.total_reads_with_adapter = 0
        self.reads_with_errors = 0
        self.total_aligned_adapter_bases = 0
        self.total_substitutions = 0
        self.total_insertions = 0
        self.total_deletions = 0

        # Cycle trackers (native lists for speed)
        self.errors_at_cycle = [0] * max_cycles
        self.observations_at_cycle = [0] * max_cycles

    def update(self, record: AlignmentRecord) -> None:
        """Accumulates metrics for a single alignment in O(1) time."""
        if not record.is_trimmed or record.overlap_length < self.min_overlap_cutoff:
            return

        self.total_reads_with_adapter += 1
        L = record.overlap_length
        m = record.read_length
        errors = record.total_errors

        self.total_aligned_adapter_bases += L
        self.total_substitutions += record.num_substitutions
        self.total_insertions += record.num_insertions
        self.total_deletions += record.num_deletions

        if errors > 0:
            self.reads_with_errors += 1

        # Track observations along read cycles covered by the adapter suffix [m - L, m - 1]
        adapter_start_cycle = max(0, m - L)
        adapter_end_cycle = min(m, self.max_cycles)
        for cycle in range(adapter_start_cycle, adapter_end_cycle):
            self.observations_at_cycle[cycle] += 1

        # Track error positions (in read cycle coordinates)
        for pos in record.error_positions:
            if 0 <= pos < self.max_cycles:
                self.errors_at_cycle[pos] += 1

    def process_stream(self, records: Iterable[AlignmentRecord]) -> dict:
        """Processes an entire generator or list of alignments in a single pass."""
        for rec in records:
            self.update(rec)
        return self.get_summary()

    def get_summary(self) -> dict:
        """Computes and returns all Task 3 metrics."""
        per_nt_rate = (
            (self.total_substitutions + self.total_insertions + self.total_deletions)
            / self.total_aligned_adapter_bases
            if self.total_aligned_adapter_bases > 0
            else 0.0
        )
        per_seq_rate = (
            self.reads_with_errors / self.total_reads_with_adapter
            if self.total_reads_with_adapter > 0
            else 0.0
        )

        total_errs = (
            self.total_substitutions + self.total_insertions + self.total_deletions
        )
        type_breakdown = {
            "substitutions": {
                "count": self.total_substitutions,
                "percentage": (self.total_substitutions / total_errs * 100)
                if total_errs > 0
                else 0.0,
            },
            "insertions": {
                "count": self.total_insertions,
                "percentage": (self.total_insertions / total_errs * 100)
                if total_errs > 0
                else 0.0,
            },
            "deletions": {
                "count": self.total_deletions,
                "percentage": (self.total_deletions / total_errs * 100)
                if total_errs > 0
                else 0.0,
            },
        }

        cycle_distribution = []
        for i in range(self.max_cycles):
            obs = self.observations_at_cycle[i]
            errs = self.errors_at_cycle[i]
            rate = (errs / obs) if obs > 0 else 0.0
            cycle_distribution.append(
                {
                    "cycle": i,
                    "errors": errs,
                    "observations": obs,
                    "error_rate": round(rate, 6),
                }
            )

        return {
            "total_adapter_reads": self.total_reads_with_adapter,
            "reads_with_errors": self.reads_with_errors,
            "total_aligned_bases": self.total_aligned_adapter_bases,
            "per_nucleotide_error_rate": round(per_nt_rate, 6),
            "per_sequence_error_rate": round(per_seq_rate, 6),
            "error_type_breakdown": type_breakdown,
            "cycle_distribution": cycle_distribution,
        }
