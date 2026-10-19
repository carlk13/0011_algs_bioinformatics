from collections import Counter
from dataclasses import dataclass, field

from preprocessor.exact_matcher import ExactAdapterTrimmer


@dataclass
class SampleStatistics:
    barcode: str
    total_reads: int = 0
    remaining_length_distribution: Counter = field(default_factory=Counter)
    sequence_frequencies: Counter = field(default_factory=Counter)


class LibraryDemultiplexer:
    """
    De-multiplexes pooled sequencing libraries based on 5' barcode tags,
    trims 3' adapter fragments, and extracts per-sample length and abundance profiles.
    """

    def __init__(
        self,
        known_barcodes: list[str],
        adapter_sequence: str,
        max_mismatches: int = 1,
    ):
        self.barcodes = known_barcodes
        self.barcode_len = len(known_barcodes[0])
        self.max_mismatches = max_mismatches
        self.trimmer = ExactAdapterTrimmer(adapter_sequence, min_overlap=1)
        self.samples = {b: SampleStatistics(barcode=b) for b in known_barcodes}
        self.unassigned_count = 0

    def match_barcode(self, read_prefix: str) -> str | None:
        """Matches a read prefix to a known barcode allowing <= max_mismatches."""
        best_barcode = None
        best_dist = self.max_mismatches + 1

        for bc in self.barcodes:
            dist = sum(1 for a, b in zip(read_prefix, bc) if a != b)
            if dist < best_dist:
                best_dist = dist
                best_barcode = bc

        if best_dist <= self.max_mismatches:
            return best_barcode
        return None

    def process_read(self, read: str) -> tuple[str | None, str, int]:
        prefix = read[: self.barcode_len]
        matched_bc = self.match_barcode(prefix)

        if not matched_bc:
            self.unassigned_count += 1
            return None, read, 0

        # Remove barcode prefix
        stripped_read = read[self.barcode_len :]
        # Trim adapter suffix
        trimmed_seq, L, _ = self.trimmer.trim_read(stripped_read)

        stats = self.samples[matched_bc]
        stats.total_reads += 1
        stats.remaining_length_distribution[len(trimmed_seq)] += 1
        stats.sequence_frequencies[trimmed_seq] += 1

        return matched_bc, trimmed_seq, L
