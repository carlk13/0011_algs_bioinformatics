from collections import Counter
from collections.abc import Generator
from dataclasses import dataclass, field

from preprocessor.exact_matcher import ExactAdapterTrimmer


@dataclass
class SampleStatistics:
    barcode: str
    total_reads: int = 0
    remaining_length_distribution: Counter[int] = field(default_factory=Counter)
    sequence_frequencies: Counter[str] = field(default_factory=Counter)


class BarcodeIdentifier:
    """
    Discovers candidate barcodes from read prefixes using frequency spectrum
    analysis and pairwise Hamming distance separation.
    """

    def __init__(self, barcode_length: int = 4, min_distance: int = 3):
        self.barcode_length = barcode_length
        self.min_distance = min_distance

    def discover_barcodes(
        self,
        read_stream: Generator[str],
        sample_limit: int | None = 100000,
        min_frequency_ratio: float = 0.05,
    ) -> list[str]:
        """
        Extracts prefixes of length B, ranks by frequency, and selects the
        dominant barcodes separated by at least min_distance Hamming distance.
        """
        prefix_counts: Counter[str] = Counter()
        total = 0
        for count, r in enumerate(read_stream):
            if sample_limit and count >= sample_limit:
                break
            if len(r) >= self.barcode_length:
                prefix_counts[r[: self.barcode_length]] += 1
                total += 1

        if not prefix_counts or total == 0:
            return []

        candidates: list[str] = []
        for prefix, count in prefix_counts.most_common():
            if count / total < min_frequency_ratio:
                break
            if all(
                sum(a != b for a, b in zip(prefix, existing)) >= self.min_distance
                for existing in candidates
            ):
                candidates.append(prefix)

        return candidates


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
        min_overlap: int = 1,
    ):
        self.barcodes = known_barcodes
        self.barcode_len = len(known_barcodes[0]) if known_barcodes else 0
        self.max_mismatches = max_mismatches
        self.trimmer = ExactAdapterTrimmer(
            adapter_sequence=adapter_sequence, min_overlap=min_overlap
        )
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
        trimmed_seq, L, is_trimmed = self.trimmer.trim_read(stripped_read)
        if not is_trimmed:
            trimmed_seq = stripped_read
            L = 0

        stats = self.samples[matched_bc]
        stats.total_reads += 1
        stats.remaining_length_distribution[len(trimmed_seq)] += 1
        stats.sequence_frequencies[trimmed_seq] += 1

        return matched_bc, trimmed_seq, L
